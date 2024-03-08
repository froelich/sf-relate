package relativeMatch

import "C"
import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"relativeMatch/internal"
	"relativeMatch/lib"
	"strconv"
	"sync"
	"time"

	"github.com/BurntSushi/toml"
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/gwas"
	"github.com/hhcho/sfgwas/mpc"
	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

var boot = "_boot"
var dec = "_dec"

type ProtocolInfo struct {
	basicProt         *BasicProtocolInfo
	simpleData        [][]float64 // samples - snps
	rowIndex          []float64
	comparisonMap     map[string][]int
	bootMap           map[string]int
	approxInt         crypto.IntervalApprox
	scaleDown         float64
	scaledownLocal    float64
	numThreads        int
	threshValue       []float64
	discretizedThresh []float64
	bucketSize        int
	blinding          bool
	single            bool // row packing, used for small
	reveal            int  // reveal = -1 = reveal all; reveal = 1 = reveal degree only; OLD: 0: reveal presence only, 1: reveal results, 2: reveal positive results only
	TestSignTest      int
	Debug_ST_flag     bool
	M                 int // number of SNPs

	simpleDataPath      string
	startingIndex       int
	numberOfColumns     int
	numberOfColumnsTest int
	npz                 bool
	n                   int
	PARA                int
	totalNbrOfRows      int
	totalNbrOfRowsTest  int
	blockLimit          int
	startKey            int
	endKey              int
	rowIndexFile        string
	columnIndexFile     string
	outFolder           string
	batchLength         int
	queryLength         int

	localTest bool
	useMPC    bool
	net       int
}

type ConfigRelative struct {
	SimpleDataPath      string           `toml:"geno_dir"`
	Separator           string           `toml:"separator"`
	NumberOfColumns     int              `toml:"number_of_columns"`
	NumberOfColumnsTest int              `toml:"number_of_columns_test"`
	StartingIndex       int              `toml:"starting_index"`
	NumberOfThreads     int              `toml:"nbr_threads"`
	ThreshValue         []float64        `toml:"thresh_value"`
	DiscretizedThresh   []float64        `toml:"discretized_thresh"`
	BucketSize          int              `toml:"bucket_size"`
	ComparisonMap       map[string][]int `toml:"comparison_map"`
	BootMap             map[string]int   `toml:"boot_map"`
	ScaleDown           float64          `toml:"scale_down"`
	TestSignTest        int              `toml:"test_sign_test"`
	// approx
	Degree     int       `toml:"degree"`
	Iter       int       `toml:"iter"`
	A          float64   `toml:"A"`
	B          float64   `toml:"B"`
	Blinding   bool      `toml:"Blinding"`
	Single     bool      `toml:"single"`
	Reveal     int       `toml:"reveal"`
	ScaleDowns []float64 `toml:"scale_down_discretize"`
	// Thrs_discretize []float64 `toml:"thrs_discretize"`
	Debug_ST_flag bool `toml:"debug_st_flag"`

	Npz                bool   `toml:"npz"`
	TotalNbrOfRows     int    `toml:"total_nbr_rows"`
	TotalNbrOfRowsTest int    `toml:"total_nbr_rows_test"`
	BlockLimit         int    `toml:"block_limit"`
	StartKey           int    `toml:"start_key"`
	EndKey             int    `toml:"end_key"`
	RowIndexFile       string `toml:"hash_table_dir"`
	ColumnIndexFile    string `toml:"sketched_snps_dir"`
	BatchLength        int    `toml:"batch_length"`
	QueryLength        int    `toml:"query_length"`
	LocalTest          bool   `toml:"local_test"`
	Haps_dir           string `toml:"haps_dir"`
	PARA               int    `toml:"PARA"`

	UseMPC bool `toml:"use_mpc"`
}

func InitializeRelativeMatchingProtocol(pid int, configFolder string, network mpc.ParallelNetworks) (relativeProt *ProtocolInfo) {
	basicProt := InitializeBasicProtocol(pid, configFolder, network)

	// Import global parameters
	configRelativeGlobal := new(ConfigRelative)
	if _, err := toml.DecodeFile(filepath.Join(configFolder, "configGlobal.toml"), configRelativeGlobal); err != nil {
		fmt.Println(err)
		return
	}
	// Import local parameters
	configRelativeLocal := new(ConfigRelative)
	if _, err := toml.DecodeFile(filepath.Join(configFolder, fmt.Sprintf("configLocal.Party%d.toml", pid)), configRelativeLocal); err != nil {
		fmt.Println(err)
		return nil
	}
	// define the local n and M here
	// read it from the file configRelativeLocal.haps_dir/data_dim.txt
	fname := filepath.Join(configRelativeLocal.Haps_dir, "data_dim.txt")
	// read the three integers, one per line
	n, M, _ := internal.ReadDataDim(fname)

	return &ProtocolInfo{
		basicProt:           basicProt,
		comparisonMap:       configRelativeGlobal.ComparisonMap,
		bootMap:             configRelativeGlobal.BootMap,
		approxInt:           crypto.IntervalApprox{A: configRelativeGlobal.A, B: configRelativeGlobal.B, Degree: configRelativeGlobal.Degree, Iter: configRelativeGlobal.Iter, InverseNew: false, ScaleDowns: configRelativeGlobal.ScaleDowns},
		TestSignTest:        configRelativeGlobal.TestSignTest,
		Debug_ST_flag:       configRelativeGlobal.Debug_ST_flag,
		scaleDown:           configRelativeGlobal.ScaleDown,
		numThreads:          configRelativeGlobal.NumberOfThreads,
		threshValue:         configRelativeGlobal.ThreshValue,
		discretizedThresh:   configRelativeGlobal.DiscretizedThresh,
		M:                   M,
		bucketSize:          configRelativeGlobal.BucketSize,
		blinding:            configRelativeGlobal.Blinding,
		single:              configRelativeGlobal.Single,
		reveal:              configRelativeGlobal.Reveal,
		simpleDataPath:      configRelativeLocal.SimpleDataPath + "all_chrs.bin",
		startingIndex:       configRelativeLocal.StartingIndex,
		PARA:                configRelativeLocal.PARA,
		numberOfColumns:     configRelativeGlobal.NumberOfColumns,
		numberOfColumnsTest: configRelativeGlobal.NumberOfColumnsTest,
		npz:                 configRelativeGlobal.Npz,
		n:                   n,
		totalNbrOfRows:      configRelativeGlobal.TotalNbrOfRows,
		totalNbrOfRowsTest:  configRelativeGlobal.TotalNbrOfRowsTest,
		blockLimit:          configRelativeGlobal.BlockLimit,
		startKey:            configRelativeGlobal.StartKey,
		endKey:              configRelativeGlobal.EndKey,
		outFolder:           configFolder + "/out/",
		rowIndexFile:        configRelativeLocal.RowIndexFile + "ID_table.npz",
		columnIndexFile:     configRelativeLocal.ColumnIndexFile + "SNPs.npz",
		batchLength:         configRelativeGlobal.BatchLength,
		queryLength:         configRelativeGlobal.QueryLength,
		localTest:           configRelativeGlobal.LocalTest,
		useMPC:              configRelativeGlobal.UseMPC,
	}
}

type BatchedCVec map[int]crypto.CipherVector

type KinshipResult struct {
	Result crypto.CipherVector
}

type GlobalPhase1Result struct {
	Result []BatchedCVec // each element a maps from BatchID to a vector of ciphertext that encrypts the result (comparison values/kinships etc)
}

type LocalPhase1Result struct {
	Result []crypto.CipherVector
}

type ComparisonResultMPC struct {
	Result           []mpc_core.RVec
	Index            mpc_core.RVec
	OtherIndex       mpc_core.RVec
	ResultControlled []mpc_core.RVec
	NbrRepeat        int
}

// ComparisonDataOther contains data to be sent to other party for comparison
type ComparisonDataOther struct {
	YClear   *mat.Dense
	Y        crypto.CipherMatrix
	YSquare  crypto.CipherVector
	Yhet     crypto.CipherVector
	YhetInv  crypto.CipherVector
	Yindex   crypto.CipherVector
	Index    []int
	NbrRows  int
	Repeated int
}

// ComparisonDataLocal contains local data prepared for comparison
type ComparisonDataLocal struct {
	X       *mat.Dense
	XSquare []float64
	Xhet    []float64
	XhetInv []float64
	Index   []float64
}

// alwaysBootstrapping nodes are always ready to bootstrap to help other parties, dedicated thread
func (pi *ProtocolInfo) alwaysBootstrapping() {
	for p := 1; p < pi.basicProt.Config.NumMainParties+1; p++ {
		if p != pi.basicProt.MpcObj[0].GetPid() {
			go func(p int) {
				for true {
					log.LLvl1("Ready to Bootstrap on ", strconv.Itoa(p)+boot)
					err := pi.basicProt.MpcObj[pi.bootMap[strconv.Itoa(p)+boot]].Network.CollectiveBootstrapWithErr(pi.basicProt.Cps, nil, p)
					if err != nil {
						break
					}
				}

			}(p)
		}
	}
}

// alwaysDecrypting nodes are always ready to decrypt to help other parties, dedicated thread
func (pi *ProtocolInfo) alwaysDecrypting() {
	for p := 1; p < pi.basicProt.Config.NumMainParties+1; p++ {
		if p != pi.basicProt.MpcObj[0].GetPid() {
			go func(p int) {
				for true {
					log.LLvl1("Ready to decrypt on ", strconv.Itoa(p)+dec)
					_, err := pi.basicProt.MpcObj[pi.bootMap[strconv.Itoa(p)+dec]].Network.CollectiveDecryptVecWithErr(pi.basicProt.Cps, nil, p)
					if err != nil {
						break
					}
				}

			}(p)
		}

	}
}

func (pi *ProtocolInfo) OpenInputData() (*os.File, []int32, []int32, [][]float64) {
	var rowIndexVector, colIndexVector []int32
	var matrixFile *os.File
	var err error

	// read number of columns from a file
	rowIndexVector, colIndexVector = internal.ReadNPZVectors(pi.rowIndexFile, pi.columnIndexFile)

	matrixFile, err = os.Open(pi.simpleDataPath)
	if err != nil {
		log.Error(err)
	}
	return matrixFile, rowIndexVector, colIndexVector, nil
}

func (pi *ProtocolInfo) readbatchOfInput(batchStart, exec int, matrixFile *os.File, rowIndexVector, colIndexVector []int32, dataset [][]float64, startKeyGlobal int) ([][]float64, []float64) {
	var batch [][]float64
	batchIndex := make([]float64, pi.batchLength)

	if pi.basicProt.MpcObj[0].GetPid() > 0 {
		// might need to fix this
		endKey := batchStart + pi.batchLength
		if endKey > startKeyGlobal+pi.totalNbrOfRowsTest {
			endKey = startKeyGlobal + pi.totalNbrOfRowsTest
		}
		log.LLvl1(pi.basicProt.MpcObj[0].GetPid(), "READING from row: ", batchStart, " to ", endKey-1)

		if pi.npz {
			readBatchMatrix := make([][]float64, pi.batchLength)
			for batchElem := 0; batchElem < pi.batchLength; batchElem++ {
				rowIndex := int32(-1)
				if batchStart+batchElem < endKey {
					rowIndex = rowIndexVector[batchStart+batchElem]
					if rowIndex == -1 {
						rowIndex = 1
					}
					batchIndex[batchElem] = float64(uint64(rowIndex)) /// 10000.0
				}
				// read row from matrix
				readBatchMatrix[batchElem] = make([]float64, len(colIndexVector))
				readBatchMatrix[batchElem] = ReadRowNpz(int(rowIndex), matrixFile, pi.totalNbrOfRows, pi.numberOfColumns, colIndexVector)

			}
			log.LLvl1("BATCH DIMENSIONS: ", len(readBatchMatrix), len(readBatchMatrix[0]))
			batch = readBatchMatrix
		} else {
			// this is not used in UKB/AoU experiments
			panic("Not Should read from npz for row/col indices!!")
		}
	}
	return batch, batchIndex
}

func ReadRowNpz(rowIndex int, matrixFile *os.File, totalNbrOfRows, numberOfColumns int, colIndexVector []int32) []float64 {
	var row []float64
	if rowIndex != -1 {
		row = internal.ReadFilteredRows(matrixFile, int(rowIndex), totalNbrOfRows, numberOfColumns, len(colIndexVector), colIndexVector)
	} else { // otherwise set to one of the rows, avoid all zeros
		row = internal.ReadFilteredRows(matrixFile, 1, totalNbrOfRows, numberOfColumns, len(colIndexVector), colIndexVector)
	}
	return row
}

// BatchProtocols executes the protocol iteratively on batches of the inputs, and each batch has max of 8192 comparisons (in CKKS)
func (pi *ProtocolInfo) BatchProtocols(configFolder string, sending bool, startKeyGlobal int, mutex *sync.Mutex) LocalPhase1Result {
	pid := pi.basicProt.MpcObj[0].GetPid()
	if pid > 0 {
		if !pi.useMPC {
			pi.alwaysBootstrapping()
			pi.alwaysDecrypting()
		}
	}

	var rowIndexVector, colIndexVector []int32
	var dataset [][]float64
	var matrixFile *os.File
	if pid > 0 {
		matrixFile, rowIndexVector, colIndexVector, dataset = pi.OpenInputData()
		// read in number of columns here
		if colIndexVector != nil {
			colIndexVector = colIndexVector[:pi.numberOfColumnsTest]
		}

	}

	var batchResult LocalPhase1Result
	numArrays := []int{1, len(pi.threshValue), len(pi.discretizedThresh), 0}[pi.reveal]
	batchResult.Result = make([]crypto.CipherVector, numArrays)
	exec := 0
	// go through batches of input
	for batchStart := pi.startKey; batchStart < pi.endKey; batchStart = batchStart + pi.batchLength {
		time_start := time.Now()

		pi.simpleData, pi.rowIndex = pi.readbatchOfInput(batchStart, exec, matrixFile, rowIndexVector, colIndexVector, dataset, startKeyGlobal)
		exec++

		currentResultHE, _ := pi.RelativeSearchProtocol(sending, batchStart, mutex)

		if pi.reveal == 0 || pi.reveal == 1 || pi.reveal == 2 {
			if pi.basicProt.MpcObj[0].GetPid() > 0 {
				// should append all result into a single result table for this batch
				// need to figure out what pid to put here
				otherPid := 3 - pid
				rst := currentResultHE[otherPid]
				// split the vector and append to the batch result
				if numArrays != len(rst.Result[0]) {
					panic("Number of arrays does not match")
				}
				for i := 0; i < numArrays; i++ {
					batchResult.Result[i] = append(batchResult.Result[i], rst.Result[0][i])
				}
			}
		} else if pi.reveal == 3 {
			if pi.basicProt.MpcObj[0].GetPid() > 0 {
				// decrypt results --- need to acquire a lock
				mutex.Lock()
				decrypted := pi.decryptVectorForDebugging(pi.basicProt.Cps, currentResultHE[3-pid].Result[0], pid)
				mutex.Unlock()
				log.LLvl1(pi.outFolder + "raw/kinship_block_" + strconv.Itoa(batchStart) + "_party" + strconv.Itoa(pid) + ".txt")
				save_array(decrypted, pi.outFolder+"raw/kinship_block_"+strconv.Itoa(batchStart)+"_party"+strconv.Itoa(pid)+".txt", false, false)
			}
		}
		log.LLvl1("Time for batch ", batchStart, " is ", time.Since(time_start))
	}
	return batchResult
}

func (pi *ProtocolInfo) RelativeSearchProtocol(sending bool, batchStart int, mutex *sync.Mutex) (map[int]LocalPhase1Result, map[int]ComparisonResultMPC) {
	pid := pi.basicProt.MpcObj[0].GetPid()
	cps := pi.basicProt.Cps

	log.LLvl1(pid, ": ", time.Now(), "Finished Setup, starts Protocol")

	compResults := make(map[int]KinshipResult, 0)
	listPartiesSent := make([]int, 0) // filled when sending, based on comparisonMap

	var Xlocal ComparisonDataLocal // local data used to send to other parties and do local comparisons

	wg := sync.WaitGroup{} // parallelization
	var mutexGlobal sync.Mutex
	// data to receive
	var Y ComparisonDataOther
	// data to send
	var X ComparisonDataOther
	var nbrofActualRepeat int

	if pid > 0 { // all nodes that have data
		log.LLvl1(pid, "prepares her data") // prepare X and the vectors for sum(X^2) and rows
		timePrepareData := time.Now()
		if pi.reveal != 0 {
			Xlocal = prepareLocalData(pid, pi.simpleData, pi.reveal, 1.0/pi.scaledownLocal)
		} else {
			// for reveal == 0 , pre-scale the het values such that no mult is needed afterwards
			Xlocal = prepareLocalData(pid, pi.simpleData, pi.reveal, 1.0/pi.scaledownLocal*pi.threshValue[0])
		}
		Xlocal.Index = pi.rowIndex
		log.LLvl1(pid, ": prepared data, TIME:", time.Since(timePrepareData))
	}

	// send encrypted data to requested pids
	if pid == 0 && pi.useMPC { // for mpc, node 0 needs to help
		panic("MPC is not implemented")
	} else if pid > 0 && sending {
		// NOTE: this implementation only works for two parties
		otherPid := 3 - pid
		log.LLvl1(pid, "prepares data to send to ", otherPid) // key in comparisonMap = node that computes
		timeSendData := time.Now()
		otherNetwork := pi.basicProt.MpcObj.GetNetworks()[pi.bootMap[strconv.Itoa(otherPid)]]
		log.LLvl1(pid, " sends its data to ", otherPid)
		listPartiesSent = append(listPartiesSent, otherPid)
		if !pi.useMPC {
			if !pi.single {
				// if number of rows is small, we already duplicate --> fill ciphertext
				nbrofPossibleRepeat := int(float64(cps.Params.Slots()) / float64(len(pi.simpleData)))
				log.LLvl1(pid, "has ", len(pi.simpleData), " rows and can duplicate it ", nbrofPossibleRepeat, "times.")
				// Maybe should scale down local here too?
				X, nbrofActualRepeat = prepareLocalDataToSend(cps, pid, Xlocal, nbrofPossibleRepeat, 1.0/pi.scaleDown, pi.bucketSize)
			} else {
				panic("1 vs n case is not implemented")
			}

			// send
			wg.Add(1)
			go func(X ComparisonDataOther, nbrofActualRepeat, otherPid int) {
				defer wg.Done()
				pi.sendComparisonDataOther(cps, otherNetwork, X, otherPid, len(pi.simpleData), nbrofActualRepeat, pi.single,
					pi.blockLimit, pi.numThreads)
			}(X, nbrofActualRepeat, otherPid)
		} else {
			panic("MPC is not implemented")
		}

		log.LLvl1(pid, ": ", "encrypted and sent her data, TIME:", time.Since(timeSendData))
	}

	// receive data from other pids
	if !sending && pid > 0 {
		log.LLvl1("Ready to receive from ", pi.comparisonMap[strconv.Itoa(pid)])
		// wg.Add(1)
		otherPid := 3 - pid
		log.LLvl1(pid, " ready to receive data from ", otherPid)
		// if not single: here the matrix Y is not received, received latter to parallelize computation and reception
		if !pi.useMPC {
			Y = pi.receiveComparisonDataOther(cps, pi.basicProt.MpcObj.GetNetworks()[pi.bootMap[strconv.Itoa(pid)]], otherPid, pi.single)
		}

		timeCompare := time.Now()
		// compute King & comparison
		log.LLvl1(pid, " computes related indicator coefficients for ", otherPid)
		var compResultsLocalTmp KinshipResult
		ownNetwork := pi.basicProt.MpcObj.GetNetworks()[pi.bootMap[strconv.Itoa(pid)]]
		ownNetworkBoot := pi.basicProt.MpcObj[pi.bootMap[strconv.Itoa(pid)+boot]].Network

		if !pi.useMPC {
			log.LLvl1(pid, " has started all comparisons for pids ", pi.comparisonMap[strconv.Itoa(pid)])
			compResultsLocalTmp = pi.computeKinshipHE(pid, ownNetworkBoot, Xlocal, Y, ownNetwork, otherPid)
			mutexGlobal.Lock()
			compResults[otherPid] = compResultsLocalTmp
			mutexGlobal.Unlock()
		} else {
			panic("MPC is not implemented")
		}

		log.LLvl1(pid, ": time to compute coeff and compare Threshold: TIME: ", time.Since(timeCompare))
	}
	if pi.useMPC {
		panic("MPC is not implemented")
	}

	// mutex.Lock()
	// pi.basicProt.MpcObj.GetNetworks().PrintNetworkLog()
	// mutex.Unlock()

	timeExchangeResults := time.Now()

	allResultsToCombine := make(map[int]LocalPhase1Result, 0)

	if pid > 0 && !sending {
		// exchange comparison results
		// send all computed results
		otherPid := 3 - pid
		otherNetwork := pi.basicProt.MpcObj.GetNetworks()[pi.bootMap[strconv.Itoa(otherPid)]]
		log.LLvl1(pid, "sends its result to ", otherPid)
		// This function should be modified to handle
		// more potential use cases if needed
		compResultToSend := createControlledOutput(cps, compResults[otherPid], pi.reveal, pi.single)

		sendControlledOutput(compResultToSend, otherPid, otherNetwork, pi.localTest, pi.reveal)

		// save result in local map
		allResultsToCombine[otherPid] = compResultToSend
	}

	// receive all computed results and save them
	for _, v := range listPartiesSent {
		wg.Add(1)
		go func(v int) {
			defer wg.Done()

			// receive
			otherPid := v
			ownNetwork := pi.basicProt.MpcObj.GetNetworks()[pi.bootMap[strconv.Itoa(pid)]]
			allResultsToCombine[otherPid] = receiveControlledOutput(cps, otherPid, ownNetwork, pi.localTest, pi.reveal)

		}(v)
	}
	wg.Wait()

	log.LLvl1(pid, "has all results to combine ", len(allResultsToCombine), " TIME:", time.Since(timeExchangeResults))
	mutex.Lock()
	pi.basicProt.MpcObj.GetNetworks().PrintNetworkLog()
	mutex.Unlock()

	// 	final processing of all results
	if pid > 0 {
		// we can use this network to locally decrypt (for debugging)
		// ownNetworkDec := pi.basicProt.MpcObj.GetNetworks()[pi.bootMap[strconv.Itoa(pid)+dec]]
		return allResultsToCombine, nil
	}
	return nil, nil
}

func hetG(mat [][]float64, reveal int, scale float64) ([]float64, []float64, float64) {
	counthet := make([]float64, len(mat))
	hets := make([]float64, len(mat))
	counthetInv := make([]float64, len(mat))
	maxHet := 0.0
	for i := 0; i < len(mat); i++ {
		for j := 0; j < len(mat[0]); j++ {
			if mat[i][j] == 1 {
				counthet[i] = counthet[i] + 1
			}
		}
		if counthet[i] == 0 {
			counthet[i] = 1
		}
		hets[i] = counthet[i] * scale
		counthetInv[i] = 1.0 / hets[i]
		if counthet[i] > maxHet {
			maxHet = counthet[i]
		}

	}
	return hets, counthetInv, maxHet
}

func prepareLocalData(pid int, data [][]float64, reveal int, scale float64) ComparisonDataLocal {
	Xlocal := ComparisonDataLocal{}
	Xlocal.X = mat.NewDense(len(data), len(data[0]), internal.Flatten(data))
	xrows, xcols := Xlocal.X.Dims()

	xSquareMatrix := mat.DenseCopyOf(Xlocal.X)
	xSquareMatrix.Apply(square, Xlocal.X)
	maxHet := 0.0
	Xlocal.Xhet, Xlocal.XhetInv, maxHet = hetG(data, reveal, scale)
	log.LLvl1(pid, " max observed HET value", maxHet)

	Xlocal.XSquare = make([]float64, len(data))
	for i := 0; i < len(data); i++ {
		Xlocal.XSquare[i] = gwas.SumF(xSquareMatrix.RawRowView(i))
	}
	log.LLvl1(pid, " local data X: ", xrows, "x", xcols, " 2 vectors", len(Xlocal.XSquare), len(Xlocal.Xhet))

	return Xlocal
}

func prepareLocalDataToSend(cps *crypto.CryptoParams, pid int, Xlocal ComparisonDataLocal, nbrofPossibleRepeat int, scale float64, bucketSize int) (ComparisonDataOther, int) {
	X := ComparisonDataOther{}

	nbrofActualRepeat := 0 // counting the actual number of repetitions, either to fill a ciphertext or to the number of rows of other party
	XSquareScaled := gwas.ScaleF(Xlocal.XSquare, scale)

	XDenseToSend := mat.DenseCopyOf(Xlocal.X)
	XDenseToSend.Scale(scale, XDenseToSend)
	hetXClearToSend := Xlocal.Xhet

	XindexToSend := Xlocal.Index
	XindexToSendInit := Xlocal.Index
	hetXInvClearToSend := Xlocal.XhetInv
	XSquareClearToSend := XSquareScaled

	if nbrofPossibleRepeat > 1 {
		var XDenseSingle *mat.Dense
		XDenseSingle = mat.DenseCopyOf(XDenseToSend)
		rowBase, colBase := XDenseToSend.Dims()
		log.LLvl1("Bucket size is ", bucketSize, " rows, so duplicating a max of ", bucketSize, " or ", nbrofPossibleRepeat, " times")

		for r := 1; r < bucketSize && r < nbrofPossibleRepeat; r++ {
			if r%10 == 0 {
				log.LLvl1("replicating local data: ", r)
			}
			// if example with large bucket --> can further optimize
			currRow, _ := XDenseToSend.Dims()
			toReceive := mat.NewDense(currRow+rowBase, colBase, nil)
			toReceive.Stack(XDenseToSend, XDenseSingle) // stack full columns
			XDenseToSend = mat.DenseCopyOf(toReceive)

			hetXClearToSend = append(hetXClearToSend, Xlocal.Xhet...)

			hetXInvClearToSend = append(hetXInvClearToSend, Xlocal.XhetInv...)
			XindexToSend = append(XindexToSend, XindexToSendInit...)
			XSquareClearToSend = append(XSquareClearToSend, XSquareScaled...)
			nbrofActualRepeat++
		}

	}

	xOrigRow, xOrigCol := Xlocal.X.Dims()
	xSendRow, xSendCol := XDenseToSend.Dims()
	log.LLvl1("original dims: ", xOrigRow, " x ", xOrigCol, " to ", xSendRow, " x ", xSendCol)

	// encrypt data
	X.YClear = XDenseToSend
	X.YSquare, _ = crypto.EncryptFloatVectorWithScale(cps, XSquareClearToSend, cps.Params.Scale())
	X.Yhet, _ = crypto.EncryptFloatVectorWithScale(cps, hetXClearToSend, float64(cps.Params.Qi()[8]))
	X.YhetInv, _ = crypto.EncryptFloatVectorWithScale(cps, hetXInvClearToSend, float64(cps.Params.Qi()[8]))
	X.Yindex, _ = crypto.EncryptFloatVectorWithScale(cps, XindexToSend, float64(cps.Params.Qi()[8]))
	log.LLvl1(pid, " repeated its local data", nbrofActualRepeat, " times, now encrypts with scale: ", cps.Params.Scale(),
		" (X^2) and ", X.Yhet[0].Scale(), " (het)")

	return X, nbrofActualRepeat
}

func prepareLocalSingleDataToSend(cps *crypto.CryptoParams, Xlocal ComparisonDataLocal, scale float64, bucketSize int) ComparisonDataOther {
	// bucketsize in this case is the number of rows for other party
	X := ComparisonDataOther{}

	XSquareScaled := gwas.ScaleF(Xlocal.XSquare, scale)
	XDenseToSend := mat.DenseCopyOf(Xlocal.X.T()) // transpose to use default column encoding as row encoding
	_, localRows := XDenseToSend.Dims()
	XDenseToSend.Scale(scale, XDenseToSend)

	hetXClearToSend := make([]float64, 0) //Xlocal.Xhet
	hetXInvClearToSend := make([]float64, 0)
	XSquareClearToSend := make([]float64, 0) //XSquareScaled

	// Replicate the vectors to simplify other node work
	// Potential improvement: limit the duplication to a single ciphertext per local row
	log.LLvl1("Other party has ", bucketSize, " rows")
	for i := 0; i < localRows; i++ {
		for r := 0; r < bucketSize; r++ {
			hetXClearToSend = append(hetXClearToSend, Xlocal.Xhet[i])
			hetXInvClearToSend = append(hetXInvClearToSend, Xlocal.XhetInv[i])
			XSquareClearToSend = append(XSquareClearToSend, XSquareScaled[i])
		}
		if localRows > 1 { // results of matrix mult will be in different ciphertexts
			for r := bucketSize; r < cps.GetSlots(); r++ { // padding
				hetXClearToSend = append(hetXClearToSend, 1.0)
				hetXInvClearToSend = append(hetXInvClearToSend, 1.0)
				XSquareClearToSend = append(XSquareClearToSend, 1.0)
			}
		}
	}

	// encrypt the data
	X.Y = crypto.EncryptDense(cps, XDenseToSend)
	X.Y = crypto.DropLevel(cps, X.Y, 6)
	X.Y = crypto.CMatRescale(cps, X.Y)
	X.YSquare, _ = crypto.EncryptFloatVectorWithScale(cps, XSquareClearToSend, cps.Params.Scale()) //*cps.Params.Scale())
	X.Yhet, _ = crypto.EncryptFloatVectorWithScale(cps, hetXClearToSend, float64(cps.Params.Qi()[8]))
	X.YhetInv, _ = crypto.EncryptFloatVectorWithScale(cps, hetXInvClearToSend, float64(cps.Params.Qi()[8]))
	return X
}

func (pi *ProtocolInfo) sendComparisonDataOther(cps *crypto.CryptoParams, net *mpc.Network, X ComparisonDataOther, otherPid, nbOfRows,
	nbrOfActualRepeat int, single bool, blockLimit, numThreads int) {

	if single {
		net.SendInt(len(X.Y), otherPid)    // # of columns
		net.SendInt(len(X.Y[0]), otherPid) // # of ciphertexts for rows
		net.SendCipherMatrix(X.Y, otherPid)
	}
	net.SendInt(len(X.Yhet), otherPid) // # of ciphers
	net.SendCipherVector(X.YSquare, otherPid)
	net.SendInt(nbOfRows, otherPid)            // # of rows
	net.SendInt(nbrOfActualRepeat+1, otherPid) // # of repeats
	// drop the level of Yhet if reveal == 0
	if pi.reveal == 0 {
		// maybe not 8 if we switch to the larger params
		X.Yhet = crypto.DropLevelVec(cps, X.Yhet, 8)
	}
	net.SendCipherVector(X.Yhet, otherPid)
	if pi.reveal != 0 {
		// used the default mode of outputing one boolean per sample
		net.SendCipherVector(X.YhetInv, otherPid)
	}
	if !single {
		net.SendCipherVector(X.Yindex, otherPid)
		// send the matrix by blocks
		log.LLvl1("Sending the matrix by blocks")
		log.LLvl1(X.YClear.Dims())
		pi.sendMatrixByBlocks(cps, X.YClear, blockLimit, numThreads, net, otherPid)
	}

}

func (pi *ProtocolInfo) sendMatrixByBlocks(cps *crypto.CryptoParams, matToSend *mat.Dense, blockLimit, numThreads int, net *mpc.Network, otherPid int) {
	// send Y by blocks
	nbrRows, nbrCols := matToSend.Dims()
	nbrOfBlocks := int(math.Ceil(float64(nbrCols) / float64(blockLimit)))
	lastblockSize := nbrCols % blockLimit
	littleWg := sync.WaitGroup{}
	for b := 0; b < nbrOfBlocks; b++ {
		var block *mat.Dense
		if b == nbrOfBlocks-1 {
			if lastblockSize == 0 {
				lastblockSize = blockLimit
			}
			block = mat.DenseCopyOf(matToSend.Slice(0, nbrRows, nbrCols-lastblockSize, nbrCols))
		} else {
			block = mat.DenseCopyOf(matToSend.Slice(0, nbrRows, b*blockLimit, (b+1)*blockLimit))
		}
		blockRows, blockCols := block.Dims()

		verbose := (blockLimit * b % 10000) == 0
		if verbose {
			log.LLvl1("Encrypting block: ", b, ": ", blockRows, blockCols)
		}
		blockEnc := crypto.EncryptDenseParallel(cps, block, numThreads)

		blockEnc = crypto.DropLevel(cps, blockEnc, 3)
		numColsToDecrypt := blockCols
		if numColsToDecrypt > 10 {
			numColsToDecrypt = 10
		}

		sbytes, cmbytes := mpc.MarshalCM(blockEnc)
		if verbose {
			log.LLvl1("Sending block: ", b, ": ", blockRows, blockCols)
		}
		littleWg.Wait()

		net.SendInt(len(blockEnc), otherPid) // # of columns
		net.SendInt(len(blockEnc[0]), otherPid)
		littleWg.Add(1)
		go func(blockEnc crypto.CipherMatrix) {
			defer littleWg.Done()
			net.SendCipherMatrixMarshalled(sbytes, cmbytes, otherPid)
		}(blockEnc)

		if verbose {
			log.LLvl1("Sent block: ", b, ": ", blockRows, blockCols)
		}
	}
}

func (pi *ProtocolInfo) receiveComparisonDataOther(cps *crypto.CryptoParams, net *mpc.Network, otherPid int, single bool) ComparisonDataOther {
	Y := ComparisonDataOther{}
	if single {
		ncolsY := net.ReceiveInt(otherPid)
		log.LLvl1("received #columns: ", ncolsY)
		nrowsCipherY := net.ReceiveInt(otherPid)
		log.LLvl1("received #nrowsCiphers: ", nrowsCipherY)
		Y.Y = net.ReceiveCipherMatrix(cps, ncolsY, nrowsCipherY, otherPid)
	}
	nvec := net.ReceiveInt(otherPid)
	log.LLvl1("received #nvecForHetCiphers: ", nvec)
	Y.YSquare = net.ReceiveCipherVector(cps, nvec, otherPid)
	log.LLvl1("Received Vector YSquare (#ciphers=", len(Y.YSquare), ")")
	Y.NbrRows = net.ReceiveInt(otherPid)
	log.LLvl1("received other # of rows: ", Y.NbrRows)
	Y.Repeated = net.ReceiveInt(otherPid)
	log.LLvl1("received other # rpetitions: ", Y.Repeated)

	Y.Yhet = net.ReceiveCipherVector(cps, nvec, otherPid)
	log.LLvl1("Received vector het (#ciphers=", len(Y.Yhet), ")")

	if pi.reveal != 0 {
		// only accumulating to boolean output
		Y.YhetInv = net.ReceiveCipherVector(cps, nvec, otherPid)
		log.LLvl1("Received vector hetInv (#ciphers=", len(Y.YhetInv), ")")
	}
	if !single {
		Y.Yindex = net.ReceiveCipherVector(cps, nvec, otherPid)
	}

	log.LLvl1("Received vector index (#ciphers=", len(Y.YhetInv), ")")

	return Y
}

func (pi *ProtocolInfo) receiveMatrixBlocksAndAggregate(cps *crypto.CryptoParams, initial crypto.CipherVector, nbrBlocks int, netMat *mpc.Network, otherPid, numThreads int,
	X ComparisonDataLocal, xRows int, bucketSize int, nbrBuckets int, nbrRepeatY int) crypto.CipherVector {
	wg := sync.WaitGroup{}
	var mutex sync.Mutex
	// err_counter := 0
	firstXIndex := 0
	for b := 0; b < nbrBlocks; b++ {
		verbose := (b * X.X.RawMatrix().Cols % 1000) == 0
		if verbose {
			log.LLvl1("ready for block: ", b)
		}
		ncolsY := netMat.ReceiveInt(otherPid)
		nrowsCipherY := netMat.ReceiveInt(otherPid)
		if verbose {
			log.LLvl1("receiving block: ", b, ":", ncolsY, nrowsCipherY)
		}
		Ymat := netMat.ReceiveCipherMatrix(cps, ncolsY, nrowsCipherY, otherPid)
		wg.Wait()
		if verbose {
			log.LLvl1("processing block: ", b, ":", ncolsY, nrowsCipherY)
		}
		vparallelize := int(math.Ceil(float64(ncolsY) / float64(numThreads)))
		for j := 0; j < ncolsY; j = j + vparallelize {
			wg.Add(1)
			go func(iT int, ncolsY int, Ymat crypto.CipherMatrix, firstXIndex int) {
				defer wg.Done()
				distanceTmpTest := gwas.NewCipherVectorAccV2(cps, int(math.Ceil(float64(xRows*bucketSize)/float64(cps.Params.Slots()))), Ymat[0][0].Level())
				for k := 0; k < vparallelize && (k+iT < ncolsY); k++ {
					distanceTmpTest = pi.compute2XYColumn(cps, X, Ymat, bucketSize, nbrBuckets, nbrRepeatY, k+iT, firstXIndex+k+iT, distanceTmpTest, 0)
				}
				distanceTmp := gwas.ModularReduceV2(cps, distanceTmpTest, Ymat[0][0].Scale()*cps.Params.Scale())
				distanceTmp = crypto.CRescale(cps, distanceTmp)

				mutex.Lock()
				initial = crypto.CSub(cps, initial, distanceTmp)
				mutex.Unlock()
			}(j, ncolsY, Ymat, firstXIndex)
		}
		firstXIndex = firstXIndex + ncolsY
		if b == nbrBlocks-1 {
			wg.Wait()
		}
	}
	return initial
}

func Frac_part(decrypte_XY []complex128) []complex128 {
	errToInt := make([]complex128, len(decrypte_XY))
	for i := 0; i < len(decrypte_XY); i++ {
		errToInt[i] = complex(real(decrypte_XY[i])-math.Round(real(decrypte_XY[i])), 0)
	}
	return errToInt
}

func (pi *ProtocolInfo) signTestComposed(cps *crypto.CryptoParams, net *mpc.Network, left crypto.CipherVector, leftPlain crypto.PlainVector, distance crypto.CipherVector, mutex *sync.Mutex, maskEncoded crypto.PlainVector, sourcePid int, saveInOut int, numIter int) crypto.CipherVector {
	if left == nil && leftPlain == nil || left != nil && leftPlain != nil {
		panic("Exactly one ciphertext left or plain left should be supplied")
	}
	// might need to change scaleInd
	tmp := pi.SignTest(cps, net, left, leftPlain, distance, mutex, maskEncoded, sourcePid, pi.approxInt, saveInOut, 2, saveInOut)
	for j := 1; j <= numIter; j++ {
		tmp = pi.SignTestImgRemoved(cps, net, tmp, mutex, maskEncoded, sourcePid, pi.approxInt, j+2, saveInOut)
	}

	return PostprcoessSignTest(cps, tmp, false)
}

func save_array(sign []complex128, filename string, imagine bool, format bool) {
	// get the directory of filename
	dir := filepath.Dir(filename)
	os.MkdirAll(dir, 0777)
	file, err := os.Create(filename)
	if err != nil {
		panic(err)
	}

	var form string
	if format {
		form = "%.30f\n"
	} else {
		form = "%e\n"
	}
	if imagine {
		for _, num := range sign {
			fmt.Fprintf(file, form, imag(num))
		}
	} else {
		for _, num := range sign {
			fmt.Fprintf(file, form, real(num))
		}
	}
	file.Close()
}

func computeMinimum(cps *crypto.CryptoParams, sign crypto.CipherVector, left crypto.CipherVector,
	right crypto.PlainVector) crypto.CipherVector {

	// minimum sign(a-b)*a + ((sign(a-b)-1)*(-b))
	signMinusOne := crypto.CAddConst(cps, sign, -1.0)
	signMinusOne = crypto.CRescale(cps, signMinusOne)
	rightSigned := crypto.CPMult(cps, signMinusOne, right)

	leftSign := crypto.CMult(cps, sign, left) //hetYInvRep)
	minimum := crypto.CAdd(cps, leftSign, rightSigned)

	return minimum
}

func (pi *ProtocolInfo) SignTestImgRemoved(cps *crypto.CryptoParams, net *mpc.Network, sub crypto.CipherVector, mutex *sync.Mutex, maskEncoded crypto.PlainVector, sourcePid int, intv crypto.IntervalApprox, idx int, saveInOut int) crypto.CipherVector {
	sub = crypto.CAdd(cps, sub, crypto.CConjugate(cps, sub))
	crypto.CMultConst(cps, sub, 0.5, true)
	sub = crypto.CRescale(cps, sub)
	return pi.SignTestCore(cps, net, sub, mutex, maskEncoded, sourcePid, intv, idx, saveInOut)
}

func (pi *ProtocolInfo) SignTestCore(cps *crypto.CryptoParams, net *mpc.Network, sub crypto.CipherVector,
	mutex *sync.Mutex, maskEncoded crypto.PlainVector, sourcePid int, intv crypto.IntervalApprox, idx int, saveInOut int) crypto.CipherVector {
	// input should have sufficient level
	subSquare := crypto.CMult(cps, sub, sub)
	subSquare = crypto.CPAdd(cps, subSquare, maskEncoded)

	mutex.Lock()
	netDec := pi.basicProt.MpcObj.GetNetworks()[pi.bootMap[strconv.Itoa(sourcePid)+dec]]
	subSquareInv := lib.CInvSqrtSquareApprox(sourcePid, net, cps, subSquare, intv, netDec)
	mutex.Unlock()

	signOne := crypto.CMult(cps, sub, subSquareInv)

	// UPDATE: easiest to bootstrap here
	for i := range signOne {
		mutex.Lock()
		net.CollectiveBootstrap(cps, signOne[i], sourcePid)
		mutex.Unlock()
	}

	sign := signOne
	return sign
}

func PostprcoessSignTest(cps *crypto.CryptoParams, result crypto.CipherVector, normalize bool) crypto.CipherVector {
	if !normalize {
		return crypto.CAddConst(cps, result, 1)
	} else {
		// maybe can add a check to see if level is too low
		return crypto.CRescale(cps, crypto.CMultConst(cps, crypto.CAddConst(cps, result, 1), 0.5, true))
	}
}

func (pi *ProtocolInfo) SignTest(cps *crypto.CryptoParams, net *mpc.Network, left crypto.CipherVector, leftPlain crypto.PlainVector, right crypto.CipherVector,
	mutex *sync.Mutex, maskEncoded crypto.PlainVector, sourcePid int, intv crypto.IntervalApprox, scale_ind int, idx int, saveInOut int) crypto.CipherVector {
	var sub crypto.CipherVector
	if left == nil {
		sub = lib.CSubPlain(cps, leftPlain, right)
	} else {
		sub = lib.CSub(cps, left, right)
	}
	if pi.Debug_ST_flag {
		// debug mode
		if left != nil {
			var dec_left []complex128
			mutex.Lock()
			dec_left = pi.decryptVectorForDebugging(cps, left, sourcePid)
			mutex.Unlock()
			save_array(dec_left, fmt.Sprintf("signouts/left%d.txt", saveInOut), false, true)
		}

		mutex.Lock()
		dec_right := pi.decryptVectorForDebugging(cps, right, sourcePid)
		mutex.Unlock()
		save_array(dec_right, fmt.Sprintf("signouts/right%d.txt", saveInOut), false, true)
		mutex.Lock()
		dec_sub := pi.decryptVectorForDebugging(cps, sub, sourcePid)
		mutex.Unlock()
		save_array(dec_sub, fmt.Sprintf("signouts/sub%d.txt", saveInOut), false, true)
	}
	if scale_ind != -1 {
		scaleDown := 1.0 / float64(intv.ScaleDowns[scale_ind])
		crypto.CMultConst(cps, sub, scaleDown, true)
		crypto.CRescale(cps, sub)
	}
	return pi.SignTestImgRemoved(cps, net, sub, mutex, maskEncoded, sourcePid, intv, idx, saveInOut)
}

// createControlledOutput computes the output according to what needs to be revealed
func createControlledOutput(cps *crypto.CryptoParams, compResult KinshipResult, reveal int, single bool) (toSend LocalPhase1Result) {
	// should do sign tests on the result if want to reveal degree only
	// reveal == 1 sends 4 boolean values indicating if each of the threhsold tests succeeded
	toSend.Result = []crypto.CipherVector{compResult.Result}
	return
}

// sendControlledOutput sends the comparison result to the other party
func sendControlledOutput(compResult LocalPhase1Result, otherPid int, net *mpc.Network, localTest bool, reveal int) {
	net.SendInt(len(compResult.Result[0]), otherPid) // size of the vector
	net.SendCipherVector(compResult.Result[0], otherPid)
}

func receiveControlledOutput(cps *crypto.CryptoParams, otherPid int, net *mpc.Network, localTest bool, reveal int) (receivedContResult LocalPhase1Result) {
	sizeVec := net.ReceiveInt(otherPid)
	receivedContResult.Result = []crypto.CipherVector{net.ReceiveCipherVector(cps, sizeVec, otherPid)}
	return receivedContResult
}

func (pi *ProtocolInfo) decryptVectorForDebugging(cps *crypto.CryptoParams, ct crypto.CipherVector, pid int) []complex128 {
	network := pi.basicProt.MpcObj.GetNetworks()[pi.bootMap[strconv.Itoa(pid)+dec]]
	return crypto.DecodeFloatVector2(cps, network.CollectiveDecryptVec(cps, ct, pid))[0]
}

func prepareXVectors(X ComparisonDataLocal, bucketSize, nbrBuckets int, scaleDown float64) ([]float64, []float64, []float64, []float64) {
	XSquareElemRepeatCol := make([]float64, 0)
	indexXRepeatCol := make([]float64, 0)
	var XSquareElem float64
	var hetXwithThreshElem float64
	var hetXInvwithThreshElem float64

	hetXwithThreshRepeatCol := make([]float64, 0)

	hetXInvwithThreshRepeatCol := make([]float64, 0)

	for buckIn := 0; buckIn < bucketSize; buckIn++ {
		XSquareElemRepeat := make([]float64, bucketSize)
		hetXwithThreshRepeat := make([]float64, bucketSize)

		hetXInvwithThreshRepeat := make([]float64, bucketSize)
		indexXElemRepeat := make([]float64, bucketSize)
		for bIndex := 0; bIndex < nbrBuckets; bIndex++ {
			XrowIndex := bIndex*bucketSize + buckIn

			hetXwithThreshElem = X.Xhet[XrowIndex]

			hetXInvwithThreshElem = X.XhetInv[XrowIndex]
			XSquareElem = X.XSquare[XrowIndex] * 1.0 / scaleDown

			for bs := 0; bs < bucketSize; bs++ {
				hetXwithThreshRepeat[bs] = hetXwithThreshElem

				hetXInvwithThreshRepeat[bs] = hetXInvwithThreshElem
				XSquareElemRepeat[bs] = XSquareElem
				indexXElemRepeat[bs] = X.Index[XrowIndex]
			}
			hetXwithThreshRepeatCol = append(hetXwithThreshRepeatCol, hetXwithThreshRepeat...)

			hetXInvwithThreshRepeatCol = append(hetXInvwithThreshRepeatCol, hetXInvwithThreshRepeat...)
			XSquareElemRepeatCol = append(XSquareElemRepeatCol, XSquareElemRepeat...)
			indexXRepeatCol = append(indexXRepeatCol, indexXElemRepeat...)
		}
	}
	return hetXwithThreshRepeatCol, hetXInvwithThreshRepeatCol, XSquareElemRepeatCol, indexXRepeatCol
}

func (pi *ProtocolInfo) compute2XYColumn(cps *crypto.CryptoParams, X ComparisonDataLocal, Y crypto.CipherMatrix, bucketSize, nbrBuckets, nbrRepeatY, columnIndexY, columnIndexX int,
	cache gwas.CipherVectorAccV2, pid int) gwas.CipherVectorAccV2 {
	XElemRepeatCol := make([]float64, 0)
	for buckIn := 0; buckIn < bucketSize; buckIn++ { // bucket
		XelemRepeat := make([]float64, bucketSize)

		for bIndex := 0; bIndex < nbrBuckets; bIndex++ {
			XrowIndex := bIndex*bucketSize + buckIn
			Xelem := X.X.At(XrowIndex, columnIndexX) * 2
			for bs := 0; bs < bucketSize; bs++ {
				XelemRepeat[bs] = Xelem
			}
			XElemRepeatCol = append(XElemRepeatCol, XelemRepeat...)
		}
	}

	// for some reason YRep always contains the same noises?
	YRep := make(crypto.CipherVector, 0)
	for yrepeat := 0; yrepeat < nbrRepeatY; yrepeat++ {
		YRep = append(YRep, Y[columnIndexY]...)
	}
	XElemRepeatEncode, _ := crypto.EncodeFloatVectorWithScale(cps, XElemRepeatCol, YRep[0].Scale())
	gwas.ToMontgomeryForm(cps, XElemRepeatEncode)

	gwas.CPMultAccWithoutMRedV2(YRep, XElemRepeatEncode, cache)
	return cache

}

func prepareYVectors(Y ComparisonDataOther, nbrRepeatY int) (crypto.CipherVector, crypto.CipherVector, crypto.CipherVector, crypto.CipherVector) {

	YSquareRep := make(crypto.CipherVector, 0)
	YhetRep := make(crypto.CipherVector, 0)

	YhetInvRep := make(crypto.CipherVector, 0)
	YindexRep := make(crypto.CipherVector, 0)
	for yrepeat := 0; yrepeat < nbrRepeatY; yrepeat++ {
		YSquareRep = append(YSquareRep, Y.YSquare...)
		YhetRep = append(YhetRep, Y.Yhet...)
		YhetInvRep = append(YhetInvRep, Y.YhetInv...)
		YindexRep = append(YindexRep, Y.Yindex...)
	}

	return YSquareRep, YhetRep, YhetInvRep, YindexRep
}

func square(i, j int, value float64) float64 {
	return value * value
}
