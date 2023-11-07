package relativeMatch

// usage directly make X or make Y or make Z on the machines
// usage: bash run_local.sh :
// export PARA=20 // reduce if there's not enough cores or memory
// export PID=1 // adjust when running other parties: 0 for a background helper helpert that does nothing (except during initialization and termination), 1 for X, 2 for Y
// export t="demo"
// export FOLDER="config_local/$t/"

import (
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"relativeMatch/internal"
	"strconv"
	"sync"
	"testing"
	"time"

	"github.com/hhcho/sfgwas/crypto"
	"go.dedis.ch/onet/v3/log"
)

var pid, _ = strconv.Atoi(os.Getenv("PID"))
var folder = os.Getenv("FOLDER")
var paraRun, _ = strconv.Atoi(os.Getenv("PARA"))

func readNumSNPs(pi *ProtocolInfo, sketch_id int) int {
	dir, _ := filepath.Split(pi.columnIndexFile)
	fname := dir + fmt.Sprintf("num_SNP%d.txt", sketch_id)
	content, err := ioutil.ReadFile(fname)
	if err != nil {
		log.Fatal(err)
	}
	numberOfColumns, err := strconv.Atoi(string(content))
	if err != nil {
		log.Fatal(err)
	}
	log.LLvl1("Number of columns: ", numberOfColumns)
	return numberOfColumns
}

func (pi *ProtocolInfo) inferParams(hard_reset bool) {
	// if pid > 0 {
	if !hard_reset {
		rowIndexVector, colIndexVector := internal.ReadNPZVectors(pi.rowIndexFile, pi.columnIndexFile)
		if pi.numberOfColumns <= 0 {
			pi.numberOfColumns = readNumSNPs(pi, 100)
		}
		// assert that numberOfColumns >= numberOfColumnsTest
		if pi.numberOfColumns < pi.numberOfColumnsTest {
			log.Fatal("numberOfColumns < numberOfColumnsTest")
		}
		log.LLvl1("Scale down factor: ", pi.scaleDown)
		if pi.numberOfColumnsTest <= 0 {
			pi.numberOfColumnsTest = len(colIndexVector)
		}
		if pi.totalNbrOfRowsTest <= 0 {
			pi.totalNbrOfRowsTest = len(rowIndexVector)
		}
		if pi.totalNbrOfRows <= 0 {
			pi.totalNbrOfRows = len(rowIndexVector)
		}
		pi.scaledownLocal = pi.scaleDown
		if pi.scaledownLocal <= 1e-6 {
			pi.scaledownLocal = float64(pi.numberOfColumnsTest) / (10.0)
		}
		pi.scaleDown = 1.0
	} else {
		rowIndexVector, colIndexVector := internal.ReadNPZVectors(pi.rowIndexFile, pi.columnIndexFile)
		pi.numberOfColumns = readNumSNPs(pi, 100)
		if pi.numberOfColumns < pi.numberOfColumnsTest {
			log.Fatal("numberOfColumns < numberOfColumnsTest")
		}
		log.LLvl1("Scale down factor: ", pi.scaleDown)
		pi.numberOfColumnsTest = len(colIndexVector)
		// pi.totalNbrOfRowsTest = len(rowIndexVector)
		pi.totalNbrOfRows = len(rowIndexVector)
		pi.scaledownLocal = float64(pi.numberOfColumnsTest) / (10.0)
		pi.scaleDown = 1.0
	}
}

func TestRelativeSearchProtocol(t *testing.T) {
	configFolder := folder

	timeTotal := time.Now()
	prot := InitializeRelativeMatchingProtocol(pid, configFolder, nil)
	prot.inferParams(false)
	reportStats(timeTotal, prot, "init protocol:")

	if prot.TestSignTest == 0 {
		// test the main setting of the protocol
		globalResult := runPhase1withTime(prot, configFolder)
		if pid > 0 && prot.reveal == 0 {
			timeTotal = time.Now()
			runPhase2WithTime(prot, globalResult)
		}
	} else {
		panic("Other testing modes are removed")
	}
	log.LLvl1("closing")
	prot.SyncAndTerminate(true)
}

func runPhase1withTime(prot *ProtocolInfo, configFolder string) ResultOutput {
	timeStart := time.Now()
	reportStats(timeStart, prot, " start phase 1: ")
	globalResult := startProtocols(prot, configFolder)
	reportStats(timeStart, prot, " finish phase 1: ")
	return globalResult
}

func runFakePhase1withTime(prot *ProtocolInfo, cps *crypto.CryptoParams) ResultOutput {
	timeStart := time.Now()
	reportStats(timeStart, prot, " start encrypting 0 in a fakse phase 1: ")
	globalResult := makeZerosForAccuTest(prot, cps, pid)
	reportStats(timeStart, prot, " finish encrypting 0 in a fakse phase 1: ")
	return globalResult
}

func runPhase2WithTime(prot *ProtocolInfo, globalResult ResultOutput) {
	timeStart := time.Now()
	reportStats(timeStart, prot, " start phase 2 of MHE : ")
	prot.accumulateByID(globalResult.AllResult, pid)
	reportStats(timeStart, prot, " perform accumulation (phase 2) protocol: ")
}

func terminateAndRestartProtocol(prot *ProtocolInfo, configFolder string) *ProtocolInfo {
	log.LLvl1("closing")
	prot.SyncAndTerminate(true)
	prot = InitializeRelativeMatchingProtocol(pid, configFolder, nil)
	prot.inferParams(false)
	return prot
}

func reportStats(timeStart time.Time, prot *ProtocolInfo, part_name string) {
	log.LLvl1(timeStart, "========= time at "+part_name, time.Since(timeStart))
	prot.basicProt.MpcObj.GetNetworks().PrintNetworkLog()
}

func makeZerosForAccuTest(prot *ProtocolInfo, cps *crypto.CryptoParams, pid int) ResultOutput {
	// make encryptions of 0 to test accumuluation (MHE-Phase 2)'s running time only
	if pid > 0 && !prot.useMPC {
		prot.alwaysDecrypting()
		prot.alwaysBootstrapping()
	}
	partitionSize := int(math.Ceil((float64(prot.totalNbrOfRowsTest)/float64(paraRun))/float64(prot.batchLength))) * prot.batchLength
	last_start := prot.startKey
	var globalResult ResultOutput
	globalResult.AllResult = make(map[int][]crypto.CipherVector)
	for exec := 0; exec < paraRun; exec++ {
		newEndKey := last_start + partitionSize
		if newEndKey > prot.startKey+prot.totalNbrOfRowsTest {
			newEndKey = prot.startKey + prot.totalNbrOfRowsTest
		}
		var partialGlobalResult ResultOutput
		partialGlobalResult.Result = make([]crypto.CipherVector, 0)
		numBlocks := int(math.Ceil(float64(newEndKey-last_start) / float64(prot.batchLength)))
		for i := 0; i < numBlocks; i++ {
			partialGlobalResult.Result = append(partialGlobalResult.Result, crypto.DropLevelVec(cps, crypto.CZeros(cps, 1), 2))
		}
		if prot.reveal == 0 {
			globalResult.AllResult[exec] = partialGlobalResult.Result
		}
		last_start += partitionSize
	}
	return globalResult
}

func startProtocols(prot *ProtocolInfo, configFolder string) (globalResult ResultOutput) {
	// partitionSize is the number of rows that each subprocess will process
	// which is the total number of rows divided by the number of subprocesses (paraRun or PARA in the bash script)
	// (rounded up to the nearest multiple of prot.batchLength = B = how many values are encrypted in each ciphertext)
	partitionSize := int(math.Ceil((float64(prot.totalNbrOfRowsTest)/float64(paraRun))/float64(prot.batchLength))) * prot.batchLength

	globalWg := sync.WaitGroup{}
	if prot.useMPC {
		prot.scaleDown = 1.0
	}

	globalWg.Add(paraRun)
	// verify that numThreads >= PARA * ((NumMainParties * 3) + 1)
	if prot.basicProt.Config.MpcNumThreads < paraRun*((prot.basicProt.Config.NumMainParties*3)+1) {
		log.Panic("numThreads must be >= paraRun * ((NumMainParties * 3) + 1)")
	}
	last_start := prot.startKey // corrected to start_key instead of 0
	globalResult.AllResult = make(map[int][]crypto.CipherVector)
	for exec := 0; exec < paraRun; exec++ {
		//execThreads := paraRun * 3
		go func(exec int, last_start int, prot *ProtocolInfo) {
			defer globalWg.Done()
			newBootMap := make(map[string]int, 0)
			for i, v := range prot.bootMap {
				// this ensures that parallel runs are not using the same networks
				newBootMap[i] = v + exec*((prot.basicProt.Config.NumMainParties*3)+1)
			}
			newEndKey := last_start + partitionSize
			if newEndKey > prot.startKey+prot.totalNbrOfRowsTest {
				newEndKey = prot.startKey + prot.totalNbrOfRowsTest
			}
			log.LLvl1(exec, "using network: ", newBootMap)
			newProt := &ProtocolInfo{ // copy of the protocol with some values changed
				basicProt:           prot.basicProt,
				comparisonMap:       prot.comparisonMap,
				bootMap:             newBootMap,
				Debug_ST_flag:       prot.Debug_ST_flag,
				approxInt:           prot.approxInt,
				scaleDown:           prot.scaleDown,
				numThreads:          prot.numThreads,
				threshValue:         prot.threshValue,
				bucketSize:          prot.bucketSize,
				blinding:            prot.blinding,
				single:              prot.single,
				reveal:              prot.reveal,
				simpleDataPath:      prot.simpleDataPath,
				resultFolder:        prot.resultFolder,
				startingIndex:       prot.startingIndex,
				numberOfColumns:     prot.numberOfColumns,
				numberOfColumnsTest: prot.numberOfColumnsTest,
				npz:                 prot.npz,
				separator:           prot.separator,
				totalNbrOfRows:      prot.totalNbrOfRows,
				totalNbrOfRowsTest:  prot.totalNbrOfRowsTest,
				blockLimit:          prot.blockLimit,
				startKey:            last_start,
				scaledownLocal:      prot.scaledownLocal,
				endKey:              newEndKey,
				rowIndexFile:        prot.rowIndexFile,
				columnIndexFile:     prot.columnIndexFile,
				batchLength:         prot.batchLength,
				queryLength:         prot.queryLength,
				localTest:           prot.localTest,
				useMPC:              prot.useMPC,
				net:                 exec,
			}
			sending := exec < paraRun/2
			if pid == 1 {
				sending = !sending
			}
			partialGlobalResult := newProt.BatchProtocols(configFolder, sending, prot.startKey)
			if prot.reveal == 0 {
				globalResult.AllResult[exec] = partialGlobalResult.Result
			}
		}(exec, last_start, prot)
		last_start += partitionSize
	}
	globalWg.Wait()
	return
}
