package relativeMatch

// usage directly make X or make Y or make Z on the machines
// usage: bash run_local.sh :
// export PARA=20 // reduce if there's not enough cores or memory
// export PID=1 // adjust when running other parties: 0 for a background helper helpert that does nothing (except during initialization and termination), 1 for X, 2 for Y
// export t="demo"
// export FOLDER="config_local/$t/"

import (
	"math"
	"os"
	"relativeMatch/internal"
	"strconv"
	"sync"
	"testing"
	"time"

	"go.dedis.ch/onet/v3/log"
)

var pid, _ = strconv.Atoi(os.Getenv("PID"))
var folder = os.Getenv("FOLDER")
var paraRun int

func (pi *ProtocolInfo) inferParams(hard_reset bool) {
	// if pid > 0 {
	if !hard_reset {
		rowIndexVector, colIndexVector := internal.ReadNPZVectors(pi.rowIndexFile, pi.columnIndexFile)
		if pi.numberOfColumns <= 0 {
			pi.numberOfColumns = pi.M
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
		// TODO: check which scaling makes more sense if number of columns is small
		pi.scaleDown = 1.0
		// pi.scaleDown = pi.scaledownLocal
		// assert that numberOfColumns >= numberOfColumnsTest
		if pi.numberOfColumns < pi.numberOfColumnsTest {
			log.Fatal("numberOfColumns < numberOfColumnsTest")
		}
	} else {
		rowIndexVector, colIndexVector := internal.ReadNPZVectors(pi.rowIndexFile, pi.columnIndexFile)
		pi.numberOfColumns = pi.M
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
	if folder == "" {
		panic("FOLDER environment variable not set")
	}
	if os.Getenv("PID") == "" {
		panic("PID environment variable not set")
	}

	timeTotal := time.Now()
	prot := InitializeRelativeMatchingProtocol(pid, configFolder, nil)
	prot.inferParams(false)
	reportStats(timeTotal, prot, "init protocol:")
	paraRun = prot.PARA

	if prot.TestSignTest == 0 {
		// test the main setting of the protocol
		globalResult := runPhase1withTime(prot, configFolder)
		if pid > 0 {
			runPhase2WithTime(prot, globalResult)
		}
	} else {
		panic("Other testing modes are removed")
	}
	log.LLvl1("closing")
	prot.SyncAndTerminate(true)
}

func runPhase1withTime(prot *ProtocolInfo, configFolder string) GlobalPhase1Result {
	timeStart := time.Now()
	reportStats(timeStart, prot, " start phase 1: ")
	globalResult := startProtocols(prot, configFolder)
	reportStats(timeStart, prot, " finish phase 1: ")
	return globalResult
}

func runPhase2WithTime(prot *ProtocolInfo, globalResult GlobalPhase1Result) {
	timeStart := time.Now()
	reportStats(timeStart, prot, " start phase 2 of MHE : ")
	if prot.reveal == 0 || prot.reveal == 1 || prot.reveal == 2 {
		// by default, accumulate every BatchCVec stored
		log.LLvl1("accumulating " + strconv.Itoa(len(globalResult.Result)) + " BatchCVecs")
		for i := 0; i < len(globalResult.Result); i++ {
			prot.accumulateByID(globalResult.Result[i], pid, i)
		}
	}
	reportStats(timeStart, prot, " perform accumulation (phase 2) protocol: ")
}

func reportStats(timeStart time.Time, prot *ProtocolInfo, part_name string) {
	log.LLvl1(timeStart, "========= time at "+part_name, time.Since(timeStart))
	prot.basicProt.MpcObj.GetNetworks().PrintNetworkLog()
}

func startProtocols(prot *ProtocolInfo, configFolder string) (globalResult GlobalPhase1Result) {
	// partitionSize is the number of rows that each subprocess will process
	// which is the total number of rows divided by the number of subprocesses (paraRun or PARA in the bash script)
	// (rounded up to the nearest multiple of prot.batchLength = B = how many values are encrypted in each ciphertext)
	partitionSize := int(math.Ceil((float64(prot.totalNbrOfRowsTest)/float64(paraRun))/float64(prot.batchLength))) * prot.batchLength

	globalWg := sync.WaitGroup{}
	mutex := sync.Mutex{}
	if prot.useMPC {
		prot.scaleDown = 1.0
	}

	globalWg.Add(paraRun)
	// verify that numThreads >= PARA * ((NumMainParties * 3) + 1)
	if prot.basicProt.Config.MpcNumThreads < paraRun*((prot.basicProt.Config.NumMainParties*3)+1) {
		log.Panic("numThreads must be >= paraRun * ((NumMainParties * 3) + 1)")
	}
	last_start := prot.startKey // corrected to start_key instead of 0
	numArrays := []int{1, len(prot.threshValue), len(prot.discretizedThresh), 0}[prot.reveal]
	globalResult.Result = make([]BatchedCVec, numArrays)
	for i := 0; i < numArrays; i++ {
		globalResult.Result[i] = make(BatchedCVec)
	}
	for exec := 0; exec < paraRun; exec++ {
		//execThreads := paraRun * 3
		go func(exec int, last_start int, prot *ProtocolInfo, mutex *sync.Mutex) {
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
				discretizedThresh:   prot.discretizedThresh,
				bucketSize:          prot.bucketSize,
				blinding:            prot.blinding,
				single:              prot.single,
				M:                   prot.M,
				reveal:              prot.reveal,
				simpleDataPath:      prot.simpleDataPath,
				startingIndex:       prot.startingIndex,
				numberOfColumns:     prot.numberOfColumns,
				numberOfColumnsTest: prot.numberOfColumnsTest,
				npz:                 prot.npz,
				totalNbrOfRows:      prot.totalNbrOfRows,
				totalNbrOfRowsTest:  prot.totalNbrOfRowsTest,
				blockLimit:          prot.blockLimit,
				startKey:            last_start,
				scaledownLocal:      prot.scaledownLocal,
				endKey:              newEndKey,
				rowIndexFile:        prot.rowIndexFile,
				columnIndexFile:     prot.columnIndexFile,
				batchLength:         prot.batchLength,
				outFolder:           prot.outFolder,
				queryLength:         prot.queryLength,
				localTest:           prot.localTest,
				useMPC:              prot.useMPC,
				net:                 exec,
			}
			sending := exec < paraRun/2
			if pid == 1 {
				sending = !sending
			}
			partialGlobalResult := newProt.BatchProtocols(configFolder, sending, prot.startKey, mutex)
			// assert that the numbers match
			if numArrays != len(partialGlobalResult.Result) {
				panic("numArrays != len(partialGlobalResult)")
			}
			// need to append based on the number of computed vectors
			for i := 0; i < len(partialGlobalResult.Result); i++ {
				globalResult.Result[i][exec] = partialGlobalResult.Result[i]
			}
		}(exec, last_start, prot, &mutex)
		last_start += partitionSize
	}
	globalWg.Wait()
	return
}
