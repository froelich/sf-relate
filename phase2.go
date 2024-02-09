package relativeMatch

import "C"
import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"strconv"
	"sync"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/sbinet/npyio/npz"
	"go.dedis.ch/onet/v3/log"
)

func (pi *ProtocolInfo) accumulateByID(allResultsToCombine map[int]crypto.CipherVector, pid int, name int) (rst crypto.CipherVector) {
	cps := pi.basicProt.Cps

	// read in id
	ids := pi.ReadIDs()
	// record where each id appears in a linked list
	idLocations := make([][]int, pi.n)
	// go over each id
	for i := range idLocations {
		idLocations[i] = make([]int, 0)
	}
	for idx := pi.startKey; idx < pi.startKey+pi.totalNbrOfRowsTest; idx++ {
		// omit the parts that are not tested
		id := ids[idx]
		if idx >= pi.startKey+pi.totalNbrOfRowsTest {
			break
		}
		if id != -1 {
			idLocations[id] = append(idLocations[id], idx)
		}
	}

	paraRun := len(allResultsToCombine)
	partitionSize := int(math.Ceil((float64(pi.totalNbrOfRowsTest)/float64(paraRun))/float64(pi.batchLength))) * pi.batchLength
	// go over each id, find the blocks that it is in ....
	idToCheck := make([]int, 0)
	compactReveal := true
	for id, loc := range idLocations {
		if len(loc) == 0 && compactReveal {
			continue
		}
		idToCheck = append(idToCheck, id)
	}

	// need to compute the correct boundary conditions of size of blocks
	numSlots := cps.Params.Slots()
	OutputBlockCnt := int(math.Ceil(float64(len(idToCheck)) / float64(numSlots)))
	lastBlockSize := len(idToCheck) - (OutputBlockCnt-1)*numSlots

	scaleDown := 1.0 / 40
	mask := make([]float64, numSlots)
	ptFive := make([]float64, numSlots)
	for i := range ptFive {
		// the value of ptFive is 3.5 if in reveal == 0 mode, otherwise it is 0.5
		if pi.reveal == 0 {
			ptFive[i] = 3.5 * scaleDown
		} else {
			ptFive[i] = 0.5 * scaleDown
		}
	}
	for i := numSlots - 1; i >= lastBlockSize; i-- {
		mask[i] = 0.0005 * scaleDown
	}
	maskEncoded, _ := crypto.EncodeFloatVector(cps, mask)
	pointFiveEncoded, _ := crypto.EncodeFloatVector(cps, ptFive)

	// this is the all zero vector
	zeroEncoded := EncodeMask(pi, nil, numSlots, 1.0)

	// processing block by block
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	log.LLvl1("Starting out blocks processing")
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	var mutex sync.Mutex
	for obid := 0; obid < OutputBlockCnt; obid++ {
		// spawn worker processes
		oblock := crypto.CZeros(cps, 1)
		obegin := obid * numSlots
		var oend int
		var maskForOut crypto.PlainVector
		if obid == OutputBlockCnt-1 {
			oend = lastBlockSize
			maskForOut = maskEncoded
		} else {
			oend = numSlots
			maskForOut = zeroEncoded
		}
		// parallelize here
		numThread := 128
		output_idx_to_process := make(chan int, numThread)
		var local_wg sync.WaitGroup
		var local_mutex sync.Mutex
		for thread_id := 0; thread_id < numThread; thread_id++ {
			go func() {
				for {
					output_idx, ok := <-output_idx_to_process
					if !ok {
						return
					}
					partialIDResult := crypto.CZeros(cps, 1)
					id := idToCheck[obegin+output_idx]
					loc := idLocations[id]
					if output_idx%50 == 0 {
						log.LLvl1("Processing the ID at index ", output_idx)
						log.LLvl1("len of loc = ", len(loc))
					}
					maskOutputEncoded := EncodeMask(pi, nil, output_idx, 1.0)
					lastStart := -1
					for _, idx := range loc {
						// remove duplicated computation for blocks
						// find which batch to read fro the allResultsToCombine
						// should be floor(idx / partitionSize)
						// WARNING: pi.StartKey should be 0 for this to work --- otherwise need to find out what ID the first block is
						if pi.startKey != 0 {
							panic("Accumulation does not work if parties do not compute kinship starting from position 0 in table")
						}
						batch_id := int(float64(idx) / float64(partitionSize))
						// need to find which block it is too
						block_id := int(float64(int(idx)-batch_id*partitionSize) / float64(pi.batchLength))
						// find the corresponding result
						rstCt := allResultsToCombine[batch_id][block_id]

						// construct a mask for id
						// namely only set mask[idx] = 1 when block[idx] == id
						// need to add the batchStart shift too
						batchStart := batch_id * partitionSize
						block_start := batchStart + block_id*pi.batchLength
						if block_start == lastStart {
							continue
						} else {
							lastStart = block_start
						}
						// End needs to be < len of ids
						// will need to append with -1
						block_end := block_start + pi.batchLength
						if block_end > len(ids) {
							block_end = len(ids)
						}
						block := ids[block_start:block_end]
						// note that scale is discarded here
						maskEncoded := EncodeMask(pi, block, id, scaleDown)

						// pt-wise multiply the two together so we get one result
						// first multiplication --- check level later
						crypto.DecodeFloatVector(cps, maskEncoded)
						extracted := crypto.CPMult(cps, crypto.CipherVector{rstCt}, maskEncoded)
						crypto.CRescale(cps, extracted)
						partialIDResult = crypto.CAdd(cps, partialIDResult, extracted)

					}
					// perform innerSum to get one result only
					// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
					foundInBlock := crypto.CipherVector{crypto.InnerSumAll(cps, partialIDResult)}
					toSave := crypto.CPMult(cps, foundInBlock, maskOutputEncoded)
					crypto.CRescale(cps, toSave)
					local_mutex.Lock()
					oblock = crypto.CAdd(cps, oblock, toSave)
					local_mutex.Unlock()
					local_wg.Done()
				}
			}()
		}
		local_wg.Add(oend)
		for output_idx := 0; output_idx < oend; output_idx++ {
			output_idx_to_process <- output_idx
		}
		close(output_idx_to_process)
		local_wg.Wait()
		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		// do a final sign test before decrypting

		ownNetworkBoot := pi.basicProt.MpcObj[pi.bootMap[strconv.Itoa(pid)+boot]].Network
		ownNetworkBoot.CollectiveBootstrap(cps, oblock[0], pid)
		oblock = pi.signTestComposed(cps, ownNetworkBoot, nil, pointFiveEncoded, oblock, &mutex, maskForOut, pid, 0, 1)
		decrypted := pi.decryptVectorForDebugging(cps, oblock, pid)

		// print out the block result here
		if pi.reveal == 0 {
			save_decryption(pi, strconv.Itoa(name), obid, pid, oend, idToCheck, obegin, decrypted)
		}
		if pi.reveal == 1 {
			save_decryption(pi, "degree"+strconv.Itoa(name), obid, pid, oend, idToCheck, obegin, decrypted)
		} else if pi.reveal == 2 {
			save_decryption(pi, "kinship"+strconv.Itoa(name), obid, pid, oend, idToCheck, obegin, decrypted)
		} else if pi.reveal == 3 {
			panic("reveal == 3 does not do accumulation")
		} else {
			panic("unknown mode of revelation")
		}
	}
	return
}

func save_decryption(pi *ProtocolInfo, name string, obid int, pid int, oend int, idToCheck []int, obegin int, decrypted []complex128) {
	log.LLvl1("saving ", obid, "block")
	folder := os.Getenv("FOLDER")
	filename := folder + name + "_" + strconv.Itoa(obid) + "_party" + strconv.Itoa(pid) + ".csv"
	csvOut, _ := os.Create(filename)
	writer := csv.NewWriter(csvOut)
	writer.Write([]string{"ID", "Detected"})
	for output_idx := 0; output_idx < oend; output_idx++ {
		writer.Write([]string{fmt.Sprintf("%d", idToCheck[obegin+output_idx]), fmt.Sprintf("%.8f", decrypted[output_idx])})
	}
	writer.Flush()
}

func EncodeMask(pi *ProtocolInfo, block []int32, id int, scale float64) crypto.PlainVector {
	mask := make([]float64, pi.batchLength)
	if block != nil {
		for idx := range mask {
			if idx < len(block) && int(block[idx]) == id {
				mask[idx] = scale
			} else {
				mask[idx] = 0.0
			}
		}
	} else {
		for idx := range mask {
			if idx == id {
				mask[idx] = scale
			} else {
				mask[idx] = 0.0
			}
		}
	}
	maskEncoded, _ := crypto.EncodeFloatVector(pi.basicProt.Cps, mask)
	return maskEncoded
}

func (pi *ProtocolInfo) ReadIDs() []int32 {
	// read rows indexes
	f, err := os.Open(pi.rowIndexFile)
	if err != nil {
		log.Fatalf("could not open npz file: %+v", err)
	}

	stat, err := f.Stat()
	if err != nil {
		log.Fatalf("could not stat npz file: %+v", err)
	}

	r, err := npz.NewReader(f, stat.Size())
	if err != nil {
		log.Fatalf("could not open npz archive: %+v", err)
	}

	fRows, rRows := f, r

	defer fRows.Close()
	var mRows []int32
	// assume a single key
	name := rRows.Keys()[0]
	log.LLvl1("Postprocess: reading file with rows index: ", name, rRows.Header(name).Descr.Shape)
	err = npz.Read(fRows, name, &mRows)
	if err != nil {
		log.Error(err)
	}

	return mRows
}
