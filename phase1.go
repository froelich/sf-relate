package relativeMatch

import (
	"math"
	"relativeMatch/lib"
	"sync"
	"time"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/gwas"
	"github.com/hhcho/sfgwas/mpc"
	"go.dedis.ch/onet/v3/log"
)

func (pi *ProtocolInfo) computeKinshipHE(sourcePid int, net *mpc.Network, X ComparisonDataLocal,
	Y ComparisonDataOther, netMat *mpc.Network, otherPid int) KinshipResult {

	cps := pi.basicProt.Cps
	xRows, xCols := X.X.Dims()

	// Prepare comparison result
	var distance crypto.CipherVector
	distance, _ = lib.InitEncryptedVector(cps, xRows*pi.bucketSize, float64(cps.Params.Qi()[cps.Params.Levels()-1])) //*cps.Params.Scale())

	var hetXRepeatCol, hetXInvRepeatCol, XSquareElemRepeatCol []float64
	var hetYRep, hetYInvRep crypto.CipherVector
	var YSquareRep crypto.CipherVector

	timeComputeDistance := time.Now()
	// number of repetitions for each row of X is bucketSize
	// number of repetitions for y
	nbrRepeatY := int(math.Ceil(float64(xRows*pi.bucketSize) / float64(Y.NbrRows*Y.Repeated)))
	vectorSize := int(math.Ceil(float64(xRows*pi.bucketSize) / float64(cps.Params.Slots())))
	nbrBuckets := int(math.Ceil(float64(Y.NbrRows) / float64(pi.bucketSize)))
	log.LLvl1("Y has to be repeated ", nbrRepeatY, " more times")
	log.LLvl1("Each vector has ", vectorSize)
	log.LLvl1("Number of buckets ", nbrBuckets)

	log.LLvl1(sourcePid, " prepare local data ")
	timePrepareX := time.Now()
	hetXRepeatCol, hetXInvRepeatCol, XSquareElemRepeatCol, _ = prepareXVectors(X, pi.bucketSize, nbrBuckets, pi.scaleDown)
	hetXInvRepeatCol = gwas.ScaleF(hetXInvRepeatCol, -1.0)

	log.LLvl1(sourcePid, " prepared local data, TIME:", time.Since(timePrepareX))

	log.LLvl1(sourcePid, " prepare other node data ")
	timePrepareY := time.Now()
	// check if hetYInvRep is correctly scaled
	YSquareRep, hetYRep, hetYInvRep, _ = prepareYVectors(Y, nbrRepeatY)
	log.LLvl1(sourcePid, " prepared other node data, TIME:", time.Since(timePrepareY))

	log.LLvl1(sourcePid, " computes XY ")
	timeComputeXY := time.Now()

	// compute distance through columns while receiving matrix
	nbrBlocks := int(math.Ceil(float64(xCols) / float64(pi.blockLimit)))
	distance = pi.receiveMatrixBlocksAndAggregate(cps, distance, nbrBlocks, netMat, otherPid, pi.numThreads, X, xRows, pi.bucketSize, nbrBuckets, nbrRepeatY)
	log.LLvl1(sourcePid, " computed XY, TIME: ", time.Since(timeComputeXY))

	// add x^2 and y^2 sums
	XSquareElemRepeatEncode, _ := crypto.EncodeFloatVectorWithScale(cps, XSquareElemRepeatCol, cps.Params.Scale())
	log.LLvl1("distance ", distance[0].Scale())
	log.LLvl1("X^2 ", XSquareElemRepeatEncode[0].Scale())

	distance = crypto.CPAdd(cps, distance, XSquareElemRepeatEncode)
	distance = crypto.CAdd(cps, distance, YSquareRep)

	//Rescale
	log.LLvl1("distance level: ", distance[0].Level())
	for i := 0; i < len(distance); i++ {
		net.CollectiveBootstrap(cps, distance[i], sourcePid)
	}
	distance = lib.CRescale(cps, distance)
	log.LLvl1("distance level: ", distance[0].Level())
	crypto.CMultConst(cps, distance, 1.0/float64(pi.scaledownLocal), true)
	crypto.CRescale(cps, distance)
	log.LLvl1(sourcePid, " computed distance, TIME:", time.Since(timeComputeDistance))
	pi.basicProt.MpcObj.GetNetworks().PrintNetworkLog()

	// repeat het if needed
	hetXRepEncode, _ := crypto.EncodeFloatVector(cps, hetXRepeatCol)
	hetXInvRepEncode, _ := crypto.EncodeFloatVector(cps, hetXInvRepeatCol)

	log.LLvl1("Sign testing ... ")
	numIter := 3

	XCompResult := computeSignTestResults(pi, net, hetXRepEncode, distance, xRows, sourcePid, numIter, hetYRep, hetYInvRep, hetXInvRepEncode)
	pi.basicProt.MpcObj.GetNetworks().PrintNetworkLog()

	return XCompResult
}

func computeSignTestResults(pi *ProtocolInfo, net *mpc.Network, hetXRepEncode crypto.PlainVector, distance crypto.CipherVector, xRows int, sourcePid int, numIter int, hetYRep crypto.CipherVector, hetYInvRep crypto.CipherVector, hetXInvRepEncode crypto.PlainVector) (XCompResult KinshipResult) {
	cps := pi.basicProt.Cps
	// compare with the
	mask := make([]float64, cps.Params.Slots()*len(hetYRep))
	for i := xRows * pi.bucketSize; i < cps.Params.Slots(); i++ {
		mask[i] = 1.0 / pi.scaledownLocal
	}
	maskEncoded, _ := crypto.EncodeFloatVectorWithScale(cps, mask, float64(cps.Params.Qi()[distance[0].Level()-1]))
	// append this comparison result
	// the list of kinship thresholds to be compared is pi.threshValue for reveal == 1
	// and [4 - 0.01 * i for i in range(1, 101)] for reveal == 2
	// to be multiplied on the left
	// append the raw kinship computed in reveal == 3
	wg := sync.WaitGroup{}
	var mutex sync.Mutex

	if pi.reveal == 0 {
		// only do two sign tests and multiply them together
		var rstHetX, rstHetY crypto.CipherVector
		wg.Add(2)
		for i := 0; i < 2; i++ {
			go func(i int) {
				defer wg.Done()
				if i == 0 {
					rstHetX = pi.signTestComposed(cps, net, nil, hetXRepEncode, distance, &mutex, maskEncoded, sourcePid, i, numIter)
				} else {
					rstHetY = pi.signTestComposed(cps, net, hetYRep, nil, distance, &mutex, maskEncoded, sourcePid, i, numIter)
				}
			}(i)
		}
		wg.Wait()

		rst := crypto.CMult(cps, rstHetX, rstHetY)
		rst = crypto.CRescale(cps, rst)
		rst = crypto.DropLevelVec(cps, rst, 2)
		rst = crypto.CRescale(cps, rst)

		XCompResult.Result = append(XCompResult.Result, rst[0])
	} else if pi.reveal == 1 || pi.reveal == 2 || pi.reveal == 3 {
		signDiffHet := pi.SignTest(cps, net, nil, hetXRepEncode, hetYRep, &mutex, maskEncoded, sourcePid, pi.approxInt, -1, 0, 20)
		// need to confirm the Levels here
		signDiffHet = PostprcoessSignTest(cps, signDiffHet, true)
		maxOfInvHet := computeMinimum(cps, signDiffHet, hetYInvRep, hetXInvRepEncode)
		kinshipCoeffs := crypto.CMult(cps, maxOfInvHet, distance)
		kinshipCoeffs = crypto.CRescale(cps, kinshipCoeffs)
		if pi.reveal == 1 || pi.reveal == 2 {
			var kinshipThres []float64
			if pi.reveal == 1 {
				kinshipThres = pi.threshValue
			} else if pi.reveal == 2 {
				kinshipThres = pi.discretizedThresh
			}

			numOutputs := len(kinshipThres)
			XCompResult.Result = make(crypto.CipherVector, numOutputs)
			for j := 0; j < numOutputs; j++ {
				wg.Add(1)
				go func(idx int) {
					compLeft := make([]float64, pi.batchLength)
					for i := 0; i < pi.batchLength; i++ {
						compLeft[i] = kinshipThres[idx]
					}
					leftEncode, _ := crypto.EncodeFloatVector(cps, compLeft)

					// bootstrap before sign test
					for i := 0; i < len(kinshipCoeffs); i++ {
						mutex.Lock()
						net.CollectiveBootstrap(cps, kinshipCoeffs[i], sourcePid)
						mutex.Unlock()
					}
					rst := pi.signTestComposed(cps, net, nil, leftEncode, kinshipCoeffs, &mutex, maskEncoded, sourcePid, idx, numIter)

					// should drop level before exchanging
					rst = crypto.DropLevelVec(cps, rst, 2)
					rst = crypto.CRescale(cps, rst)
					XCompResult.Result[idx] = rst[0]
					wg.Done()
				}(j)
			}
			wg.Wait()
		} else if pi.reveal == 3 {
			// reveal the raw kinship coefficients
			XCompResult.Result = append(XCompResult.Result, kinshipCoeffs[0])
		}
	} else {
		panic("reveal mode not supported")
	}
	return
}
