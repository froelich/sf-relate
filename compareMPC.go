package relativeMatch

import (
	"math"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

// ComparisonResultWithOutDec contains the decrypted result
type ComparisonResultWithOutDec struct {
	Result              [][]float64
	OtherIndexEncrypted []float64
	NbrRepeat           int
	ResultControlled    [][]float64
	IndexControlled     []float64
}

func reportTime(timeStart time.Time, prot *ProtocolInfo, part_name string) {
	log.LLvl1(timeStart, "========= time at "+part_name, time.Since(timeStart))
}

func (pi *ProtocolInfo) compareMPC(Input ComparisonDataLocal, role string, net int) ComparisonResultMPC {
	nbrBuckets := int(math.Ceil(float64(pi.batchLength) / float64(pi.bucketSize)))
	rType := pi.basicProt.MpcObj[0].GetRType()

	squareRepLoc := make([]float64, 0)
	hetRepLoc := make([]float64, 0)
	var inputT *mat.Dense

	nbrOfComparisons := pi.batchLength * pi.bucketSize

	// initialize values for MPC
	distanceFloat := make([]float64, nbrOfComparisons)
	distance := mpc_core.FloatToRVec(rType, distanceFloat, pi.basicProt.MpcObj[0].GetFracBits())
	locSquareSS := mpc_core.FloatToRVec(rType, distanceFloat, pi.basicProt.MpcObj[0].GetFracBits())
	otherSquareSS := mpc_core.FloatToRVec(rType, distanceFloat, pi.basicProt.MpcObj[0].GetFracBits())
	timeStart := time.Now()
	reportTime(timeStart, pi, " start MPC: ")

	// prepare local vectors
	if role == "comparator" {
		log.LLvl1("Comparator prepares local vectors")
		hetRepLoc, _, squareRepLoc, _ = prepareXVectors(Input, pi.bucketSize, nbrBuckets, pi.scaleDown)
		locSquareSS = mpc_core.FloatToRVec(rType, squareRepLoc, pi.basicProt.MpcObj[0].GetFracBits())
	} else if role == "comparee" {
		log.LLvl1("Comparee prepares local vectors")
		// bucket size = 1
		squareRepLoc = append(squareRepLoc, Input.XSquare...)
		hetRepLoc = append(hetRepLoc, Input.Xhet...)
		inputT = mat.DenseCopyOf(Input.X.T())
		otherSquareSS = mpc_core.FloatToRVec(rType, squareRepLoc, pi.basicProt.MpcObj[0].GetFracBits())
	} else {
		log.LLvl1("Helper prepares local vectors")
		squareRepLoc = make([]float64, nbrOfComparisons)
		hetRepLoc = make([]float64, nbrOfComparisons)
	}
	// optimized matrix operation version
	log.LLvl1("The input matrix is ", nbrOfComparisons, " by ", pi.numberOfColumnsTest)
	locMat := mpc_core.InitRMat(pi.basicProt.MpcObj[0].GetRType().Zero(), pi.numberOfColumnsTest, nbrOfComparisons)
	otherMat := mpc_core.InitRMat(pi.basicProt.MpcObj[0].GetRType().Zero(), pi.numberOfColumnsTest, nbrOfComparisons)

	// prepare matrix
	for j := 0; j < pi.numberOfColumnsTest; j = j + 1 {
		locColumn := make([]float64, 0)
		if role == "comparator" {
			XelemRepeat := make([]float64, pi.bucketSize)

			for bIndex := 0; bIndex < nbrBuckets; bIndex++ {
				XrowIndex := bIndex
				Xelem := Input.X.At(XrowIndex, j) * (-2)
				XelemRepeat[0] = Xelem
				locColumn = append(locColumn, XelemRepeat...)
			}
		} else if role == "comparee" {
			locColumn = append(locColumn, inputT.RawRowView(j)...)
		}

		locElemRepeatColSS := mpc_core.FloatToRVec(rType, locColumn, pi.basicProt.MpcObj[0].GetFracBits())
		otherElemRepeatColSS := mpc_core.InitRVec(pi.basicProt.MpcObj[0].GetRType().Zero(), nbrOfComparisons)
		// optimized matrix operation version
		if role == "comparator" || role == "comparee" {
			locMat[j] = locElemRepeatColSS
			otherMat[j] = otherElemRepeatColSS
		}
	}

	sending := role == "comparator"
	if role == "comparator" {
		locElemRepeatColSSR, locElemRepeatColSSM := pi.basicProt.MpcObj[net].BeaverPartitionMat(locMat)
		otherElemRepeatColSSR, otherElemRepeatColSSM := pi.basicProt.MpcObj[net].BeaverPartitionMat(otherMat)
		distance = pi.basicProt.MpcObj[net].BeaverMultElemMat(locElemRepeatColSSR, locElemRepeatColSSM, otherElemRepeatColSSR, otherElemRepeatColSSM, sending).Sum(1)
	} else if role == "comparee" {
		locElemRepeatColSSR, locElemRepeatColSSM := pi.basicProt.MpcObj[net].BeaverPartitionMat(otherMat)
		otherElemRepeatColSSR, otherElemRepeatColSSM := pi.basicProt.MpcObj[net].BeaverPartitionMat(locMat)
		distance = pi.basicProt.MpcObj[net].BeaverMultElemMat(otherElemRepeatColSSR, otherElemRepeatColSSM, locElemRepeatColSSR, locElemRepeatColSSM, sending).Sum(1)
	} else {
		locElemRepeatColSSR, locElemRepeatColSSM := pi.basicProt.MpcObj[net].BeaverPartitionMat(locMat)
		otherElemRepeatColSSR, otherElemRepeatColSSM := pi.basicProt.MpcObj[net].BeaverPartitionMat(otherMat)
		distance = pi.basicProt.MpcObj[net].BeaverMultElemMat(locElemRepeatColSSR, locElemRepeatColSSM, otherElemRepeatColSSR, otherElemRepeatColSSM, sending).Sum(1)
	}
	distance = pi.basicProt.MpcObj[net].BeaverReconstructVec(distance)
	locSquareSS.Add(otherSquareSS)
	distance = pi.basicProt.MpcObj[net].TruncVec(distance, pi.basicProt.MpcObj[0].GetDataBits(), pi.basicProt.MpcObj[0].GetFracBits())
	distance.Add(locSquareSS)

	timeStart = time.Now()
	// sign test
	locHet := mpc_core.InitRVec(rType.Zero(), nbrOfComparisons)
	otherHet := mpc_core.InitRVec(rType.Zero(), nbrOfComparisons)
	if role == "comparator" {
		locHet = mpc_core.FloatToRVec(rType, hetRepLoc, pi.basicProt.MpcObj[0].GetFracBits())
	} else if role == "comparee" {
		otherHet = mpc_core.FloatToRVec(rType, hetRepLoc, pi.basicProt.MpcObj[0].GetFracBits())
	}

	signDiffHet := pi.basicProt.MpcObj[net].LessThan(locHet, otherHet, true)
	signDiffHetFlip := pi.basicProt.MpcObj[net].FlipBit(signDiffHet)

	// minimum (choose between hetx or hety)
	xHetFilt := pi.basicProt.MpcObj[net].SSMultElemVec(signDiffHet, locHet)
	yHetFilt := pi.basicProt.MpcObj[net].SSMultElemVec(signDiffHetFlip, otherHet)
	xHetFilt.Add(yHetFilt)

	rst := pi.basicProt.MpcObj[net].LessThan(distance, xHetFilt, true)
	// to facilitate decryption
	rst.MulScalar(pi.basicProt.MpcObj[0].GetRType().FromFloat64(1.0, pi.basicProt.MpcObj[0].GetFracBits()))
	if role == "comparator" || role == "comparee" {
		rstRevealed := pi.basicProt.MpcObj[net].RevealSymVec(rst)
		rsts := rstRevealed.ToFloat(pi.basicProt.MpcObj[0].GetFracBits())
		reportTime(timeStart, pi, " ended MPC: ")
		return ComparisonResultMPC{Result: rsts}

	}

	return ComparisonResultMPC{}

}
