package internal

func KingDegree(val float64) int {
	threshValue := []float64{1.8232, 1.6464, 1.2928, 0.5857}
	if val < threshValue[3] {
		return 0
	} else if val < threshValue[2] {
		return 1
	} else if val < threshValue[1] {
		return 2
	} else if val < threshValue[0] {
		return 3
	} else {
		return 4
	}
}

func ComputeKingClear(rowLocal, rowOther []float64, scale float64) (float64, float64, float64) {
	kingdistance := 0.0
	hetLocal := 0.0
	hetOther := 0.0
	xy := 0.0
	for c := 0; c < len(rowLocal); c++ {
		kingdistance = kingdistance + (rowLocal[c]-rowOther[c])*(rowLocal[c]-rowOther[c])
		xy += -2.0 * rowLocal[c] * rowOther[c]
		if rowLocal[c] == 1 {
			hetLocal = hetLocal + 1
		}
		if rowOther[c] == 1 {
			hetOther = hetOther + 1
		}
	}
	min := hetLocal
	if hetOther < hetLocal {
		min = hetOther
	}
	kinship := 1.0/2.0 - (1.0/4.0)*(kingdistance/min)
	return kinship, xy, kingdistance / float64(len(rowLocal)) * scale
}
