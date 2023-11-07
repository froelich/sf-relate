package internal

import (
	"io"
	"os"

	"github.com/sbinet/npyio/npz"
	"go.dedis.ch/onet/v3/log"
)

func ReadNPZVectors(rowIndexFile, columnIndexFile string) ([]int32, []int32) {
	// read rows indexes
	fRows, rRows := dataReader(rowIndexFile)
	defer fRows.Close()
	var mRows []int32
	// assume a single key
	name := rRows.Keys()[0]
	log.LLvl1("Reading file with rows index: ", name, rRows.Header(name).Descr.Shape)
	err := npz.Read(fRows, name, &mRows)
	if err != nil {
		log.Error(err)
	}
	// read column indexes
	fCols, rCols := dataReader(columnIndexFile)
	defer fCols.Close()
	var mCols []int32
	// assume a single key
	nameCols := rCols.Keys()[0]

	err = npz.Read(fCols, nameCols, &mCols)
	if err != nil {
		log.Error(err)
	}
	log.LLvl1("Reading file with cols index: ", nameCols, rCols.Header(nameCols).Descr.Shape)
	return mRows, mCols //[:1000]

}
func ReadNPZMatrix(matrixFile string) [][]float64 {
	var f *os.File
	var r *npz.Reader
	f, r = dataReader(matrixFile)
	defer f.Close()
	// assume single key
	name := r.Keys()[0]
	// read  and put into float matrix
	var m []int64
	err := npz.Read(f, name, &m)
	if err != nil {
		log.Error(err)
	}

	mfloat := make([][]float64, r.Header(name).Descr.Shape[0])
	for i := 0; i < len(mfloat); i++ {
		mfloat[i] = make([]float64, r.Header(name).Descr.Shape[1])
		for j := 0; j < len(mfloat[i]); j++ {
			mfloat[i][j] = float64(m[j*len(mfloat)+i])
		}
		if i < 10 {
			log.LLvl1(mfloat[i][:10])
		}
	}
	return mfloat
}

func ReadNPZMatrixSubset(rowIndex, columnIndex []int64, matrix [][]float64, startingIndex int, batchSize int) [][]float64 {
	matrixSubset := make([][]float64, batchSize)

	// vector for row selection
	rowIndex = rowIndex[startingIndex : startingIndex+batchSize]
	for i, v := range rowIndex {
		matrixSubset[i] = make([]float64, len(columnIndex))
		for j, w := range columnIndex {
			matrixSubset[i][j] = matrix[v][w]
		}
	}
	return matrixSubset
}
func ReadFilteredRows(file *os.File, rowIndex int, nrowsTotal int, ncols int, numSamples int, sampleFilter []int32) []float64 {

	if rowIndex >= nrowsTotal {
		panic("Row index out of range")
	}
	file.Seek(int64(rowIndex)*int64(ncols), 0)

	data := make([]byte, ncols)
	_, err := io.ReadFull(file, data)
	if err != nil {
		panic(err)
	}

	out := make([]float64, numSamples)

	idx := 0
	for i := range sampleFilter {
		//if sampleFilter[i] {
		if sampleFilter[i] == -1 {
			sampleFilter[i] = 1
		}
		out[idx] = float64(int8(data[sampleFilter[i]]))
		idx++
		//}
	}

	return out
}

func dataReader(filePath string) (*os.File, *npz.Reader) {
	f, err := os.Open(filePath)
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

	return f, r
}
