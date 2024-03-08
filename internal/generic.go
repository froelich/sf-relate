package internal

import (
	"fmt"
	"os"

	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

func LoadDataset(path string, separator rune, header bool) ([][]float64, error) {

	data := TabulatedMatrixLoader{
		Path: path,
	}

	dataset, err := data.Load(separator, header)
	if err != nil {
		return nil, err
	}
	return dataset, nil
}

type TabulatedMatrixLoader struct{ Path string }

// Load loads the dataset at the given Path
func (c TabulatedMatrixLoader) Load(separator rune, header bool) ([][]float64, error) {
	ds, err := CSV{
		Path:      c.Path,
		Separator: separator,
		Header:    header,
	}.Load(false)

	return ds, err
}

// Loader loads a DataSet
// onlyInfo is to get the header of a DataSet, usually the number of column
type Loader interface {
	Load(onlyInfo bool) ([][]float64, error)
}

func LoadMatrix(loader Loader) (*mat.Dense, error) {
	dataset, err := loader.Load(false)

	if err != nil {
		return nil, err
	}

	return mat.NewDense(len(dataset), len(dataset[0]), Flatten(dataset)), nil
}

func Flatten(array [][]float64) []float64 {
	ret := make([]float64, len(array)*len(array[0]))
	lineLen := len(array[0])

	for i, col := range array {
		for j, el := range col {
			ret[i*lineLen+j] = el
		}
	}

	return ret
}

func ToFloat(array []uint64) []float64 {
	ret := make([]float64, len(array))
	for i, col := range array {
		ret[i] = float64(col)
	}

	return ret
}

// Side is Left or Right
type Side bool

func ReadDataDim(filePath string) (int, int, int) {
	f, err := os.Open(filePath)
	if err != nil {
		log.Fatal("Unable to read input file "+filePath, err)
	}
	defer f.Close()

	// read three integers
	var n, M, m int
	_, err = fmt.Fscanf(f, "%d\n%d\n%d", &n, &M, &m)
	if err != nil {
		log.Fatal("Unable to read input file "+filePath, err)
	}

	return n, M, m
}
