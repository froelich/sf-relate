package internal

import (
	"encoding/csv"
	"errors"
	"fmt"
	"go.dedis.ch/onet/v3/log"
	"io"
	"os"
	"strconv"
)

// CSV is a File splitting at the given Side
type CSV struct {
	Path      string
	Separator rune
	Split     Side
	Header    bool
}

var _ Loader = CSV{}

// Load loads the dataset found in the File
func (f CSV) Load(onlyInfo bool) ([][]float64, error) {
	matrix, err := f.loadMatrix(onlyInfo)
	if err != nil {
		return nil, fmt.Errorf("loading matrix: %v", err)
	}

	return matrix, err
}

func (f CSV) loadMatrix(onlyInfo bool) ([][]float64, error) {
	file, err := os.Open(f.Path)
	if err != nil {
		return nil, fmt.Errorf("opening file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = f.Separator
	reader.TrimLeadingSpace = true

	readHeader := false
	var records [][]float64
	for {
		record, err := reader.Read()
		if record == nil && errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading record: %w", err)
		}

		if f.Header && !readHeader {
			readHeader = true
			continue
		}

		line := make([]float64, len(record))
		for i, v := range record {
			parsed, err := strconv.ParseFloat(v, 64)
			if err != nil {
				return nil, fmt.Errorf("parsing as float at line %v, element %v: %w", len(records), i, err)
			}
			line[i] = parsed
		}

		records = append(records, line)

		if onlyInfo {
			return records, nil
		}
	}

	return records, nil
}

func ReadCsvFile(filePath string) [][]string {
	f, err := os.Open(filePath)
	if err != nil {
		log.Fatal("Unable to read input file "+filePath, err)
	}
	defer f.Close()

	csvReader := csv.NewReader(f)
	records, err := csvReader.ReadAll()
	if err != nil {
		log.Fatal("Unable to parse file as CSV for "+filePath, err)
	}

	return records
}
