# SF-Relate

Software for secure and federated genetic relatives detection, as described in:

**Secure Discovery of Genetic Relatives across Large-Scale and Distributed Genomic Datasets**\
Matthew Man-Hou Hong, David Froelicher, Ricky Magner, Victoria Popic, Bonnie Berger, and Hyunghoon Cho,
Under review, 2023

This repository contains a set of scripts for generating test cases for testing sf-relate.
- In the branch `1KG`, we demonstrate using the software to detect related samples in the publicly available [1000 Genomes Project](https://pubmed.ncbi.nlm.nih.gov/36055201/).
- In the branch `UKB`, we store scripts specialized to generate the inputs from the __access-limited__ [UK-Biobank](https://www.ukbiobank.ac.uk/).
- For usages on other datasets, refer to [SF-Relate Usage](#sf-relate-usage)

## Installation

### Dependencies

SF-Relate requires that `go` and `python3` are available in the exec path in shell. Here are the links for installation:

- [Go](https://go.dev/doc/install) (>=1.18.3)
- Python (>=3.9.2) with [NumPy](https://numpy.org/install/), [joblib](https://joblib.readthedocs.io/en/stable/), [pandas](https://pandas.pydata.org/).

### Install SF-Relate

To install SF-Relate, clone the repository and try building as follows. It will automatically fetch dependencies.
```
git clone https://github.com/froelich/sf-relate.git
cd sf-relate
go get relativeMatch
go build
```

If `go build` produces an error, run commands suggested by Go and try again. If the build
finishes without any output, the package has been successfully configured.

## Running an automatic pipeline for on [1000 Genomes](https://pubmed.ncbi.nlm.nih.gov/36055201/) Data
We provide a series of bash scripts for building and testing the pipeline using 1000 Genome samples. [Extra dependencies](#extra-dependencies) are required to process 1000 Genomes data.

To run them, first make sure all dependencies are installed, 
then execute the following.
```
bash 0_prepare_1KG.sh
bash 1_hashing.sh
bash 2_sketch.sh
bash 3_run_MHE.sh
bash 4_verify_output.sh
```

## Preparation of Test Data from [1000 Genomes](https://pubmed.ncbi.nlm.nih.gov/36055201/)
The following describes how to generate example test data based on 1000 Genomes phase 3, phased data hosted on [PLINK2](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg).

The generated test data are split between two parties. Party `i`'s local data is stored in
the `notebooks/data/2party_{n}/party{i}`,
where `n` is the number of samples on each party, and `i` is 0-indexed.
By default, we use `n=1601` (i.e. all samples from 1000 Genomes).
The file `notebooks/param.sh` stores the default parameters used in the generation process and can be modified when needed.

### Extra Dependencies 
- [PLINK 2, build 2023.11.23](https://www.cog-genomics.org/plink/2.0/) is required for manipulating the 1000 Genome data.
- [Pgenlib v0.9.1](https://pypi.org/project/Pgenlib/) is required for reading the PLINK2 `.pgen` files in Python.
- Linux commandline tools `wget` and `unzip` (install via `sudo apt-get install wget unzip` if you have admin rights).

### Test Case Preparation Detail
The script `0_prepare_1KG.sh` runs the following:
- `download_1KG.sh` downloads the 1KG dataset from [PLINK2](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg) using `wget`.
Note that this script also downloads PLINK2 to fetch and manipulate `.pgen` files.
- `filt.sh` performs a PCA to select a subset of SNPs to be used in the KING inference, following [UK Biobank's Relatedness Inference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4380465/). It also thins the dataset down to 1M loci (enough for relatedness inference).
- `step0_split_1KG_to_tests.py` randomly generate a split dataset in `notebooks/data/2party_1601/party{i}` for subsequent tests.

## SF-Relate Usage
To run SF-Relate on more than 2 parties, run it between every pair of parties.

### Step 1 --- Hashing and bucketing
The script `1_hashing.sh` reads the parameters from `param.sh` and runs `notebooks/step1_encode_and_bucketing.py`, which creates the output directory `notebooks/data/2party_1601/table/mode4cM4len160k8/`.
Note that the path may change if parameters are changed.
Logs are stored at `hash.log` in the same directory.
For party `i`, the list of buckets are stored at `tables/party{i + 1}/ID_table.npz`

For real-world usages when datasets are on two machines, the Python scripts should be updated accordingly.

#### Input Format
In the following files, `mf` in the file names (e.g. 0.01) signifies that SNPs with minor allele frequencies (MAF) less than mf are removed.

The genetic map files need to be placed at
- `notebooks/data/maps/chr[1-22].b(37|38).gmap.gz`. One can retrieve these files from [shapeit4](https://github.com/odelaneau/shapeit4/tree/master/maps) or other public resources.
- `notebooks/data/pos/chr[1-22].txt` contains the list of physical positions of each SNPs on the haplotypes.
Note that the build 37/38 are used in UK Biobank data/1000 Genomes data, repectively. One needs to ensure that these physical positions match the build version used in the genetic maps for a customized case.

For a test dataset of `n` individuals on each side, by default, it should be placed at `notebooks/data/2party_{n}`, with the following contents.
- `party{i}/haps/chr[1-22].npy` stores the phased haplotypes of party `i`, each being a numpy matrix of `2n` rows and `m` columns. 
Haplotypes (rows) are encoded as vector of integers (`np.int8` or `np.int16`), and different bytes are considered different variants by default.
- `geno/party{i+1}/all_chrs_maf{mf}.bin` stores the genotype count matrix of size `n` by `M`, where `M <= m` is the number of SNPs on which KING (see [Manichaikul et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025716/)) is computed. Each count (0/1/2) is encoded as individual bytes in the binary file, in row major. 

Other files generated by the test preparation scripts (but not needed for a real distributed use case) are:
- `party{i}/maf/{mf}/chr[1-22]_numsnps.txt` stores the number of SNPs in each haplotype.
- `ground_truths/all_king_maf{mf}.npy` stores all pairwise KING coefficients between the samples, computed on the subset of `M` SNPs.
- `ground_truths/KING.dat` stores all related pairs across the parties, where each row contains `P0	ID0	P1	ID1	Kinship	deg` with the following meaning
    - `P0` and `P1` specifies where the samples are located (party 0 or 1)
    - `ID0` and `ID1` specifies which row of the genotype matrix the sample is located on.
    - `Kinship` is the value of the KING coefficient.
    - `deg` is relatedness degree, determined from `Kinship` based on the recommended thresholds in [Manichaikul et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025716/).

### Step 2 --- MHE
#### Subsampling SNPs for Faster KING Computations
The script `2_sketching.sh` reads in the number of SNPs `M` from the matrix `geno/party{i+1}/all_chrs_maf{mf}.bin` and run `notebooks/step2_subsample_SNPs.py`, which randomly subsample the loci of SNPs to be used and save them in `sketched/`.

#### Setting the Configuration
Example configuration files are provided in `config/demo/`. 
There are both global config parameters shared by all parties and party-specific parameters (that specifies the directories of the input data).
The input directories are set to the tables generated with the default parameters in `notebooks/param.sh`.
In order to run MHE on other tables, update the directories in the local configurations accordingly.

#### Running the MHE 
The example `Makefile` contains the example `go` commands for running the __MHE-Phase 1__ and __MHE-Phase 2__ in SF-Relate.

To run the two machines locally, use `make -j3`.  This script spawns 3 processes on the same machine---one for each of the two data-contributing parties (`PID=1` and `PID=2`)and the third for the an auxiliary party that helps synchronize the computation(`PID=0`). 

In practice, each party runs their process on their own machine and provides the IP addresses of other parties in the configuration for network communication. The auxiliary party (`PID=0`) can be run on the same machine as party 1. In other words, using the same `Makefile`, on two separate machines, party 1 should run `make X Z -j2` and party 2 should run `make Y`, with the ports in `config/demo` correctly adjusted.

#### Output
Once SF-Relate finishes, it stores its output at `out/demo/`, with the following files
- `[0-floor(n/8192)]_party{i}.csv` stores the indicator for each local sample on party {i}, specifying whether they have a relative. The order in which the boolean value appears correspond to the order in which the haplotype/genotypes appear in the input.
- `[X,Y,Z]/test.txt` stores the log of each party's execution.

### Optional --- Verifying the output
The optional python script `python3 step4_verify_output.py` reports the recall and precision of the test and print them on `stdout`.

## References
- The 1000 Genomes Project Consortium. A global reference for human genetic variation. _Nature_ __526__, 68–74 (2015).
- Byrska-Bishop M, Evani US, Zhao X, et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. _Cell_ __2022;185(18):3426-3440.e19__. doi:10.1016/j.cell.2022.08.004

## Contact for Questions
Matthew Man-Hou Hong, matthong@mit.edu; 
David Froelicher, dfroelic@mit.edu; 
Hoon Cho, hoon.cho@yale.edu
