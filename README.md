# SF-Relate

Software for secure and federated genetic relatives detection, as described in:

**Secure Discovery of Genetic Relatives across Large-Scale and Distributed Genomic Datasets**\
Matthew Man-Hou Hong, David Froelicher, Ricky Magner, Victoria Popic, Bonnie Berger, and Hyunghoon Cho,
Under review, 2023

This repository contains a set of scripts specialized to generate the inputs from the UK-Biobank dataset for testing the tool.
We are working on developing a simple demonstration of the software based on the (publicly available) 1000 Genomes Project for convenience.

## Installation

### Dependencies

SF-Relate requires that `go` and `python3` are available in the exec path in shell. Here are the links for installation:

- [Go](https://go.dev/doc/install) (>=1.18.3)
- Python (>=3.9.2) with [NumPy](https://numpy.org/install/), [joblib](https://joblib.readthedocs.io/en/stable/), [pandas](https://pandas.pydata.org/) and [bgen-reader](https://pypi.org/project/bgen-reader/).

### Install SF-Relate

To install SF-Relate, clone the repository and try building as follows. It will automatically fetch dependencies.
```
git clone https://github.com/froelich/sf-relate.git
cd sf-relate
go get github.com/froelich/sf-relate
go build
```

If `go build` produces an error, run commands suggested by Go and try again. If the build
finishes without any output, the package has been successfully configured.

## Preparation of Test Data from [1000 Genomes](https://pubmed.ncbi.nlm.nih.gov/36055201/)
The following describes how to generate example test data based on 1000 Genomes phase 3, phased data hosted on [PLINK2](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg).

The generated test data are split between two parties. Party `i`'s local data is stored in
the `notebooks/data/2party_{n}/party{i}`,
where `n` is the number of samples on each party, and `i` is 0-indexed.
By default, we use `n=1601` (i.e. all samples from 1000 Genomes).
The file `notebooks/param.sh` stores the default parameters used in the generation process and can be modified when needed.

### Extra Dependencies 
- [PLINK 2, build 2023.11.23](https://www.cog-genomics.org/plink/2.0/) is required for manipulating the 1000 Genome data.
- [Pgelib v0.9.1](https://pypi.org/project/Pgenlib/) is required for reading the PLINK2 `.pgen` files.

### Usage
```
bash 0_prepare_1KG.sh
```

## SF-Relate Usage
To run SF-Relate on more than 2 parties, run it between every pair of parties.

### Step 1 --- Hashing and bucketing
To run the bucketing step, run
```
bash 1_hashing.sh
```
With the default parameters, for party `i`, the list of samples in buckets are stored at `notebooks/data/2party_1601/table/mode4cM4len160k8/tables/party{i + 1}/ID_table.npz` (note the party's indices are shifted).
For real-case usages when datasets are on two machines, the directories in the scripts should be updated accordingly.

#### Input Format
In the following files, `mf` in the file names (e.g. 0.01) signifies that SNPs with minor allele frequencies (MAF) under mf are removed from UKB in the preparation.

By default, a generated test case is stored at `notebooks/data/2party_{n}`. All input files are under this directory.
- `party{i}/haps/chr[1-22].npy` stores the phased haplotypes of party `i`, each being a numpy matrix of `2n` rows. 
Haplotypes (rows) are encoded as vector of bytes, and different bytes are considered different variants by default.
- `geno/party{i+1}/all_chrs_maf{mf}.bin` stores the genotype count matrix of size `n` by `M`, where `M` is the number of SNPs on which KING (see [Manichaikul et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025716/)) is computed. Each count (0/1/2) is encoded as bytes in the binary file, in row major. 

Additionally, the genetic map files need to be placed at
- `notebooks/data/maps/chr[1-22].b37.gmap.gz`. One can retrieve these files from [shapeit4](https://github.com/odelaneau/shapeit4/tree/master/maps) or other public resources.
- `notebooks/data/pos/chr[1-22].txt` contains the list of physical positions of each SNPs on the haplotypes.
Note that the build 37 is used in UK Biobank data, and one needs to ensure that these physical positions match the build version used in the genetic maps.

Other files generated from the [preparation step](#preparation-of-test-data-from-1000-genomes) are:
- `party{i}/maf/{mf}/chr[1-22]_numsnps.txt` stores the number of SNPs in each haplotype.
- `ground_truths/all_king_maf{mf}.npy` stores all pairwise KING coefficients between the samples, computed on the subset of `M` SNPs.
- `ground_truths/KING.dat` stores all related pairs across the parties, where each row contains `P0	ID0	P1	ID1	Kinship	deg` with the following meaning
    - `P0` and `P1` specifies where the samples are located (party 0 or 1)
    - `ID0` and `ID1` specifies which row of the genotype matrix the sample is located on.
    - `Kinship` is the value of the KING coefficient.
    - `deg` is relatedness degree, based on the recommended thresholds in [Manichaikul et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025716/).

### Step 2 --- MHE
#### Subsampling SNPs for Faster KING Computations
Use the following to subsample SNPs for Step 2 in SF-Relate.
```
bash 2_sketching.sh
```

#### Setting the Configuration

Example configuration files are provided in `config/demo/`. 
There are both global config parameters shared by all parties and party-specific parameters (that specifies the directories of the input data).
The input directories are set to the tables generated with the default parameters in `notebooks/param.sh`.
In order to run MHE on other tables, update the directories in the local configurations accordingly.

#### Running the MHE 
Use the following to run the two machines locally.
```
bash 3_run_MHE.sh
```

In more detail, an example make script `Makefile` shows how to run the online phases, MHE-Phase 1 and MHE-Phase 2 in SF-Relate.
To run all parties, use `make -j3`.
Otherwise, on two separate machines, party 1 should run `make X Z -j2` and party 2 should run `make Y`, with the ports in `config/demo` correctly adjusted.
The script spawns 3 processes on the same machine---one for each of the two data-contributing parties (`PID=1` and `PID=2`)and the third for the an auxiliary party that helps synchronize the computation(`PID=0`). 
In practice, each party runs their process on their own machine and provides the IP addresses of other parties in the configuration for network communication. The auxiliary party (`PID=0`) can be run on the same machine as party 1.

#### Output
Once SF-Relate finishes, it stores its output at `out/demo/`, with the following files
- `[0-floor(n/8192)]_party{i}.csv` stores the indicator for each local sample on party {i}, specifying whether they have a relative. The order in which the boolean value appears correspond to the order in which the haplotype/genotypes appear in the input.
- `[X,Y,Z]/test.txt` stores the log of each party's execution.

### Step 4 --- Verifying the output
Use the following script to report the recall and precision of the test.
```
bash 4_verify_output.sh
```

## Contact for Questions
Matthew Man-Hou Hong, matthong@mit.edu; 
David Froelicher, dfroelic@mit.edu; 
Hoon Cho, hoon.cho@yale.edu
