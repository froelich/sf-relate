# SF-Relate

Software for secure and federated genetic relatives detection, as described in:

**Secure Discovery of Genetic Relatives across Large-Scale and Distributed Genomic Datasets**\
Matthew Man-Hou Hong, David Froelicher, Ricky Magner, Victoria Popic, Bonnie Berger, and Hyunghoon Cho,
Under review, 2023

This repository contains a set of scripts for generating test cases for testing sf-relate.
- In the branch `sfkit` (default), we provide the software for use on a single party on a machine in a federated study.
- In the branch `1KG`, we demonstrate the software in detection of related samples in the __publicly available__ [1000 Genomes Project](https://pubmed.ncbi.nlm.nih.gov/36055201/). (See [An Automatic Pipeline For Testing SF-Relate](#an-automatic-pipeline-for-1000-genomes-data))
- In the branch `UKB`, we provide scripts for generating the test cases based on the __access-limited__ [UK-Biobank](https://www.ukbiobank.ac.uk/). The UK Biobank files need to be stored at the correct paths for these scripts.

READMEs on the corresponding branches detail the different usages.

## Install SF-Relate

### Dependencies

SF-Relate requires that `go` and `python3` are available in the exec path in shell. Here are the links for installation:

- [Go](https://go.dev/doc/install) (>=1.18.3)
- Python (>=3.9.2) with [NumPy](https://numpy.org/install/), [joblib](https://joblib.readthedocs.io/en/stable/), [pandas](https://pandas.pydata.org/), [tomlkit](https://pypi.org/project/tomlkit/0.4.6/) and [Pgenlib](https://pypi.org/project/Pgenlib/).

### Instructions
To install SF-Relate, clone the repository and try building as follows. It will automatically fetch dependencies.
```
git clone https://github.com/froelich/sf-relate.git
cd sf-relate
go get relativeMatch
go build
go test -c -o goParty
```
If `go build` produces an error, run commands suggested by Go and try again. If the build
finishes without any output, the package has been successfully configured.

## Usage
### Input data
Input to this workflow consists of the following files:
- `all_chrs.[pgen|psam|pvar]` - Phased haplotype and metadata files encoding sample and variant information in the standard [PLINK2 genotype format](https://www.cog-genomics.org/plink/2.0/input#pgen).
- `chr[1-22].gmap.gz` - Gzipped genetic maps. The first line of the file contains `pos\tchr\tcM`, and each following line contains the bp location, the chromosome ID and the corresponding genetic location (separated by tabs). One can retrieve these files from [shapeit4](https://github.com/odelaneau/shapeit4/tree/master/maps) or other public resources, but should be careful to make sure the genome build positions match the ones in `haps/chr[1-22].pvar`.
- `snpsForKING.txt` - This file lists the RSIDs, one on each line, on which the KING estimator is computed. The RSIDs should be a subset of the variants in the .pvar file.
- `chr[1-22].txt` - Allele frequency files. In the order of appearance in the `all_chrs.pvar` file, each line in the file stores a floating point number denoting the minor allele frequency of the base pair.

### Pipeline
#### Exchanging secret keys
In the example configuration folder `$FOLDER=config/demo` we have the following:
- `$FOLDER/keys/shared_key_$i_$j.bin` and `$FOLDER/keys/shared_key_global.bin` store example 32-byte secrets in a binary format (`i` and `j` denote the party's indices).
However, it is necessary to use secure key exchange protocols to acquire shared key materials to derive cryptographic keys in real deployments.
- `$FOLDER/keys/seed.bin` stores an example 16-byte random seed for locally generating shared parameters for hashing.
Our script does not use this for cryptographic purposes, but it derives the shared seed from the shared secrets and cache it in this file.

Specify the input in the [configurations files](#configurations) in `FOLDER`, before running any of the subsequent steps.
#### Preparing additional input files
We provided a helper script `notebooks/pgen_to_npy.py` that creates intermediate cache files for the pipeline. 
```bash
python3 notebooks/pgen_to_npy.py -PARTY $i -FOLDER $FOLDER
```
- It caches intermediate genotype and haplotype as `haps/chr$i.npy` and `geno/all_chrs.bin` files for faster Python and Go processing.
- It extracts the list of base pairs positions as `pos/chr$i.txt`.
#### Example script for running the entire pipeline
We provided a Makefile that reads the path to the configuration files defined in `test_param.sh` and triggers the corresponding party's execution.
```bash
# for party 1
make party1 -j2
# for party 2
make party2
```
The two `make` commands should be executed on each of the two data-contributing parties (`PID=1` and `PID=2`) separately.
There is also a placeholder party 0 (`PID=0`) that is executed on party 1, and the first job `party1` spawns 2 go jobs for this.
The exact `go` and `python` commands also detailed in the following, can be found in `[X,Y,Z]_local.sh`.

##### Step 0: Sampling Shared Randomness 
```bash
python3 step0_sample_shared_randomness.py -PARTY $i -FOLDER $FOLDER
```
It generates the shared parameters by applying a non-cryptographic random number generator to the seed `$FOLDER/keys/seed.bin` and saves them in `shared_param_dir` and `sketched_snps_dir`.
##### Step 1: Hashing
Both parties locally execute step 1 to hash the input samples into buckets.
They do it locally because this step handles sensitive data.
```bash
python3 step1_hashing.py -PARTY $i -FOLDER $FOLDER
```
For party `i`, the list of buckets used in [step 2](#step-2-mhe) are stored as `{hash_table_dir}/ID_table.npz`.

##### Step 2: MHE
Step 2 performs the kinship evaluation and accumulation under encryption over networks.
Use the following to run step 2
```bash
CONFIG=config/demo PARTY=1 ./goParty
```
Note that the binary executable re-runs step 0 to step 1.
##### Step 3: Post-process output
Step 3 locally clean up the decrypted results from step 2 into human-readable formats (see [output](#output)).
```bash
python3 notebooks/step3_post_process.py -PARTY $i -FOLDER $FOLDER
```

### Output
The parties only get their corresponding outputs.
Once the script finishes, it stores the output at `output_dir`:
- For `mode == 0`, the output is `boolean_party{i}.tsv`, storing the indicator for each local sample on party `i`, specifying whether they have a 3rd degree relative. 
- For `mode == 1`, the output is `degree_party{i}.tsv`, storing the closest degree for each local sample on party `i`.
- For `mode == 2`, the output is `bin_party{i}.tsv`, storing the index of the non-zero bin corresponding to the greatest kinship for each local sample on party `i`.
- For `mode == 3`, the output is `kinship_block{i}.tsv`, storing the computed distributed kinship computed.

For all those files, the header line contains `ID\tResult`, and each subsequent line contains two numbers separated by `\t`, containing the sample's IID (from `.psam`) and the decrypted result.
Moreover, `output_dir/[X,Y,Z]/test.txt` stores the log of each party's execution.
Finally, `output_dir/raw/` stores The raw contents of the decrypted ciphertexts.


## Configuration
The path to the configuration files is `$FOLDER`.
We show example configurations in `config/demo/`.
There are both global configuration parameters (`$FOLDER/configGlobal.toml`) shared by all parties and local parameters (`$FOLDER/configLocal.Party[0-2].toml`).

Note party 0 is run on the same machine as party1, but there still needs to be a config file for party0 (with the same contents as the one for party 1) for the go code to work.
#### Customizing the location of the configuration files.
The helper script `test_param.sh` defines where to look for configuration files.
```bash
export t="demo"
export FOLDER="config/$t/"
```

### Local Parameters (`configLocal.Party[0-2].toml`)
```toml
PARA = 1 # Number of parallel processes to use. Set to 20 for the UKB dataset with 100K individual * 90K SNPs on the Google Cloud machine with 128 cores and 576GB memory. Should be set as large as possible to utilize all CPUs and memory. Exact value depends on the machine and dataset sizes. Users can provide reasonable parameters like 5 and retry with a smaller one if it fails due to memory constraints.

# input directories
haps_dir = "notebooks/trial/party1/haps/" # containing all_chrs.[pgen|pvar|psam]
snp_list = "notebooks/trial/snps_king.txt"
gmap_dir = "notebooks/trial/gmap/" # containing chr[1-22].gmap.gz

# where to save intermediate outputs
geno_dir = "notebooks/trial/party1/geno/" # genotypes (reconstructed from haplotypes in haps_dir)
sketched_snps_dir = "notebooks/trial/sketched/" ## share this with the other party
shared_param_dir = "notebooks/trial/params/" ## share this with the other party
hash_table_dir = "notebooks/trial/party1/table/" ## local hash table directory
```

### Global Experimental Parameters (`configGlobal.toml`)
#### Configuring [step 0: Sampling shared randomness](#step-0-sampling-shared-randomness)
This step samples shared hashing parameters and sketching parameters, respectively. 
```toml
# =================== EXPERIMENT CONFIGURATIONS ==================

## ==================== STEP 0 (BASIC) ===========================
N = 204928 # size of hash tables (recommend: 64 * total number of individuals on both parties)
# Increase N to improve recall, 
# at the cost of inicreasing number of comparisons and thus runtime in step 2

## ==================== STEP 0 (ADVANCED) ========================
# The following advanced parameters come with default values. 
# While the default values should work well on most datasets, 
# users can modify them based on their needs, notably when specific IBD structures are known about the datasets.
enclen = 80 # the number of snps in each encoded split haplotype segment
seglen = 8.0 # centi-Morgan length of each split haplotype segment
steplen = 4.0 # centi-Morgan spacing between the beginning of each split haplotype segment
k = 8 # number of SNPs in each kSNP token for hashing
l = 4 # number of hash tokens to construct every hash index
maxL = 6 # max number of repetitive hashing; increase and retry if table saturation is low (should be larger than the argument to step_1_hashing.py)
s = 0.7 # subsampling rate (i.e. num(outputSNPs)/num(SNPs in .pvar))
```
#### Configuring [step 1: hashing](#step-1-hashing)
```toml
## ======================== STEP 1 ===============================
L = 3 # number of repetitive hashing; increase and retry if table saturation is low; if not enough repetition hash keys are sampled in step0 (default number of keys is maxL=6), redo step0 with a larger maxL.
```
Note that logs which include information about the table saturation is printed to `stdout`. 
The other parameters (`L`) need to be adjusted according to the local statistics (e.g. parties with a smaller number of individuals might need a larger `L` to saturate the local table).

#### Configuring [step 2: MHE](#step-2-mhe)
```toml
## ======================== STEP 2 ===============================
# select output modes
reveal = 0
# reveal = 1
# reveal = 2
# reveal = 3
# reveal = 0 is SF-Related's default (one indicator per individual per party) that computes whether the max kinship is larger than the degree 3 threshold
# reveal = 1 computes 4 indicators per individual per party, each meaning whether the max kinship is larger than the corresponding degree threshold in thres_value 
thresh_value= [1.8232, 1.6464, 1.2928, 0.5857] # correspond to degree 3, 2, 1, 0
# reveal = 2 computes 16 indicators per individual per party using the following thresholds (can be customized, but the runtime of phase 2 will scale accordingly.)
discretized_thresh= [2.0,1.9375,1.875,1.8125,1.75,1.6875,1.625,1.5625,1.5,1.375,1.25,1.125,1.0,0.75,0.5,0.25] # each sub-interval between the degrees are is split into 4
# reveal = 3 reveals all computed intermediate kinship and decrypt them.
```
#### Machine configurations.
```toml
# ================== MACHINE CONFIGURATIONS ==================
# num threads --- should be set to about 10 * num of cores
nbr_threads = 1280 

# Ports for listening to incoming connections from the other parties
# Party 0 & 1 are hosted on the same machine
[servers.party0]
ipaddr = "127.0.0.1" #should be local IP of GCP machine
ports  = {party1 = "5110", party2 = "7320"}

[servers.party1]
ipaddr = "127.0.0.1"
ports  = {party2 = "9210"}

[servers.party2]
ipaddr = "127.0.0.1"
ports  = {}
```

## Contact for Questions
Matthew Man-Hou Hong, matthong@mit.edu; 
David Froelicher, dfroelic@mit.edu; 
Hoon Cho, hoon.cho@yale.edu
