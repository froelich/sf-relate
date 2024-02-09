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
- Python (>=3.9.2) with [NumPy](https://numpy.org/install/), [joblib](https://joblib.readthedocs.io/en/stable/), [pandas](https://pandas.pydata.org/).

### Instructions
To install SF-Relate, clone the repository and try building as follows. It will automatically fetch dependencies.
```
git clone https://github.com/froelich/sf-relate.git
cd sf-relate
go get relativeMatch
go build
```

If `go build` produces an error, run commands suggested by Go and try again. If the build
finishes without any output, the package has been successfully configured.

## SF-Relate Parameters
### Command line arguments for [step 0](#step-0-sampling-shared-parameters) of SF-Relate 
This step samples shared hashing parameters and sketching parameters, respectively. 
There are two scripts for this step:
* `notebooks/step0a_sample_hash_randomness.py` samples shared random parameters that need to be the same across the two parties. We can let party 1 run this script and send the output files to party 2. It takes in the following parameters:
```
  -N N              size of hash tables (recommend: 64 * total number of
                    individuals on both parties)
  -pos POS          the directory storing the bp positions. Files are named
                    chr{i}.txt for i = 1..22. They contain the list of
                    physical positions of each SNPs on the haplotypes (one
                    number per line).
  -gmap GMAP        the directory storing the genetic maps. Files are named
                    chr{i}.gmap.gz (in gz format) for i = 1..22. The first
                    line of the file contains `pos\tchr\tcM`, and each following
                    line contains the bp location, the chromosome ID and the
                    corresponding genetic location (separated by tabs). One
                    can retrieve these files from [shapeit4](https://github.co
                    m/odelaneau/shapeit4/tree/master/maps) or other public
                    resources, but should be careful to make sure the correct
                    genome build locations is used.
  -maf MAF          the directory storing the maf files. Files are named
                    chr{i}.txt for i = 1..22. Each line in the file stores a
                    floating point number denoting the MAF of that base pair.
  -out OUT          output directory to store the parameters
```
It also takes in the following advanced parameters which come with default values. While the default values should work well on most datasets, users can modify them based on their needs, notably when specific IBD structures are known about the datasets.
```
  -enclen ENCLEN    the number of snps in each encoded split haplotype segment
                    (default: 80).
  -seglen SEGLEN    centi-Morgan length of each split haplotype segment
                    (default: 8.0)
  -steplen STEPLEN  centi-Morgan spacing between the beginning of each split
                    haplotype segment (default: 4.0)
  -k K              number of SNPs in each kSNP token for hashing (default: 8)
  -l L              number of hash tokens to construct every hash index
                    (default: 4)
  -maxL MAXL        max number of repetitive hashing; increase and retry if
                    table saturation is low (default: 6; should be larger than
                    the argument to step_1_hashing.py)
```
* `notebooks/step0b_sample_SNPs.py` samples a subset of SNPs to be used on [step2](#step-2-mhe). It takes in the following parameters:
```
-M M        number of SNPs on which KING is computed
-s S        subsampling rate (i.e. size(outputSNPs)/size(total))
-out OUT    output directory
```
### Command line arguments for [step 1](#step-1-hashing)
* `notebooks/step1_hashing.py` builds the hash tables using the shared parameters. It takes in the following parameters.
```
  -n N          number of individuals in the local dataset (i.e. the haplotype
                file is of dim 2n * M)
  -param PARAM  directory to the shared params output by step0a.
  -hap HAP      directory to the haplotype files, named chr$i.npy for i=1..22.
                Haplotypes are stored as a 2D numpy integer matrix of dim 2n *
                M, with different integers treated as different SNPs.
  -L L          number of repetitive hashing; increase and retry if table
                saturation is low; if not enough repetition hash keys are
                sampled in step0a (default maxL=6), redo step0a with a larger
                maxL.
  -out OUT      output table directory
```
Note that logs which include the table saturation is printed to `stdout`. The parameters stored in `-param` are the same as the output by `step0a_sample_hash_randomness.py` and need to be the same on both parties, but the other parameters (`n`, and `L`) need to be adjusted according to the local parties statistics (e.g. parties with a smaller `n` might need a larger `L` to saturate the local table).
### Example configuration files for [step 2](#step-2-mhe) (in `config/demo/` and `test_param.sh`).
*The current implementation for step2 assumes the two parties have the same number of individuals.*

For [step 2](#step-2-mhe), there are both global config parameters (`config/demo/configGlobal.toml`) shared by all parties and party-specific parameters (`config/demo/configLocal.Party[0-2].toml`).
#### The experimental parameters are stored in `configGlobal.toml`.
Users should modify them on their customized datasets.
```toml
# ================== EXPERIMENT CONFIGURATIONS ==================
# Number of individuals in the dataset (note that it was n in python code)
N = 1601
M = 145181 # total number of SNPs in the full panel

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
It also contains the following machine parameters, which should be updated.
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
####  The local parameters file (`config/demo/configLocal.Party[0-2].toml`) stores the following information for party i=0..2:
```toml
# This is the directory the binary genotype 
# file of all individuals. Each genotype is 
# encoded as a byte in individual-major order.
# So there are a total of n * M bytes.
simple_data_path = "notebooks/trial/party1/geno/all_chrs.bin"
# this is the output by the Python step1 script
row_index_file= "notebooks/trial/party1/table/ID_table.npz"
# this is the output by the Python step0b script.
column_index_file="notebooks/trial/sketched/SNPs.npz"
```
Note that party 0 is run on the same machine as party1, but there still needs to be a config file for party0 (with the same contents as the one for party 1) for the go code to work.
Also note that the number `M` is the number of SNPs on which the KING kinship estimator (see [Manichaikul et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025716/)) is computed. Each count (0/1/2) is encoded as individual bytes in the binary file.
#### Finally, `test_param.sh` defines some environmental variables that defines the config and output directories, along with the process-level of parallelization of the MHE.
```bash
# Number of parallel processes to use
# Set to 20 for the UKB dataset with 100K individual * 90K SNPs on the Google Cloud machine with 128 cores and 576GB memory.
# Should be set as large as possible to utilize all CPUs and memory. Exact value depends on the machine and dataset sizes. Users can provide reasonable parameters like 5 and retry with a smaller one if it fails due to memory constraints.
export PARA=1
export t="demo"
# Directory to the config files
export FOLDER="config/$t/"
# Output directories
export OUT_FOLDER="out/$t/"
```
## SF-Relate Usage
### Step 0: Sampling Shared Parameters
To run step 0, which samples shared parameters, follow the example instructions in `example.sh`:
```bash
python3 step0a_sample_hash_randomness.py -out trial -maf trial/maf -pos trial/pos -gmap trial/gmap -enclen 80 -N 204928 -k 8 -seglen 8 -steplen 4 -l 4
python3 step0b_sample_SNPs.py -M 145181 -s 0.7 -out trial/sketched
```
Be sure to share the parameters across the machines before executing [step 1](#step-1-hashing).
### Step 1: Hashing
Step 1 hashes the input samples into buckets, and needs to be executed on the two data parties locally as they handle sensitive data.
```bash
# REQUIRE: already run step0.sh and have the same parameters.
# Run on one machine
python3 step1_hashing.py -n 1601 -param trial -out trial/party1/table -hap trial/party1/haps -L 3
# # Run on the other machine
# python3 step1_hashing.py -n 1601 -param trial -out trial/party2/table -hap trial/party2/haps -L 3
```
#### Intermediate output: hash tables
For each party `i`, the list of buckets to be used in [step 2](#step-2-mhe) are stored at `{OUT}/ID_table.npz` (`trial/party{i}/table/` in the example).

### Step 2: MHE
Step 2 performs the kinship evaluation and accumulation under encryption over networks.
We provided a Makefile that reads in `test_param.sh` and triggers the corresponding parties.
```bash
# for party 1
make party1 -j2
# for party 2
make party2
```
The two `make` commands should be executed on each of the two data-contributing parties (`PID=1` and `PID=2`) separately.
There is also a placeholder party 0 (`PID=0`) that is executed on party 1, and the first job `party1` spawns 2 go jobs for this.
The exact `go` commands can be found in `[X,Y,Z]_local.sh`.

### Output
The parties only get their corresponding outputs.
Once the script finishes, it stores the output at `$OUT_FOLDER` (specified in `test_param.sh`), with the following files
- For `mode == 0`, the outputs are `[0-floor(n/8192)]_party{i}.csv`, storing the indicator for each local sample on party `{i}`, specifying whether they have a 3rd degree relative. 
- For `mode == 1`, the outputs are `degree{d}[0-floor(n/8192)]_party{i}.csv`, storing the indicator for each local sample on party `{i}`, specifying whether the max kinship passes the d-th relatedness threshold.
- For `mode == 2`, the outputs are `kinship{d}[0-floor(n/8192)]_party{i}.csv`, storing the indicator for each local sample on party `{i}`, specifying whether the max kinship passes the d-th customized kinship threshold.

For all those files, the header line contains `ID,Detected`, and each subsequent line contains two numbers separated by `,`, representing the sample's location in the file (0-indexed) and the decrypted result (a complex number).
When the second number is very close to 0, it signals that the kinship passes the detection threshold.
- Moreover, `$OUT_FOLDER/[X,Y,Z]/test.txt` stores the log of each party's execution.

## Contact for Questions
Matthew Man-Hou Hong, matthong@mit.edu; 
David Froelicher, dfroelic@mit.edu; 
Hoon Cho, hoon.cho@yale.edu