# Number of parallel processes to use
# Set to 20 for the UKB dataset with 100K individual * 90K SNPs on the Google Cloud machine with 128 cores and 576GB memory.
# Should be set as large as possible to utilize all CPUs and memory. Exact value depends on the machine and dataset sizes. Users can provide reasonable parameters like 5 and retry with a smaller one if it fails due to memory constraints.
export PARA=1
export t="demo"
# Directory to the config files
export FOLDER="config/$t/"
# Ouput directories
export OUT_FOLDER="out/$t/"
