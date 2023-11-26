cd notebooks/
source param.sh
mkdir data/2party_"$n"/sketched/
python3 step2_subsample_SNPs.py | tee -a data/2party_"$n"/sketched/sketch.log