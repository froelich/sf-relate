cd notebooks
source param.sh
mkdir data/2party_"$n"/table/"$tid"/ -p
python3 step1_encode_and_bucketing.py | tee -a data/2party_"$n"/table/"$tid"/hash.log