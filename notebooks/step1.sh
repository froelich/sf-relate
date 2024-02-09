# REQUIRE: already run step0.sh and have the same parameters.
# Run on one machine
python3 step1_hashing.py -n 1601 -param trial -out trial/party1/table -hap trial/party1/haps -L 3
# # Run on the other machine
# python3 step1_hashing.py -n 1601 -param trial -out trial/party2/table -hap trial/party2/haps -L 3