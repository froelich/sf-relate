cd notebooks
cd data/1KG
bash download_1KG.sh
bash filt.sh
cd ../..
python3 step0_split_1KG_to_tests.py
mkdir data/maps
cd data/maps
wget https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz
tar -xzf genetic_maps.b38.tar.gz