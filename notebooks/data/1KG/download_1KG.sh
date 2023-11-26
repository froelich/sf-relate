wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20231123.zip
unzip plink2_linux_avx2_20231123.zip
rm plink2_linux_avx2_20231123.zip

# fetch 1KG data from PLINK2 
wget https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
mv all_hg38.pgen.zst?dl=1 all_hg38.pgen.zst
wget "https://www.dropbox.com/scl/fi/id642dpdd858uy41og8qi/all_hg38_rs_noannot.pvar.zst?rlkey=sskyiyam1bsqweujjmxqv1h55&dl=1"
mv "all_hg38_rs_noannot.pvar.zst?rlkey=sskyiyam1bsqweujjmxqv1h55&dl=1" all_hg38.pvar.zst
wget https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1
mv hg38_corrected.psam?dl=1 all_hg38.psam

# decompress the genotype file
./plink2 --zd all_hg38.pgen.zst >all_hg38.pgen
rm all_hg38.pgen.zst