# Keep a smaller set of variants
## remove rare variants 
./plink2 -pfile vzs all_hg38 --mind 0.1 --geno 0.1 --maf 0.01 --max-maf 0.99 --chr 1-22  --min-alleles 2 --max-alleles 2 --make-pgen vzs --out maf_0.01_all_hg
## thin data
./plink2 -pfile vzs maf_0.01_all_hg --thin-count 1000000 --make-pgen vzs --out thinned
## exclude HLD region
./plink2 -pfile vzs thinned --exclude bed1 exclude_HLD_b38.bed --king-cutoff 0.044 --make-pgen vzs --out unrelated

# Select a set of SNPs for KING inference
## LD-prune for PCA
./plink2 -pfile vzs unrelated --rm-dup force-first --indep-pairwise 1000 80 0.1 --out LD_filt
./plink2 -pfile vzs unrelated --extract LD_filt.prune.in --make-pgen vzs --out LD_pruned
## PCA
./plink2 -pfile vzs LD_pruned --pca biallelic-var-wts 8 --out pcavar
## Run Jupyter notebook for selection
python3 filterForKING.py
## Extract the selected SNPs for KING
./plink2 -pfile vzs thinned --extract snpsForKING.txt --make-pgen vzs --out snpsForKING
## compute the inferred relateives
./plink2 -pfile vzs snpsForKING --make-king-table --king-table-filter 0.044 --out king

# Decompress the variant file for later use
./plink2 -zd thinned.pvar.zst > thinned.pvar
./plink2 -zd snpsForKING.pvar.zst > snpsForKING.pvar