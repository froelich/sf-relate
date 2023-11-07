import os
import subprocess
import numpy as np
from bgen_reader import open_bgen
from my_util import read_snp_qc, run_command

for id in range(1, 23):
    maf_fname = f'data/maf/chr{id}.txt'
    # only keep lines that are not imputed (i.e. last entry = 1 in the row)
    print(f"filtering mafs on chr {id} ... ")
    open(maf_fname,'w').writelines(line for line in open(f"data/maf/ukb_mfi_chr{id}_v3.txt") if line.split('\t')[-1] == '1\n')
    print(f"========== done =========")

# prepare SNP QC information
# use the UKB official SNP QC file to output a list of filter for subsets of SNPs to be used in the KING computation
snp_qc = read_snp_qc()
for chromo in range(1, 23):
    print(f"working on {chromo}")
    snp_qc_now = snp_qc.get_group(chromo)
    with open_bgen(f"data/raw/chr{chromo}.bgen", verbose=False) as bgen:
        positions = bgen.positions.copy()
        ids = bgen.rsids.copy()

    np.save(f"data/maf/chr{chromo}", snp_qc_now[np.isin(snp_qc_now['position'], positions)]['in_Relatedness'].to_numpy().astype(np.bool8))