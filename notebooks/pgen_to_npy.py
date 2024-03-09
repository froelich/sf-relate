import pgenlib
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
import tomlkit
from my_util import *
import os

parser = argparse.ArgumentParser(description='Sample Hashing Randomness')
parser.add_argument('-PARTY', type=int, help='party id', required=True)
parser.add_argument('-FOLDER', type=str, help='path to the configuration folder', required=True)
args = parser.parse_args()
party_id = args.PARTY
FOLDER = args.FOLDER

args_local = tomlkit.load(open(f"{FOLDER}/configLocal.Party{party_id}.toml"))
args_global = tomlkit.load(open(f"{FOLDER}/configGlobal.toml"))

# mkdir all intermediate paths if not exist
os.makedirs(args_local['geno_dir'], exist_ok=True)
os.makedirs(args_local['pos_dir'], exist_ok=True)
os.makedirs(args_local['sketched_snps_dir'], exist_ok=True)
os.makedirs(args_local['shared_param_dir'], exist_ok=True)
os.makedirs(args_local['hash_table_dir'], exist_ok=True)


# input: (1) PLINK2 (.pgen, .pvar, .psam) files, (2) SNPs for KING (.txt with rsids)
fprefix = args_local['haps_dir'] + "all_chrs"
pgen_file = fprefix + ".pgen" # need to be phased

psam, snps, var_info = read_pgen_metadata(args_local)

def read_pgen(data_path):
    reader = pgenlib.PgenReader(data_path.encode("utf-8"))
    n = 2 * reader.get_raw_sample_ct()
    m = reader.get_variant_ct()
    phase_present = np.zeros((m, n), dtype=np.int8)
    buf = np.zeros((m, n), dtype=np.int32)
    reader.read_alleles_and_phasepresent_list(np.arange(0, m, dtype=np.uint32), buf, phase_present)
    return n, m, buf.transpose().astype(np.int8)

n, m, buf = read_pgen(pgen_file)

# recode pgen to numpy files
# and also output the list of pos
# and the dimension of the matrix

# split by chromosome
for i in range(1, 23):
    snp_range = var_info['#CHROM'] == i
    var_info_chr = var_info[snp_range]
    np.savetxt(f"{args_local['pos_dir']}/chr{i}.txt", var_info_chr.loc[snp_range, 'POS'], fmt="%d")

maf_LB = 0.01

def recode_genomes(chromo):
    allele_range = var_info['#CHROM'] == chromo
    print(f"Saving haplotype {chromo}")
    print(f"chromo {chromo}, party {party_id}")
    var_info_chr = var_info[allele_range]
    nvariants=len(var_info_chr)

    def alleles2id(m):
        # Map genotype probs into ACGT encoding
        return np.unique(np.stack([var_info_chr['REF'], var_info_chr['ALT']]), return_inverse=True)[1] \
            .astype(np.int16).reshape((m, 2))

    allele_ids = alleles2id(nvariants).transpose()
    fname_hap = f"{args_local['haps_dir']}/chr{chromo}"

    # select the corresponding rows in buf
    hap_chr = buf[:, allele_range]
    haps = hap_chr * allele_ids[0][np.newaxis,:] + (1 - hap_chr) * allele_ids[1][np.newaxis,:]
    haps.shape
    np.save(fname_hap, haps)

for chromo in range(1, 23):
    recode_genomes(chromo)

def to_geno_count(x):
    return x[::2] + x[1::2]
geno_matrix = to_geno_count(buf)


king_snps = np.loadtxt(args_local['snp_list'], dtype=str)
king_filt = var_info['ID'].isin(king_snps)

geno = geno_matrix[:, king_filt].astype(np.int8)
geno.tofile(f"{args_local['geno_dir']}/all_chrs.bin")

nn = len(psam)
mm = geno.shape[1]
mmm = len(var_info)
print(f"Number of samples: {nn}, number of SNPs: {mm}, number of variants: {mmm}")
print(f"Number of rows in pgen: {n}, number of SNPs in pgen: {m}")
