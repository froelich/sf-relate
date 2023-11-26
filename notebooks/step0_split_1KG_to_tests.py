import numpy as np
import pgenlib
from my_util import *
import pandas as pd
from user_config import *
from joblib import Parallel, delayed

def read_pgen(data_path):
    reader = pgenlib.PgenReader(data_path.encode("utf-8"))
    n = 2 * reader.get_raw_sample_ct()
    m = reader.get_variant_ct()
    phase_present = np.zeros((m, n), dtype=np.int8)
    buf = np.zeros((m, n), dtype=np.int32)
    reader.read_alleles_and_phasepresent_list(np.arange(0, m, dtype=np.uint32), buf, phase_present)
    return n, m, buf.transpose().astype(np.int8)

n, m, buf = read_pgen("data/1KG/thinned.pgen")
n2, m2, buf_king = read_pgen("data/1KG/snpsForKING.pgen")
assert(n == n2)

def to_geno_count(x):
    return x[::2] + x[1::2]
geno_count_king = to_geno_count(buf_king)
geno_count = to_geno_count(buf)

psam = pd.read_csv("data/1KG/snpsForKING.psam", sep='\t')
king_table = pd.read_csv("data/1KG/king.kin0", sep='\t')
print(n, m, m2)

# %%
# all_samples = pd.read_csv(f"data/1KG/", sep=' ', skiprows=2, header=None)[0]
all_samples = psam
bbank_n = len(all_samples)
n = bbank_n // 2 
k = 2
test_dir = f"data/{k}party_{n}/"
list_names = [test_dir + f'party{i}/list' for i in range(k)]

def get_var_info(name='thinned'):
    try:
        var_info = pd.read_csv(f"data/1KG/{name}.var", sep='\t')
    except FileNotFoundError:
        # remove all lines with ## in the beginning of the file
        run_command(['sed', '-i', '/^##/d', f"data/1KG/{name}.pvar"])
        var_info = pd.read_csv(f"data/1KG/{name}.pvar", sep='\t', usecols=['#CHROM', 'ID', 'POS', 'REF', 'ALT'])
        var_info.to_csv(f"data/1KG/{name}.var", sep='\t', index=False)
        run_command(['rm', '-r', f"data/1KG/{name}.pvar"])
    
    return var_info

maf = np.sum(geno_count, axis=0) / (2 * n)
var_info=get_var_info()
var_info_king=get_var_info("snpsForKING")

run_command(['mkdir', f'{test_dir}/maf/0.01', '-p'])
run_command(['mkdir', f'{test_dir}/pos/0.01', '-p'])
# split by chromosome
for i in range(1, 23):
    snp_range = var_info['#CHROM'] == i
    var_info_chr = var_info[snp_range]
    var_info_chr.to_csv(f"{test_dir}/maf/0.01/chr{i}.var", sep='\t', index=False)
    maf_chr = maf[snp_range]
    np.savetxt(f"{test_dir}/maf/0.01/chr{i}.txt", maf_chr)
    np.savetxt(f"{test_dir}/pos/0.01/chr{i}.txt", var_info_chr.loc[snp_range, 'POS'])
    write_number(f"{test_dir}/maf/0.01/chr{i}_numsnps.txt", len(var_info_chr))

# %%
# this part confirms the kinship coefficients computed make sense
# king_computed = []
# king_subset = []
# for tup in king_table.itertuples():
#     a = tup[1]
#     b = tup[2]
#     i, j= np.where(psam['#IID'] == a)[0][0],np.where(psam['#IID'] == b)[0][0]
#     king_computed.append(king_formula(geno_count[i, :], geno_count[j, :]))
#     king_subset.append(tup[-1])
# df = pd.DataFrame({'king_raw': king_computed, 'king_UKB_processed': king_subset})
# df.to_csv("data/1KG/king_compare.csv", index=False)
# print(df)
#%%

# only compute KING when we want to split by placing relatives on two sides first
# recompute_KING(test_dir, k)
list_ppl = [[], []]
def split_pgen(test_dir, listnames, all_samples, n, k=2):
    run_command(['mkdir', f'{test_dir}'])
    print("Generating test sets")
    used = np.zeros(len(all_samples), dtype=bool)
    for tup in king_table.itertuples():
        a = tup[1]
        b = tup[2]
        i, j= np.where(all_samples == a)[0][0],np.where(all_samples == b)[0][0]
        b = np.random.randint(2)
        if b == 0:
            i, j = j, i
        if used[i] or used[j]:
            continue
        used[i] = True
        used[j] = True
        list_ppl[0].append(i)
        list_ppl[1].append(j)

    unused = np.nonzero(~used)[0]
    half = len(unused) // 2
    print(len(unused))
    for i in range(k):
        list_ppl[i] += unused[half * i: half * (i + 1)].tolist()

    print(len(list_ppl[i]))

    for i in range(k):
        np.random.shuffle(list_ppl[i])

    for i in range(len(listnames)):
        run_command(['mkdir', f'{test_dir}/party{i}'])
        with open(listnames[i] + '.txt', "w") as file:
            for ID in list_ppl[i]:
                # ID = int(ID)
                file.write(str(ID) + ' ' + str(ID) + '\n')

    print(f"Sampled a total of {len(np.unique(np.concatenate((list_ppl[0], list_ppl[1]))))} individuals")

    # recode the haplotypes
    print("Recode the haplotypes")

split_pgen(test_dir, list_names, all_samples['#IID'].values, n, k)

from joblib import Parallel, delayed
maf_LB = 0.01
haps_dir = 'haps/'
geno_dir = 'geno/'
for i in range(k):
    party_dir = f"{test_dir}party{i}/"
    run_command(['mkdir', f'{party_dir}{haps_dir}'])
    run_command(['mkdir', f'{party_dir}{geno_dir}'])
def recode_genomes(chromo, i):
    allele_range = var_info['#CHROM'] == chromo
    print(f"Recoding {chromo, i}")
    print(f"chromo {chromo}, party {i}")
    party_dir = f"{test_dir}party{i}/"
    var_info_chr = var_info[allele_range]
    allele_range_king = var_info_king['#CHROM'] == chromo
    nvariants=len(var_info_chr)

    def alleles2id(alleles, m):
        # Map genotype probs into ACGT encoding
        return np.unique(np.stack([var_info_chr['REF'], var_info_chr['ALT']]), return_inverse=True)[1] \
            .astype(np.int16).reshape((m, 2))

    allele_ids = alleles2id(var_info, nvariants).transpose()
    print(allele_ids.min(), allele_ids.max())

    fname_hap = f"{party_dir}{haps_dir}chr{chromo}_maf{maf_LB:.2f}"
    fname_gen = f"{party_dir}{geno_dir}chr{chromo}_maf{maf_LB:.2f}"

    # parse haplotypes data
    def find_row_index(idx):
        # first duplicate idx, but keep the order
        idx = 2 * np.vstack([idx, idx]).reshape(-1, order='F')
        idx[1::2] += 1
        return idx


    # select the corresponding rows in buf
    indices = find_row_index(list_ppl[i])
    hap_chr = buf[indices][:, allele_range]
    haps = hap_chr * allele_ids[0][np.newaxis,:] + (1 - hap_chr) * allele_ids[1][np.newaxis,:]
    haps.shape
    np.save(fname_hap, haps)

    genotype_scores = geno_count_king[list_ppl[i]][:, allele_range_king].astype(np.int8)
    np.save(fname_gen, genotype_scores)

Parallel(n_jobs=8)(delayed(recode_genomes)(chromo, i) for chromo in (range(22, 0, -1)) for i in range(k))

chr_range = range(22, 0, -1)
def piece_all_chromo(i):
    all_geno = [] 
    fname_all_geno = f"{test_dir}party{i}/{geno_dir}all_chrs_maf{maf_LB:.2f}"
    for chromo in chr_range:
        print(f"reading {chromo}", end='...')
        all_geno.append(np.load(test_dir + f"party{i}/{geno_dir}chr{chromo}_maf{maf_LB:.2f}.npy"))
        print("done")
    whole = np.hstack(all_geno)
    print("HI")
    np.save(fname_all_geno, whole)
    run_command(['mkdir', f"{test_dir}{geno_dir}/party{i + 1}", '-p'])
    # convert the npy files to binary files for MHE
    whole.tofile(f"{test_dir}{geno_dir}/party{i + 1}/all_chrs_maf{maf_LB:.2f}.bin")
    print("BYE")

for i in range(k):
    piece_all_chromo(i)

exec(param_from_bash())
args = eval(ARG_STR)
rid, rep = read_cmd_args(conf, args=args)
k = 2
ver = 0
start = start_timer()
out_dir, ratio = f'data/2party_{n}/', 0.01
chr_range = range(22, 0, -1)
ratio = 0.01
# batch_splits = [0] + [4] * 17 + [2] * 3 + [1] * 2
batch_cnt = 8
read_suff = ''
out_suff = ''
permuted = False 
batch_size = int(np.ceil(n / batch_cnt))
        
# # split all chromo into batches so that we can compute matrix product in memory
party_range = [0, 1]
for i in party_range:
    fname_all_geno = f"{out_dir}party{i}/geno/all_chrs_maf{ratio:.2f}{read_suff}.npy"
    mat = np.load(fname_all_geno)
    for j in range(batch_cnt):
        b_begin = j * batch_size
        b_end = (j + 1) * batch_size
        np.save(out_dir + f"party{i}/geno/all_chrs_maf{ratio:.2f}_{j}{read_suff}_p{i}", mat[b_begin:b_end])

def compute_king(p1, local=None):
    party0 = 0
    party1 = 1
    if local is not None:
        party0 = local
        party1 = local
    p1file = f"{out_dir}party{party0}/geno/all_chrs_maf{ratio:.2f}_{p1}{read_suff}_p{party0}.npy"
    outfile = f"{out_dir}king_{p1}{out_suff}"

    x = np.load(p1file)

    x2 = np.nansum(x.astype(np.float32) ** 2, axis = 1).ravel()
    print("done reading gX")
    kings = np.zeros((x.shape[0], n), dtype=np.float32)

    for p2 in range(batch_cnt):
        p2file = f"{out_dir}party{party1}/geno/all_chrs_maf{ratio:.2f}_{p2}{read_suff}_p{party1}.npy"
        y = np.load(p2file)
        b_begin = p2 * batch_size
        b_end = (p2+1) * batch_size
        print(f"working on p1, p2 = {p1, p2} done reading gY")

        y2 = np.nansum(y.astype(np.float32) ** 2, axis = 1).ravel()
        print("A", end='..')

        numer = x2[:,np.newaxis] -2 * x.astype(np.float32) @ y.astype(np.float32).T
        numer += y2[np.newaxis,:]
        print("B", end='..')

        hetx = np.nansum(x == 1, axis = 1, dtype=np.int32).ravel()
        hety = np.nansum(y == 1, axis = 1, dtype=np.int32).ravel()

        print("C", end='..')

        kings[:,b_begin:b_end] = 0.5 - 0.25 * numer /np.minimum(hetx[:,np.newaxis], hety[np.newaxis,:]) #denom
        print("done computing")
    np.save(outfile, kings)
    
Parallel(n_jobs=8)(delayed(compute_king)(p1, None) for p1 in range(0, batch_cnt))
stop_timer(start)

suff = out_suff
start = start_timer()
import numpy as np
def decode_submatrix(p1):
    return np.load(f"{out_dir}king_{p1}{suff}.npy")

all_king = np.vstack([decode_submatrix(p1) for p1 in range(batch_cnt)])
run_command(['mkdir', f'{out_dir}ground_truths'])
np.save(f"{out_dir}ground_truths/all_king_maf{ratio:.2f}{suff}", all_king)
stop_timer(start)

start = start_timer()
out_dir, ratio = f'data/2party_{n}/', 0.01
ratio = 0.01
permuted = False 
output_kinship = KING_THRESHOLDS[-1]
print(output_kinship)

records = []
sample_from_thr(all_king, records, output_kinship, target=None, verbose=True)
fname = f"{out_dir}ground_truths/KING{suff}{'_all' if output_kinship == 0 else ''}.dat"
write_king_file(fname, records)

# %%
start = start_timer()
out_dir, ratio = f'data/2party_{n}/' + (f'v{ver}' if ver != 0 else ''), 0.01
ratio = 0.01
permuted = False 
output_kinship = KING_THRESHOLDS[-1]
print(output_kinship)

records = []
sample_from_thr(all_king, records, output_kinship, target=None, verbose=True)
fname = f"{out_dir}ground_truths/KING{suff}{'_all' if output_kinship == 0 else ''}.dat"
write_king_file(fname, records)

# remove the intermediate files
for i in party_range:
    fname_all_geno = f"{out_dir}party{i}/geno/all_chrs_maf{ratio:.2f}{read_suff}.npy"
    for j in range(batch_cnt):
        run_command(["rm", out_dir + f"party{i}/geno/all_chrs_maf{ratio:.2f}_{j}{read_suff}_p{i}.npy"])

# remove tmp files
for p1 in range(0, batch_cnt):
    run_command(['rm', f"{out_dir}king_{p1}{out_suff}.npy"])
