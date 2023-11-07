#%%
from my_util import *
from user_config import conf
import pandas as pd
import numpy as np

%cd notebooks

if __name__ == "__main__":
    conf['n_job'] = 16
    exec(param_from_bash())
    args = eval(ARG_STR)
    rid, rep = read_cmd_args(conf, args=args)
    test_dir = conf['test_dir']
    king_file = f"{test_dir}/ground_truths/KING.dat"
    rep = 1
    related = read_king_file_raw(king_file)
    rel_size = related.shape[0]
    print(f"num (+) samples = {related.shape[0]}")
    lists = [related['ID' + str(i)].to_numpy() for i in range(2)]
    init_conf(conf, n)
    conf['no_hashing'] = True
    print_summary(conf)
    M = 0
    for chromo in range(1, 23):
        filt = np.load(f"data/maf/chr{chromo}.npy")
        M += filt.sum()

    print(f"#SNP = {M}")
    snp_ratios = [0.1, 0.4, 0.7, 0.8, 0.9, 1]
    
    save_path = conf['test_dir'] + "sketched/"
    run_command(['mkdir', save_path])
    for ratio in snp_ratios:
        snp_len = int(ratio * M)
        print(snp_len)
        snp_range = np.random.choice(M, (snp_len), replace=False)
        snp_range = np.sort(snp_range).astype(np.int32)
        np.savez(save_path + f"SNPs_{ratio * 100:.0f}.npz", snp_range = snp_range)
        with open(save_path + f"num_SNP{ratio*100:.0f}.txt", "w") as f:
            f.write(f"{snp_len}")
# %%
