from kmers import HaplotypeArray, KmerBuilder
from my_util import *
from pathlib import Path
import multiprocessing as mp
import logging
import concurrent
import concurrent.futures
from default_config import conf
import ctypes
import tomlkit
import os
import pickle

info = mp.get_logger().info

import numpy as np


def Hamming_LSH(k, x):
    assert (k[0] < len(x))
    if len(x.shape) == 1:
        return x[k[0]]
    else:
        return x[:, k[0]]

def exact_LSH(k, x):
    return x

def minHash_LSH(k, x):
    a = k[0]
    b = k[1]
    # compute the minhash value of each row in x
    if len(x.shape) == 1:
        return np.min((a * x + b) % mod_Lehmer)
    else:
        return np.min((a * x + b) % mod_Lehmer, axis=1)


def hash_batch(arg):
    conf, chr_id, seg_id, beg, end, indices, h_keys = arg

    hap_array = HaplotypeArray(conf, chr_id, seg_id, permuted=permuted, beg=beg, end=end)

    # window + k-mers
    tid = find_overall_trial_ID(conf, chr_id, seg_id)

    LSH_function = minHash_LSH if conf['LSH_mode'] == 'minHash' else Hamming_LSH if conf['LSH_mode'] == 'Hamming' else exact_LSH
    trim_cap = conf['trim_cap'] 

    builder = KmerBuilder(conf)
    all_table = list()
    L = conf['L']

    for trial_id in range(L):
            
        # Can make them share the same kmers if we don't project every time
        rand_seed = tid * L + trial_id
        uhash_key = h_keys[trial_id][1]
        kmers = builder.get_kmers(hap_array, rand_seed, indices[trial_id], h_key=uhash_key)
        mini_n2, _ = kmers.shape

        bits = np.column_stack([LSH_function(h_key, kmers) for h_key in h_keys[trial_id][0]])
        ell = bits.shape[1]
        hash_values = \
            apply_uhash(bits, conf['hash_mode'], ell, ell, trial_id, uhash_key) % conf['key_space_size']
        
        table = np.vstack([np.int32(np.arange(mini_n2) / 2), hash_values.T])

        # only keep the first unique items
        table = table.transpose()
        hashes_cnts = np.unique(table[:, 1], return_counts=True)[1]
        inverse = np.unique(table[:, 1], return_inverse=True)[1]
        table = table[hashes_cnts[inverse] <= trim_cap] 
        np.random.shuffle(table)
        # got the first few unique hashes
        to_adds = []
        for _ in range(conf['cap_each']):
            unique_hash_locs = np.unique(table[:, 1], return_index=True)[1]
            trim_cap = conf['trim_cap'] 
            to_adds.append(table[unique_hash_locs])
            table= table[~np.isin(np.arange(0, len(table)), unique_hash_locs)]
        to_add_items = np.vstack(to_adds)

        all_table.append(to_add_items)

    return tid, all_table


def thread_init():
    pass


def ad_hoc_order(max_num_trial, flip):
    first_half = np.arange(0, max_num_trial, 2)
    merge_order_now = first_half
    if flip == 1:
        merge_order_now = merge_order_now[::-1]
    return merge_order_now


def build_LSH_table(conf_old, indices, h_keys):
    global buck_sat, LL, tid_total_count
    buck_sat = 0
    tid_total_count = 0
    conf = conf_old.copy()

    LL = conf['L']
    cap_merged = conf['cap_merged']
    L = conf['L']
    max_num_trial = conf['max_num_trial']
    actual_num_trial = conf['actual_num_trial']
    if merge_order is None:
        merge_order_now = ad_hoc_order(max_num_trial, 0)
        # need to subtract back the empty segments
        actual_num_trial = len(merge_order_now) - sum(trial is None for trial in indices[:, 0])
    chr_range = conf['chr_range']

    '''
    Compute expected # items in bucket
    '''
    start = start_timer()

    logger = mp.log_to_stderr()
    logger.setLevel(logging.INFO)
    N = conf['key_space_size']
    table = np.full((N, cap_merged), -1, dtype=np.int32)
    table_cnt = np.zeros((N), dtype=ctypes.c_int8)
    not_full = np.ones((N), dtype=np.bool_)

    print(f"actual # trials with repeat = {actual_num_trial * L}")
    # enough space to hold everything
    partial_tables = [[None for _ in range(L)] for _ in range(max_num_trial)]

    def collect_table(rst):
        global buck_sat, LL, tid_total_count
        rst = rst.result()
        tid = rst[0]
        for trial_id in range(conf['L']):
            if buck_sat / (N * cap_merged) > 0.99995:
                return
            partial_tables[tid][trial_id] = rst[1][trial_id]

        def collect_table_based_on_tid(now):
            global LL, buck_sat, tid_total_count
            if buck_sat / (N * cap_merged) > 0.99995:
                LL = 0
                return
            if 0 < LL == len(partial_tables[now]):
                small_table = partial_tables[now][-1]
                for _ in range(cap_merged):
                    bucket_index = small_table[:, 1]
                    rows_to_add = np.unique(bucket_index, return_index=True)[1]
                    bucket_to_add = small_table[rows_to_add, 1]
                    # mostly using this branch
                    # add items one by one
                    to_add_not_full = not_full[bucket_to_add]
                    to_add_items = small_table[rows_to_add][to_add_not_full]
                    bucket_index = to_add_items[:, 1]
                    cnts = table_cnt[bucket_index]
                    table[bucket_index,cnts] = to_add_items[:, 0]
                    table_cnt[bucket_index] += 1
                    not_full[bucket_index] = table_cnt[bucket_index] < cap_merged
                    small_table = small_table[~np.isin(np.arange(len(small_table)), rows_to_add)]

                buck_sat = np.sum(table_cnt)

                if now % verbose == 0:
                    info(f"collecting table {now, LL - 1} with {buck_sat/(N * cap_merged) * 100:.2f}% buckets saturated")

                partial_tables[now].pop()

                tid_total_count += 1


        # all batches are done
        collect_table_based_on_tid(tid)

        # decrease global counter on repetitions too if all exp are done
        # and then need to collect until done
        while LL > 1 and tid_total_count == actual_num_trial:
            LL -= 1
            tid_total_count = 0
            # change the order in which to hash the chromosome to shared random order?
            if merge_order is None:
                merge_order_now = ad_hoc_order(max_num_trial, (L - LL) % 2)
            else:
                merge_order_now = merge_order[LL]
            for now in merge_order_now:
                collect_table_based_on_tid(now)

    cMs = HaplotypeArray.read_cM(conf_old, chr_range)
    params = [None for _ in range(max_num_trial)]
    with concurrent.futures.ProcessPoolExecutor(initializer=thread_init,
                                            initargs=()) as executor:
        for chr_id in chr_range:
            max_seg = conf['max_segs'][chr_id]
            for seg_id in range(max_seg):
                tid = find_overall_trial_ID(conf, chr_id, seg_id)
                if indices[tid, 0] is None:
                    # not enough data for this segment
                    info(f"skipping {chr_id, seg_id}")
                    continue
                beg, end = HaplotypeArray.get_boundaries(conf, seg_id, cMs[chr_id])
                params[tid] = [conf, chr_id, seg_id, beg, end, indices[tid, :], h_keys[tid]]

        if merge_order is None:
            merge_order_now = ad_hoc_order(max_num_trial, (L - LL) % 2)
        else:
            merge_order_now = merge_order[LL - 1]
        for tid in merge_order_now:
            if params[tid] is not None:
                executor\
                    .submit(hash_batch, params[tid])\
                    .add_done_callback(collect_table)

    stop_timer(start)

    B = cap_merged
    print(f"========== Table Summary ==========")
    print(f"max B = {B}")
    print(f"buck_sat = {100 * buck_sat / (N * cap_merged):.2f}%")
    save_table(conf_old, table)


def save_table(conf, hash_table):
    out_dir = conf['out_dir']
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    N, B = hash_table.shape
    id_table = hash_table.flatten()
    np.savez(out_dir + f"/ID_table", ID_table=id_table, N=N, B=B)

def load_LSH_params(param, no_merge_order=False):
    Path(param).mkdir(parents=True, exist_ok=True)
    params = np.load(param+ f"/LSH_params.npz", allow_pickle=True)
    if not no_merge_order:
        return params['indices'], params['h_keys'], params['merge_order']
    else:
        return params['indices'], params['h_keys']


def load_conf(param_dir):
    # read conf from pickled file
    conff = pickle.load(open(param_dir + f"/LSH_params_conf.pkl", "rb"))
    # re-enter some args from command line if the old one is outdated.
    return conff

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build Local Hash Table')
    parser.add_argument('-PARTY', type=int, help='party id', required=True)
    parser.add_argument('-FOLDER', type=str, help='path to the configuration folder', required=True)
    args = parser.parse_args()
    party_id = args.PARTY
    FOLDER = args.FOLDER

    args_local = tomlkit.load(open(f"{FOLDER}/configLocal.Party{party_id}.toml"))
    args_global = tomlkit.load(open(f"{FOLDER}/configGlobal.toml"))
    args_local['param'] = args_local['shared_param_dir']
    args_local['out'] = args_local['hash_table_dir']
    args_local['hap'] = args_local['haps_dir']
    args_local['pos'] = args_local['pos_dir']
    args_local['gmap'] = args_local['gmap_dir']
    args_local['pos'] = args_local['pos_dir']
    args_local['maf'] = args_local['maf_dir']
    namespace = argparse.Namespace()
    namespace.__dict__.update(args_global)
    namespace.__dict__.update(args_local)
    # find out which line in pvar_file start with #CHROM
    psam, snps, pvar = read_pgen_metadata(args_local)
    n, M, m = len(psam), len(snps), len(pvar)
    namespace.__dict__.update({'n': n, 'm': m})
    # save n, M, m to the hap directory
    # in txt
    with open(f"{args_local['hap']}/data_dim.txt", "w") as f:
        f.write(f"{n}\n")
        f.write(f"{M}\n")
        f.write(f"{m}\n")

    param = namespace.param
    conf = load_conf(param)
    indices, h_keys, merge_order = load_LSH_params(param)

    init_conf(conf, namespace)
    print_summary(conf)
    print("=" * 50)
    print("*" * 2 + " " * 3 +  "SF-Relate Step 1: Build Hash Tables" + " " * 3 + "*" * 2)
    build_LSH_table(conf, indices, h_keys)
    print(f"Output saved in {conf['out_dir']}" + ("/" if conf['out_dir'][-1] != "/" else ""))


#%% code for sanity checks
# import numpy as np
# # %ls notebooks/
# a = np.load("notebooks/trial/party0/table/ID_table.npz")['ID_table']
# b = np.load("notebooks/trial/party1/table/ID_table.npz")['ID_table']
# from default_config import conf
# conf['test_dir'] = 'notebooks/data/2party_1601/'
# a,b

# import pandas as pd
# n = 1601
# num_degree = 5
# print(num_degree)
# t = [a,b]
# related = pd.read_csv("notebooks/data/2party_1601/ground_truths/KING.dat", sep="\t")
# related
# related['pair'] = related['ID0'].astype(int) + related['ID1'].astype(int) * 10000
# related.set_index('pair', inplace=True)
# filt = (a != -1) & (b != -1)
# pairs = a[filt] + 10000 * b[filt]
# related.loc[pairs[np.isin(pairs, related.index)], 'found'] = True
# related['found'].fillna(False, inplace=True)
# related.groupby('deg').agg({'found': 'mean'})