from kmers import HaplotypeArray, KmerBuilder
from my_util import *
from pathlib import Path
import pandas as pd
import multiprocessing as mp
import logging
import concurrent
import ctypes
from datetime import datetime
import pickle
from joblib import Parallel, delayed

info = mp.get_logger().info

import numpy as np

verbose = 30
buck_sat = 0
LL = 0
chr_cc = 0
permuted = False
mod_Lehmer = 4294967311 # 2^32 + 15 (from lLASH implementation)
coef_uHash_range = 536870912 + 1 # 2^29 + 1 (from l^2LSH manual)
save_kmer = False
ENT_HI = np.inf
ENT_LO = 0
tmp_table = True


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


def hash_batch(arg, half_mode=False):
    assert(not half_mode) # only implemented in cached_kmer mode
    start = start_timer(lambda x: x)
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
        if save_kmer and trial_id == 0:
            np.save(f"{conf['kmer_dir']}/{tid}.npy", kmers)
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


def hash_worker_directly_load_kmer(arg, half_mode=False):
    conf, tid, h_keys = arg
    LSH_function = minHash_LSH if conf['LSH_mode'] == 'minHash' else Hamming_LSH if conf['LSH_mode'] == 'Hamming' else exact_LSH
    trim_cap = conf['trim_cap'] 

    all_table = list()
    L = conf['L']

    if not half_mode:
        kmers = np.load(f"{conf['kmer_dir']}/{tid}.npy", mmap_mode='r')
    else:
        kmers = np.load(f"{conf['kmer_dir']}/{tid}.npy", mmap_mode='r')[:conf['n'],:]
    for trial_id in range(L):
        # for this split version, we are using the same subsampled indices every time
        uhash_key = h_keys[trial_id][1]
        mini_n2, _ = kmers.shape
        bits = np.column_stack([LSH_function(h_key, kmers) for h_key in h_keys[trial_id][0]])
        ell = bits.shape[1]
        hash_values = \
            apply_uhash(bits, conf['hash_mode'], ell, ell, trial_id, uhash_key) % conf['key_space_size']
        
        # this gives a table with two rows, first being 0 0 1 1 2 2
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


def find_params(conf, party):
    conf = conf.copy()
    if party is not None:
        conf['test_dir'] += f'party{party}/'
        conf['kmer_dir'] += f'party{party}/'

    conf['num_stdvar'] = 10

    return conf

def find_overall_trial_ID(conf, chr_id, seg_id):
    max_segs = conf['max_segs']
    s = sum(max_segs[conf['chr_begin']:chr_id])
    return s + seg_id

def sample_LSH_params(conf):
    print(f"sampling LSH params")
    chr_begin = conf['chr_begin']
    chr_end = conf['chr_end']
    L = conf['L']

    max_segs = conf['max_segs']
    max_num_trial = conf['max_num_trial']
    actual_num_trial = max_num_trial

    indices = np.full((max_num_trial, L), None, dtype=np.object_)
    h_keys = np.full((max_num_trial, L), None, dtype=np.object_)

    chr_range = conf['chr_range']
    cMs = HaplotypeArray.read_cM(conf, chr_range)

    builder = KmerBuilder(conf)
    for chr_id in range(chr_begin, chr_end):
        max_seg = max_segs[chr_id]
        for seg_id in range(max_seg):
            tid = find_overall_trial_ID(conf, chr_id, seg_id)
            try:
                conf_here = conf.copy()
                conf_here['test_dir'] += f'/party{0}/'
                hap_array = HaplotypeArray(conf_here, chr_id, seg_id, cMs[chr_id])
                mafs = hap_array.read_mafs()
                cM_begin = seg_id * conf['step_cM'] if conf['mode'] >= 2 else None
                for trial_id in range(L):
                    indices[tid, trial_id] = builder\
                        .sample_indices(mafs, hap_array.get_cM(), cM_begin).astype(np.int32)
                    kmer_enc_len = indices[tid, trial_id].shape[0] // conf['inc']
                    if kmer_enc_len == 0:
                        raise ValueError(f"not enough SNPs in this segment")
                    # sample hash keys according to LSH_mode in conf
                    if conf['LSH_mode'] == 'Hamming':
                        h_keys[tid, trial_id] = ([(np.random.randint(0, kmer_enc_len, dtype=np.int32),) for _ in range(conf['l'])],
                                                 [np.random.randint(1, coef_uHash_range, dtype=np.uint32)
                                                 for _ in range(max(conf['l'], conf['k']))] if conf['hash_mode'] == 2 else None)
                                                  
                    elif conf['LSH_mode'] == 'minHash':
                        # every key is a tuple of (a, b) where a, b are random integers in [0, p-1]
                        h_keys[tid, trial_id] = ([(np.random.randint(mod_Lehmer, dtype=np.uint64), 
                                                  np.random.randint(mod_Lehmer, dtype=np.uint64)) for _ in range(conf['l'])],
                                                 [np.random.randint(1, coef_uHash_range, dtype=np.uint32)
                                                 for _ in range(max(conf['l'], conf['k']))] if conf['hash_mode'] == 2 else None)
                    elif conf['LSH_mode'] == 'exact':
                        h_keys[tid, trial_id] = ([None for _ in range(conf['l'])], 
                                                 [np.random.randint(1, coef_uHash_range, dtype=np.uint32)
                                                 for _ in range(max(conf['l'], conf['k']))] if conf['hash_mode'] == 2 else None)
                    else:
                        raise(NameError(f"LSH_mode {conf['LSH_mode']} not supported"))
            except ValueError as e:
                print(f"failed to sample LSH params for chr {chr_id} seg {seg_id} with tid={tid}", file=sys.stderr)
                # print(f"might be because not enough SNPs in this segment", file=sys.stderr)
                print(e, file=sys.stderr)
                indices[tid, 0] = None
                actual_num_trial -= 1
                continue
            except ZeroDivisionError as e:
                print(f"failed to sample LSH params for chr {chr_id} seg {seg_id} with tid={tid}", file=sys.stderr)

    conf['actual_num_trial'] = actual_num_trial 
    out_dir = conf['out_dir']
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    if save_kmer:
        for i in range(2):
            Path(out_dir + f"kmer/party{i}/").mkdir(parents=True, exist_ok=True)

    # optionally random order for merging the segments
    # in paper, the default preference is to choose from lowest to highest indices
    merge_order = np.full((L), None, dtype=np.object_)
    for i in range(L):
        if i % 2 == 0:
            merge_order[i] = np.random.permutation(max_num_trial)
        else:
            merge_order[i] = merge_order[i-1][::-1]  


    np.savez(out_dir + f"LSH_params", indices=indices, h_keys=h_keys, merge_order=merge_order)
    # pickle file for later use
    # change ndarray to list for dict()
    conf['max_cMs'] = conf['max_cMs'].tolist()
    pickle.dump(conf, open(out_dir + f"LSH_params_conf.pkl", "wb"))
    with open(out_dir + f"LSH_params_conf.txt", "w") as f:
        f.write(print_summary(conf, out=False))
    

def ad_hoc_order(max_num_trial, flip):
    first_half = np.arange(0, max_num_trial, 2)
    merge_order_now = first_half
    if flip == 1:
        merge_order_now = merge_order_now[::-1]
    return merge_order_now


def build_LSH_table(conf_old, indices, h_keys, party=0, merge_order=None, hash_batch_func=hash_batch, hash_with_encode=True, tmp="", half_mode=False):
    global buck_sat, LL, tid_total_count
    buck_sat = 0
    tid_total_count = 0
    conf = find_params(conf_old, party)

    n_job = conf['n_job']
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
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_job, initializer=thread_init,
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
                if hash_with_encode:
                    params[tid] = [conf, chr_id, seg_id, beg, end, indices[tid, :], h_keys[tid]]
                else:
                    params[tid] = [conf, tid, h_keys[tid]]

        if merge_order is None:
            merge_order_now = ad_hoc_order(max_num_trial, (L - LL) % 2)
        else:
            merge_order_now = merge_order[LL - 1]
        for tid in merge_order_now:
            if params[tid] is not None:
                executor\
                    .submit(hash_batch_func, params[tid], half_mode)\
                    .add_done_callback(collect_table)

    stop_timer(start)

    B = cap_merged
    print(f"========== Table {party } Summary ==========")
    print(f"max B = {B}")
    print(f"buck_sat = {100 * buck_sat / (N * cap_merged):.2f}%")
    save_table(conf_old, party, table, tmp)
    return None


def save_table(conf, party, hash_table, tmp=""):
    out_dir = conf['out_dir'] + f"tables{tmp}/party{party + 1}/"
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    N, B = hash_table.shape
    id_table = hash_table.flatten()
    np.savez(out_dir + f"ID_table", ID_table=id_table, N=N, B=B)


def find_per_sample_recall(conf, t, B, half_mode=False):
    n = conf['n']
    if half_mode:
        n //= 2
    num_degree = len(KING_THRESHOLDS) - 1
    print(num_degree)
    def find_max_deg(related):
        rst = [num_degree * np.ones((n), int) for _ in range(2)]
        for (deg, data) in related:
            data = data[(data['ID0'] < n) & (data['ID1'] < n)]
            for j in range(2):
                ids = data[f'ID{j}'].values
                rst[j][ids] = np.minimum(rst[j][ids], deg)
        return rst
    gt_lists = find_max_deg(read_king_file(f"{conf['test_dir']}{conf['ground_truth']}"))
    t0 = t[0]
    t1 = t[1]
    related = read_king_file_raw(f"{conf['test_dir']}{conf['ground_truth']}")
    related_pairs = (n * related.ID0.astype(np.int64) + related.ID1).values
    lists = [np.zeros((n), dtype=bool) for _ in range(2)]
    for i in range(B):
        for j in range(B):
            # use n * a + b to represent (a, b)
            a = t0[:, i]
            b = t1[:, j]
            # must deal with -1
            gg = (a == -1) | (b == -1)
            a = a[~gg]
            b = b[~gg]
            pairs = a.astype(np.int64) * n + b
            found = np.isin(pairs, related_pairs)
            a = a[found]
            b = b[found]
            lists[0][a] = True
            lists[1][b] = True

    rst = np.zeros((num_degree), dtype=np.float32)
    total_cnt = 0
    total_n = 0
    for j in range(2):
        now_list = gt_lists[j]
        for deg in range(num_degree):
            gt_locs = np.nonzero(now_list == deg)[0]
            n_deg = len(gt_locs)
            total_n += n_deg / 2
            cnt = lists[j][gt_locs].sum()
            total_cnt += cnt / 2
            found_rate = cnt/n_deg
            if j == 1:
                print(bcolors.OKBLUE, end='')
            else:
                print(bcolors.OKGREEN, end='')
            if n_deg == 0:
                print(f"** No deg {deg} is present in dataset")
            else:
                print(f"** deg {deg} found rate = {found_rate * 100:.4f}% = {cnt}/{n_deg} **")
            rst[deg] += (found_rate) / 2
    print(f"!!    ** total = {total_cnt / total_n * 100:.4f}% = {total_cnt}/{total_n} ** {bcolors.ENDC}")
    return rst


def load_LSH_params(conf, no_merge_order=False):
    out_dir = conf['out_dir']
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    params = np.load(out_dir + f"LSH_params.npz", allow_pickle=True)
    if not no_merge_order:
        return params['indices'], params['h_keys'], params['merge_order']
    else:
        return params['indices'], params['h_keys']


def sample_cols(conf, hash_cols=True):
    indices, h_keys = load_LSH_params(conf, True)
    if hash_cols is False:
        for i, indice in enumerate(indices): 
            h_keys[i] = np.arange(indice.shape[0])
    return h_keys


def sample_hashes(conf, party=None, seg_range=None, h_key_idx=0, cols=None):
    # cols = Parallel(n_jobs=conf['n_job'])(delayed(sample_cols)(conf, 0, i) for i in range(conf['max_num_trial']))
    if seg_range is None:
        seg_range = range(conf['max_num_trial'])
    if cols is None:
        cols = sample_cols(conf)

        seg_used = np.array([cols[i][0] is not None for i in seg_range])
    else:
        seg_used = cols is not None
    # read in kmers
    if party is None:
        kmers = [Parallel(n_jobs=conf['n_job'],verbose=10)(delayed(compute_hash_values)(conf, party, i, cols[i], seg_used[i], h_key_idx) for i in seg_range) for party in range(2)]
    else:
        if type(seg_used) == bool:
            kmers = compute_hash_values(conf, party, seg_range[0], cols, seg_used, 0)
        else:
            kmers = Parallel(n_jobs=conf['n_job'], verbose=10)(delayed(compute_hash_values)(conf, party, i, cols[i], seg_used[i], h_key_idx) for i in seg_range)
    return seg_used,kmers


def compute_hash_values(conf, party, i, h_keys, used, h_key_idx=0):
    # info(i)
    LSH_function = minHash_LSH if conf['LSH_mode'] == 'minHash' else Hamming_LSH if conf['LSH_mode'] == 'Hamming' else exact_LSH
    n = conf['n']
    if used:
        kmers = np.load(f"{conf['kmer_dir']}party{party}/{i}.npy", mmap_mode='r')
        assert (conf['hash_mode'] == 1)
        # // there is a problem here
        bits = np.column_stack([LSH_function(h_key, kmers) for h_key in h_keys[0][h_key_idx]])
        ell = bits.shape[1]
        return apply_uhash(bits, conf['hash_mode'], ell, ell, i, h_keys[1]) % conf['key_space_size']
    else:
        return np.zeros((2 * n), dtype=np.uint32)


def load_conf(conf, args, cap=1, tau=64, half_mode=False):
    out_dir = conf['out_dir']
    # read conf from pickled file
    conff = pickle.load(open(out_dir + f"LSH_params_conf.pkl", "rb"))
    # re-enter some args from command line if the old one is outdated.
    read_cmd_args(conff, args)
    conff['cap_each'] = conff['cap_merged'] = cap
    if cap > 2:
        conff['L'] = 1
    elif cap > 1:
        conff['L'] = 2
    
    conff['n_job'] = 32
    conff['key_space_ratio'] = tau
    conff['key_space_size'] = int(np.ceil(conf['n'] * 2 * tau))
    if half_mode:
        conff[f'key_space_size'] = int(np.ceil(conf['n'] * conff['key_space_ratio']))

    return conff

def build_tables(conf, parties,hash_batch_func=hash_batch, hash_with_encode=True, args=None, no_merge_order=False, cap=1, tmp="", tau=64, half_mode=False):
    indices, h_keys, merge_order = load_LSH_params(conf)
    conf_to_use = load_conf(conf, args, cap, tau, half_mode)
    print(print_summary(conf_to_use, out=False), file=sys.stderr)

    untrimmed_cnts = []
    for party in parties:
        if not no_merge_order:
            untrimmed_cnts.append(build_LSH_table(conf_to_use, indices, h_keys, party, merge_order, hash_batch_func, hash_with_encode, tmp, half_mode))
        else:
            untrimmed_cnts.append(build_LSH_table(conf_to_use, indices, h_keys, party, None, hash_batch_func, hash_with_encode, tmp, half_mode))
    return untrimmed_cnts

def run_two_party_exp(conf, n, parties=None, found=None, no_param=False, no_check=False, args=None, refresh_param=False, tmp="", cap=1, tau=64, half_mode=False):
    conf = conf.copy()
    start = start_timer()
    parties = range(conf['n_party'])
    if not no_param:
        sample_LSH_params(conf)
    if no_param:
        untrimmed_cnts = build_tables(conf, parties, hash_worker_directly_load_kmer, False, args=args, no_merge_order=True, cap=cap, tmp=tmp, tau=tau, half_mode=half_mode)
    else:
        untrimmed_cnts = build_tables(conf, parties, args=args, cap=cap, tmp=tmp, tau=tau, half_mode=half_mode)
    print("========== Summary ==========")
    conff = find_params(conf, 0)
    print_summary(conff)
    if no_check:
        stop_timer(start)
        return 0, untrimmed_cnts
    # report false negative rates here
    recalls = check_table(conf, tmp, half_mode=half_mode)
    return recalls, untrimmed_cnts

def check_table(conf, tmp="_tmp", cap=1, half_mode=False):
    start = start_timer()
    bb, N, B = load_tables(conf, tmp)
    for i in range(2):
        bb[i] = bb[i][:int(len(bb[i]) * cap)]
    conff = find_params(conf, 0)
    recalls = find_per_sample_recall(conf, bb, B, half_mode=half_mode)
    print(f"max B = {B}")
    n = conff['n']
    if half_mode:
        n //= 2
    key_space_size = bb[0].shape[0]
    print(f"# cmp w/padding = {B * B * key_space_size}/{n * n} = {B * B * key_space_size / (n * n) * 100:.3f}%")
    print("----------------------------------------")
    stop_timer(start)
    return recalls

if __name__ == '__main__':
    # the actual test
    from user_config import conf
    exec(param_from_bash())
    args = eval(ARG_STR)
    rid, rep = read_cmd_args(conf, args=args)
    assert conf['hash_mode'] <= 1
    n = conf['n']
    init_conf(conf, n, hash_id=conf['hid'], create_out=True)
    # if half_mode, use the first 50K in each file to get the UKB-50K dataset
    half_mode=False
    if half_mode:
        conf['key_space_size'] = int(np.ceil(conf['key_space_size'] / 2))
    recalls, table_cnts = run_two_party_exp(conf, n, no_param=False, args=args, tmp='_tmp', half_mode=half_mode, tau=64)
