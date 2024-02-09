import numpy as np
import time
import sys
import pandas as pd
from pathlib import Path
from datetime import datetime
import argparse
import pytz

verbose = 30
buck_sat = 0
LL = 0
chr_cc = 0
permuted = False
mod_Lehmer = 4294967311 # 2^32 + 15 (from lLASH implementation)
coef_uHash_range = 536870912 + 1 # 2^29 + 1 (from l^2LSH manual)

def find_overall_trial_ID(conf, chr_id, seg_id):
    max_segs = conf['max_segs']
    s = sum(max_segs[conf['chr_begin']:chr_id])
    return s + seg_id

def apply_uhash(data, hash_mode, k, inc, rand_seed, h_key = None):
    FNV_offset_basis = 0x811c9dc5
    FNV_prime = 0x01000193
    if len(data.shape) == 2:
        n2, m = data.shape
        enc_len = len(range(0, m, inc))
        if hash_mode != 1:
            raise NotImplementedError("Other hash functions are removed in the implmentation for simplicity")
        # use FNV-1a
        FNV_first_hash = ((FNV_offset_basis ^ rand_seed) * FNV_prime) & 0xFFFFFFFF

        hash = np.full((n2, enc_len), FNV_first_hash, dtype=np.int64)
        for p in range(k):
            byte = data[:, p:enc_len * inc:inc]
            hash[:, 0:byte.shape[1]] ^= byte
            hash = (hash * FNV_prime) & 0xFFFFFFFF

        return hash.astype(np.uint32)
    else:
        raise ValueError("data should be 2D or 3D array")


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def print_summary(conf, out=True):
    output = f"""
    {bcolors.OKBLUE}============================== test set:
    |   (dirs) = {conf['hap_dir'] if 'hap_dir' in conf else ""} {conf['pos_dir'], conf['gmap_dir'], conf['maf_dir'], conf['out_dir']}
    {bcolors.HEADER}============================== parameters:
    |   (chr) = {conf['chr_range']}
    |   mode = {'naive maf' if conf['mode'] == 1 else 'gmap' if conf['mode'] == 2 else 'cM+fixed' if conf['mode'] == 3 else 'cM+target' if conf['mode'] == 4 else 'random'}
    """
    output += f"|   sampling (target_len, len_seg_cM, step_cM) = {conf['target_len'], conf['len_seg_cM'], conf['step_cM']}"
    output += f"""
    |   (splits) = {conf['max_segs']}
    |   = total splits = {conf['max_num_trial']}
    |   encoding (k, inc) = {conf['k'], conf['inc']}
    |   LSH/Uhash = {conf['LSH_mode'], 'FNV' if conf['hash_mode'] < 2 else "uHash"}
    """
    if 'n' in conf:
        n = conf['n']
        key_space_size = conf['key_space_size']
        if 'actual_num_trial' not in conf:
            conf['actual_num_trial'] = conf['max_num_trial'] = sum(conf['max_segs'])
        output += f"""| hash (l, L) = {conf['l'], conf['L']}
|   num trials = {conf['actual_num_trial'] * conf['L']}"
|   actual/max # segs = {conf['actual_num_trial']}/{conf['max_num_trial']}"
|   key_space (ratio, size) = {conf['key_space_ratio'], key_space_size}"
|   (empty_cap, cap_each, cap_merged) =  {conf['trim_cap']}, {conf['cap_each']}, {conf['cap_merged']}"
"""
        eff = conf['cap_merged'] ** 2 * key_space_size / n ** 2
        output += f"\ntable_size = {key_space_size} efficiency = {eff * 100:.2f}%"
    output += f"{bcolors.ENDC}"
    if out:
        print(output)
    else:
        return output


def report_time(a, b, output=False):
    # global verbose
    start = time.time()
    try:
        if len(b) == 2:
            H = a(b[0], b[1])
        elif len(b) == 3:
            H = a(b[0], b[1], b[2])
        elif len(b) == 4:
            H = a(b[0], b[1], b[2], b[3])
        elif len(b) == 5:
            H = a(b[0], b[1], b[2], b[3], b[4])
        else:
            raise TypeError
    except TypeError:
        # print("handled")
        H = a(b)
    m, s = divmod(time.time() - start, 60)
    h, m = divmod(m, 60)
    s,m,h = int(round(s, 0)), int(round(m, 0)), int(round(h, 0))
    # if verbose:
    if output:
        print("Execution Time: " + "{0:02d}:{1:02d}:{2:02d}".format(h, m, s))
    else:
        print("Execution Time: " + "{0:02d}:{1:02d}:{2:02d}".format(h, m, s), file=sys.stderr)
    return H

def init_conf(conf, namespace, create_out=True):
    try :
        conf['n'] = namespace.n
        conf['param_dir'] = namespace.param
        conf['hap_dir'] = namespace.hap
    except AttributeError:
        conf['key_space_size'] = int(namespace.N)
        conf['pos_dir'] = namespace.pos
        conf['maf_dir'] = namespace.maf
        conf['gmap_dir'] = namespace.gmap
        conf['SNP_CNT_MIN'] = namespace.enclen
        # default to  half-overlapping segments
        assert conf['mode'] == 4
        conf['len_seg_cM'] = namespace.seglen  
        conf['step_cM'] = namespace.steplen
        conf['k'] = namespace.k
        # if 'step_cM' not in conf:
        #     conf['step_cM'] = conf['len_seg_cM'] / 2
    try:
        conf['L'] = namespace.L
    except AttributeError:
        conf['L'] = namespace.maxL

    conf['max_cMs'] = np.asarray(conf['max_cMs'])
    conf['max_segs'] = [0] + [int(a) for a in np.floor(((conf['max_cMs'] - conf['len_seg_cM'])/conf['step_cM']))]
    conf['out_dir'] = namespace.out
    conf['chr_cnt'] = len(conf['chr_range'])
    conf['chr_begin'] = conf['chr_range'].start
    conf['chr_end'] = conf['chr_range'].stop

    conf['max_num_trial'] = sum(conf['max_segs'][conf['chr_begin']:conf['chr_end']])

    if 'inc' not in conf:
        conf['inc'] = conf['k']

    # check if out_dir exists
    if not create_out and not Path(conf['out_dir']).exists():
        raise FileNotFoundError(f"{conf['out_dir']} does not exist")


KING_THRESHOLDS = np.asarray([1, .5 ** 1.5, .5 ** 2.5, .5 **3.5, .5 ** 4.5])
# for deg-4 tests use the following instead
# KING_THRESHOLDS = np.asarray([1, .5 ** 1.5, .5 ** 2.5, .5 **3.5, .5 ** 4.5, .5 ** 5.5])

def start_timer(info=print):
    start = datetime.now(pytz.timezone('US/Eastern'))
    info(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    info(f"start = {start}")
    return start

def stop_timer(start, info=print, strs=""):
    stop = datetime.now(pytz.timezone('US/Eastern'))
    total_time = str(stop - start)
    info(f"{strs} total = {total_time}, end = {stop}")
    info(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

def write_number(fname, M):
    with open(fname, "w") as f:
        f.write(f"{M}\n")

def jaccard_distance(a, b):
    return 1.0 - len(np.intersect1d(a, b)) / len(np.union1d(a, b))

def hamming_distance(a, b):
    return np.count_nonzero(a != b) / len(a)

def load_table(conf, party, tmp=""):
    conf = conf.copy()
    out_dir = conf['out_dir'] + f"tables{tmp}/party{party + 1}/"
    loaded = np.load(out_dir + f"ID_table.npz")
    N = loaded['N']
    B = loaded['B']
    id_table = loaded['ID_table']
    return id_table.reshape((N, B)), N, B