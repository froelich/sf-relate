import numpy as np
import time
import sys
import pandas as pd
from pathlib import Path
from datetime import datetime
import subprocess
import pytz

def run_command(command, silent = False, cwd=None):
    result = subprocess.run(command, capture_output=True,text=True, cwd=cwd)
    if silent:
        return
    if result.stderr != '':
        print(result.stderr)
    if result.returncode != 0:
        print(result.returncode)
        print(result.stdout)
    return result


def king_formula(x, y, hx=None, hy=None):
    sub = x - y
    z = np.nansum(sub * sub)
    if hx is None:
        hx = np.count_nonzero(x == 1)
    else:
        hx = np.count_nonzero(hx ==1)
    if hy is None:
        hy = np.count_nonzero(y == 1)
    else:
        hy = np.count_nonzero(hy ==1)
    return 1 / 2 - 1 / 4 * z / min(hx, hy)


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
    |   (dirs) = {conf['test_dir'], conf['out_dir']}
    |   (ground_truth) = {conf['ground_truth']}
    {bcolors.HEADER}============================== parameters:
    |   (chr) = {conf['chr_range']}
    |   mode = {'naive maf' if conf['mode'] == 1 else 'gmap' if conf['mode'] == 2 else 'cM+fixed' if conf['mode'] == 3 else 'cM+target' if conf['mode'] == 4 else 'random'}
    """
    if conf['mode'] == 2:
        output += f"|   sampling (w_cM, len_seg_cM, step_cM) = {conf['w_cM'], conf['len_seg_cM'], conf['step_cM']}" 
    elif conf['mode'] == 1:
        output += f"|   sampling (w, len_seg, step) = {conf['w'], conf['len_seg'], conf['step']}"
    elif conf['mode'] == 3:
        output += f"|   sampling (w, len_seg_cM, step_cM) = {conf['w'], conf['len_seg_cM'], conf['step_cM']}"
    elif conf['mode'] == 4:
        output += f"|   sampling (target_len, len_seg_cM, step_cM) = {conf['target_len'], conf['len_seg_cM'], conf['step_cM']}"
    output += f"""
    |   (splits) = {conf['max_segs']}
    |   = total splits = {conf['max_num_trial']}
    |   encoding (k, inc) = {conf['k'], conf['inc']}
    |   LSH/Uhash = {conf['LSH_mode'], 'FNV' if conf['hash_mode'] < 2 else "uHash"}
    """
    key_space_size = conf['key_space_size']
    n = conf['n']
    if "no_hashing" not in conf:
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

def init_conf(conf, n, hash_id=0, create_out=True):
    conf['n'] = n
    conf['test_dir'] = conf['test_dir'].format(n)
    conf['chr_cnt'] = len(conf['chr_range'])
    conf['chr_begin'] = conf['chr_range'].start
    conf['chr_end'] = conf['chr_range'].stop
    conf['out_dir'] = conf['out_dir_tmp'].format(conf['test_dir'], hash_id)
    conf['kmer_dir'] = conf['out_dir']+'kmer/'
    # default to  half-overlapping segments

    if conf['mode'] >= 2:
        if 'step_cM' not in conf:
            conf['step_cM'] = conf['len_seg_cM'] / 2
        conf['max_segs'] = [0] + [int(a) for a in np.floor(((conf['max_cMs'] - conf['len_seg_cM'])/conf['step_cM']))]
    # if conf['mode'] == 1 or conf['mode'] == 3:
    if conf['mode'] == 1:
        if 'step' not in conf:
            conf['step'] = int(conf['len_seg'] / 2)
        conf['max_segs'] = [0]
        for chrom in range(1, 23):
            # read numbers of snps in each chromosome in the file f"chr{chrom}_maf0.01_numsnps.txt"
            with open(f"{conf['test_dir']}/maf/{conf['maf_LB']:.2f}/chr{chrom}_numsnps.txt", 'r') as f:
                num_snps = int(f.read())
            conf['max_segs'].append(int(np.ceil(((num_snps - conf['len_seg'])/conf['step']))))
         

    conf['max_num_trial'] = sum(conf['max_segs'][conf['chr_begin']:conf['chr_end']])

    if 'inc' not in conf:
        conf['inc'] = conf['k']

    conf['key_space_size'] = int(conf['n'] * 2 * conf['key_space_ratio'])
    # check if out_dir exists
    if not create_out and not Path(conf['out_dir']).exists():
        raise FileNotFoundError(f"{conf['out_dir']} does not exist")


import argparse
def read_cmd_args(conf, args=None):
    parser = argparse.ArgumentParser( description='Build two LSH tables using the parameters')

    parser.add_argument('-h_mode', type=int, help='universal hash mode: 0 for FNV-1, 1 for FNV-1a, 2 for random linear combinations')

    parser.add_argument('-n', type=int, help='input size', required=True)
    parser.add_argument('-mode', type=int, default=4, help='encoding mode')
    parser.add_argument('-len_seg', type=int, help='len of segment')
    parser.add_argument('-k', type=int, default=8, help='kmer length')
    parser.add_argument('-len_seg_cM', default=7, type=float, help='len of segment (in cM)')
    parser.add_argument('-w', type=int, help='width of sampling window')
    parser.add_argument('-l', type=int, default=4, help='AND composition of LSH')
    parser.add_argument('-table_ratio', type=float, default=64, help='table ratio')
    parser.add_argument('-L', type=int, default=4, help='num of repetition')
    parser.add_argument('-target_len', type=int, default=80, help='target length of kmer vectors')
    parser.add_argument('-LSH_mode', type=str, default='Hamming', help='LSH mode: Hamming or minHash')
    parser.add_argument('-hid', type=str, default='0', help='test id: out_dir.')
    parser.add_argument('-rid', type=str, default='0', help='king id: rid') 
    parser.add_argument('-rep', type=int, default='1', help='num rep of the test')
    parser.add_argument('-snpmin', type=int, default='50', help='num min snps in each segment')
    parser.add_argument('-trim_cap', type=int, default='5', help='removal threshold for large buckets')

    if args is not None:
        namespace = parser.parse_args(args=args)
    else:
        namespace = parser.parse_args()
    conf['test_dir'] = conf['test_dir'].format(namespace.n)
    mode = conf['mode'] = namespace.mode
    if mode == 1:
        conf['len_seg'] = namespace.len_seg
    elif mode == 3:
        conf['len_seg_cM'] = namespace.len_seg_cM
        conf['w'] = namespace.w
    elif mode == 4:
        conf['len_seg_cM'] = namespace.len_seg_cM
        conf['step_cM'] = namespace.len_seg_cM  / 2
        conf['target_len'] = namespace.target_len

    conf['hid'] = namespace.hid
    conf['n'] = namespace.n
    conf['hash_mode'] = namespace.h_mode
    conf['k'] = namespace.k
    conf['L'] = namespace.L
    conf['LSH_mode'] = namespace.LSH_mode
    conf['l'] = namespace.l
    conf['key_space_ratio'] = namespace.table_ratio 
    conf['key_space_size'] = int(namespace.table_ratio * namespace.n * 2)
    conf['SNP_CNT_MIN'] = namespace.snpmin
    conf['trim_cap'] = namespace.trim_cap
    return namespace.rid, namespace.rep


KING_THRESHOLDS = np.asarray([1, .5 ** 1.5, .5 ** 2.5, .5 **3.5, .5 ** 4.5])
# for deg-4 tests use the following instead
# KING_THRESHOLDS = np.asarray([1, .5 ** 1.5, .5 ** 2.5, .5 **3.5, .5 ** 4.5, .5 ** 5.5])
KING_THRESHOLDS_REV = np.asarray(KING_THRESHOLDS[::-1])
def relativedegree(kinship):
    return len(KING_THRESHOLDS) - 1 - np.searchsorted(KING_THRESHOLDS_REV, kinship) 


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

def read_king_file_raw(fname):
    related = pd.read_csv(fname, sep='\t',
                          dtype={'P0':np.int8, 'ID0':np.int32, 'P1':np.int8, 'ID1':np.int32, 'Kinship':np.float32, 'deg':np.int8})
    return related

def load_samples(test_dir, k, chromo=21):
    from bgen_reader import open_bgen
    samples = [[], []]
    for i in range(k):
        with open_bgen(f"{test_dir}party{i}/bgen/chr{chromo}.bgen", verbose=False) as bgen:
            samples[i] = np.asarray(bgen.samples.copy(), dtype=np.int64)
    return samples

def load_samples_from_raw():
    from bgen_reader import open_bgen
    with open_bgen(f"data/raw/chr21.bgen", verbose=False) as bgen:
        samples = np.asarray(bgen.samples.copy(), dtype=np.int64)
    return samples

def read_king_file(fname):
    related = read_king_file_raw(fname)
    return related.groupby(lambda idx: related.loc[idx]['deg'])

#================================================== 
# snp qc files are only used when computing the set of SNPs on which to compute KING
def read_snp_qc_raw():
    snp_qc = pd.read_csv("data/maf/ukb_snp_qc.txt", sep='\s+', usecols=['rs_id', 'chromosome', 'PC1_loading', 'PC2_loading', 'PC3_loading', 'in_Relatedness', 'position', 'in_Phasing_Input'])
    snp_qc_autosome = snp_qc[snp_qc['chromosome'] <= 22]
    return snp_qc_autosome

def read_snp_qc():
    snp_qc = read_snp_qc_raw()
    return snp_qc.groupby(lambda idx: snp_qc.loc[idx]['chromosome'])
# end of snp qc files
#================================================== 

def write_number(fname, M):
    with open(fname, "w") as f:
        f.write(f"{M}\n")

def jaccard_distance(a, b):
    return 1.0 - len(np.intersect1d(a, b)) / len(np.union1d(a, b))

def hamming_distance(a, b):
    return np.count_nonzero(a != b) / len(a)


def write_king_file(fname, data):
    file = open(fname, 'w')
    file.write(f"P0\tID0\tP1\tID1\tKinship\tdeg\n")
    data.sort(key=lambda x:x[-1])
    file.write("\n".join(["\t".join([str(x) for x in row]) for row in data]))
    file.close()

def sample_from_thr(all_king, data, thr, target=5000, up_thr=None, remaining_pairs=None, remaining_kinship=None, verbose=False):
    if up_thr is None:
        up_thr = 1
    if remaining_pairs is None:
        rst_pairs = np.argwhere((all_king >= thr) & (all_king < up_thr))
        print(f"between {thr} and {up_thr} samples = {len(rst_pairs)}")
        if target is None or target > len(rst_pairs):
            target = len(rst_pairs)
        sampled = np.random.choice(len(rst_pairs), target, replace=False)
        rst_pairs = rst_pairs[sampled]
        rst_kinship = all_king[rst_pairs[:,0], rst_pairs[:,1]]
        remaining_pairs = np.argwhere(all_king >= up_thr)
        remaining_kinship = all_king[(all_king >= up_thr)]
    else:
        rst_pairs = remaining_pairs[remaining_kinship < up_thr]
        print(f"between {thr} and {up_thr} samples = {len(rst_pairs)}")
        if target is None or target > len(rst_pairs):
            target = len(rst_pairs)
        sampled = np.random.choice(len(rst_pairs), target, replace=False)
        rst_pairs = rst_pairs[sampled]
        rst_kinship = remaining_kinship[remaining_kinship < up_thr][sampled]
        remaining_pairs = remaining_pairs[remaining_kinship >= up_thr]
        remaining_kinship = remaining_kinship[remaining_kinship >= up_thr]

    if all_king is not None:
        del all_king
        all_king = None

    print(f"adding {len(rst_pairs)} entries")
    rel_deg = [0] * (len(KING_THRESHOLDS))
    for i in range(len(rst_pairs)):
        a, b = rst_pairs[i]
        kin = rst_kinship[i]
        deg = relativedegree(kin)
        rel_deg[deg] += 1
        if deg == -1:
            deg = len(KING_THRESHOLDS)
        data.append((0, a, 1, b, kin, deg))

    if verbose:
        for i in range(len(rel_deg)):
            print(f"have {rel_deg[i]} many deg {i} individuals in subset.")

    return remaining_pairs, remaining_kinship


def param_from_bash():
    qq = '''
minHash="minHash"
exact="exact"
Hamming="Hamming"
'''
    qq += ''.join(open('param.sh').readlines()[:-3])
    qq += "\ntid=f\'"+ open('param.sh').readlines()[-1].replace("$", "").split("=")[1].strip() + "\'"
    # print(qq)
    return qq

ARG_STR = 'f"-h_mode 1 -n {n} -mode {samp_mode} -k {k} -l {l} -table_ratio 64 -L {L} -len_seg_cM {len_seg_cM} -target_len {target_len} -LSH_mode {LSH} -rid {rid} -rep {rep} -hid {tid}".split(' ')'


#==================================================================================================== 
# the following are only used for post-analysis of experiments
# not used in the main experiment
def find_max_king(conf, half_mode=False):
    tmp = "_half" if half_mode else ""
    try:
        max_kings = np.load(f"plots/data/max_kings{tmp}.npy")
        degrees = np.load(f"plots/data/degrees{tmp}.npy")
    except FileNotFoundError:
        all_king = None
        if not half_mode:
            all_king = np.load(f"{conf['test_dir']}/ground_truths/all_king_maf0.01.npy")
        else:
            all_king = np.load(f"{conf['test_dir']}/ground_truths/all_king_maf0.01.npy")[:conf['n']//2, :conf['n']//2]
        # find max kinship per row/column
        max_kings = np.asarray([np.max(all_king, axis = 1- i) for i in range(2)])
        degrees = np.asarray([relativedegree(max_king) for max_king in max_kings])
        np.save(f"plots/data/max_kings{tmp}.npy", max_kings)
        np.save(f"plots/data/degrees{tmp}.npy", degrees)
    return max_kings, degrees

def get_data(conf, test, numBlock=13):
    ob_range = range(0, numBlock)
    IDs_parties = []
    scores_parties = []
    for party_id in range(2):
        IDs = []
        scores = []
        for i in ob_range:
            with open(f"{conf['test_dir']}/out/{test}/party{party_id + 1}/{i}_party{party_id + 1}.csv") as f:
                i = 0
                for line in f:
                    i += 1
                    if i == 1:
                        continue
                    a, b = line.split(',')
                    IDs.append(int(a))
                    b = b.replace('i', 'j')
                    scores.append(complex(b[1:-2]))
        IDs_parties.append(np.asarray(IDs))
        scores_parties.append(np.asarray(scores))
    return IDs_parties, scores_parties

def get_classes(IDs_parties, scores_parties, max_kings, degrees):
    avg_precision = 0
    avg_recall = 0
    watches = []
    caughts = []
    misses = []
    extras = []
    for party_id in range(2):
        IDs = np.asarray(IDs_parties[party_id])
        scores = np.asarray(scores_parties[party_id])
        real_scores = np.real(scores)
        watch = IDs[np.round(real_scores).astype(int) == 0]
        watches.append(watch)
        print(f"Party {party_id}")
        ground_truth = np.unique(np.nonzero(degrees[party_id] <= 3)[0])
        # ground_truth = np.unique(related[f'ID{party_id}'])
        good = np.intersect1d(ground_truth, watch)
        missed = ground_truth[~np.isin(ground_truth, watch)]
        ttl = len(np.unique(watch))
        extra = watch[~np.isin(watch, ground_truth)]
        caughts.append(good)
        misses.append(missed)
        extras.append(extra)
        recall = len(good)/len(ground_truth)
        precision = len(good)/ttl
        avg_recall += recall
        avg_precision += precision
        print(f"recall = {recall}")
        print(f"precision = {precision}")

    avg_recall /= 2
    avg_precision /= 2

    print(f"============ Overall ==========")
    print(f"recall = {avg_recall}")
    print(f"precision = {avg_precision}")

    # find_average discovery rate per class
    def report_kinship(list, i, label, good='error'):
        global xmin, xmax, ymin, ymax
        king_list = max_kings[i][list[i]]
        degree = degree_class[i][list[i]]
        return pd.DataFrame({"Kinship":king_list, "Deg":degree, "class":label, "Status": good, 'party' : i, 'ID' : list[i]})
        
    degree_class = degrees.astype(str)
    degree_class[degrees > 3] = 'unrelated'
    to_cat = []
    for i in range(2):
        to_cat += [report_kinship(extras, i, "found"), report_kinship(misses, i, "missed")]
        to_cat.append(report_kinship(caughts, i, "found", "ok"))
    data = pd.concat(to_cat)

    return data, caughts, missed, extras

def get_suf_ith(i):
    suf = ['st', 'nd', 'rd']
    try:
        return suf[i - 1]
    except IndexError:
        return 'th'
# end of analysis scripts
#==================================================================================================== 



def load_table(conf, party, tmp=""):
    conf = conf.copy()
    out_dir = conf['out_dir'] + f"tables{tmp}/party{party + 1}/"
    loaded = np.load(out_dir + f"ID_table.npz")
    N = loaded['N']
    B = loaded['B']
    id_table = loaded['ID_table']
    return id_table.reshape((N, B)), N, B


def load_tables(conf, tmp):
    bb = []
    N = B = 0
    for i in range(conf['n_party']):
        tbb, tN, tB = load_table(conf, i, tmp)
        if N != 0: assert (N == int(tN))
        if B != 0: assert (B == int(tB))
        N = int(tN)
        B = int(tB)
        bb.append(tbb)
    return bb, N, B
