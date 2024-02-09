import numpy as np
from default_config import conf
import argparse
import pickle
from my_util import *
from kmers import *

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=
                                   '''Python script to sample shared random parameters for hashing.

                                   Please make sure the parameters are the same on the two parties and provide the path to step_1_hashing.py
                                   ''')
    parser.add_argument('-N', type=str, help='size of hash tables (recommend: 64 * total number of individuals on both parties)', required=True)
    parser.add_argument('-pos', type=str, help='the directory storing the bp positions. Files are named chr{i}.txt for i = 1..22. They contain the list of physical positions of each SNPs on the haplotypes (one number per line).', required=True)
    parser.add_argument('-gmap', type=str, help='the directory storing the genetic maps. Files are named chr{i}.gmap.gz (in gz format) for i = 1..22. The first line of the file contains `pos\tchr\tcM`, and each following line contains the bp location, the chromosome ID and the corresponding genetic location (separated by tabs). One can retrieve these files from [shapeit4](https://github.com/odelaneau/shapeit4/tree/master/maps) or other public resources, but should be careful to make sure the correct genome build locations is used.', required=True)
    parser.add_argument('-maf', type=str, help='the directory storing the maf files. Files are named chr{i}.txt for i = 1..22. Each line in the file stores a floating point number denoting the MAF of that base pair.', required=True)
    parser.add_argument('-out', type=str, help='output directory to store the parameters', required=True)
    parser.add_argument('-enclen', type=int, help='the number of snps in each encoded split haplotype segment (default: 80).', default=80)
    parser.add_argument('-seglen', type=float, help='centi-Morgan length of each split haplotype segment (default: 8.0)', default=8.0)
    parser.add_argument('-steplen', type=float, help='centi-Morgan spacing between the beginning of each split haplotype segment (default: 4.0)', default=4.0)
    parser.add_argument('-k', type=int, help='number of SNPs in each kSNP token for hashing (default: 8)', default=8)
    parser.add_argument('-l', type=int, help=' number of hash tokens to construct every hash index (default: 4)', default=4)
    parser.add_argument('-maxL', type=int, help='max number of repetitive hashing; increase and retry if table saturation is low     (default: 6; should be larger than the argument to step_1_hashing.py)', default=6)
    namespace = parser.parse_args()
    init_conf(conf, namespace)

    chr_begin = conf['chr_range'].start
    chr_end = conf['chr_range'].stop
    max_segs = conf['max_segs']

    print("=" * 50)
    print("*" * 2 + " " * 3 +  "SF-Relate Step 0a: Sample Hashing Randomness" + " " * 3 + "*" * 2)

    L = namespace.maxL
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
    Path(conf['out_dir']).mkdir(parents=True, exist_ok=True)

    # optionally random order for merging the segments
    # in paper, the default preference is to choose from lowest to highest indices
    merge_order = np.full((L), None, dtype=np.object_)
    for i in range(L):
        if i % 2 == 0:
            merge_order[i] = np.random.permutation(max_num_trial)
        else:
            merge_order[i] = merge_order[i-1][::-1]  


    np.savez(conf['out_dir']+ f"/LSH_params", indices=indices, h_keys=h_keys, merge_order=merge_order)
    # pickle file for later use
    # change ndarray to list for dict()
    conf['max_cMs'] = conf['max_cMs'].tolist()
    pickle.dump(conf, open(conf['out_dir']+ f"/LSH_params_conf.pkl", "wb"))
    with open(conf['out_dir']+ f"/LSH_params_conf.txt", "w") as f:
        f.write(print_summary(conf, out=False))

    print(f"Output saved in {conf['out_dir']}" + ("/" if conf['out_dir'][-1] != "/" else ""))