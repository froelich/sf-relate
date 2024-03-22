import numpy as np
from default_config import conf
import argparse
import pickle
from kmers import *
from my_util import *
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
import tomlkit
import os 

if __name__ == "__main__":
    # read environment variables in python
    parser = argparse.ArgumentParser(description='Sample Hashing Randomness')
    parser.add_argument('-PARTY', type=int, help='party id', required=True)
    parser.add_argument('-FOLDER', type=str, help='path to the configuration folder', required=True)
    args = parser.parse_args()
    party_id = args.PARTY
    FOLDER = args.FOLDER
    
    args_local = tomlkit.load(open(f"{FOLDER}/configLocal.Party{party_id}.toml"))
    args_global = tomlkit.load(open(f"{FOLDER}/configGlobal.toml"))
    args_local['param'] = args_local['shared_param_dir']
    args_local['out'] = args_local['shared_param_dir']
    args_local['hap'] = args_local['haps_dir']
    args_local['pos'] = args_local['pos_dir']
    args_local['gmap'] = args_local['gmap_dir']
    args_local['pos'] = args_local['pos_dir']
    args_local['maf'] = args_local['maf_dir']
    namespace = argparse.Namespace()
    namespace.__dict__.update(args_global)
    namespace.__dict__.update(args_local)

    with open(args_global['shared_keys_path'] + "/seed.bin", "rb") as f:
        seed = int.from_bytes(f.read(), 'big')
    print(seed)
    gen = np.random.default_rng(seed)

    # find out which line in pvar_file start with #CHROM
    psam, snps, pvar = read_pgen_metadata(args_local)
    n, m = len(psam), len(pvar)
    namespace.__dict__.update({'n': n, 'm': m})
    init_conf(conf, namespace)
    conf['L'] = namespace.maxL

    # =============================================================================
    # Step 0a: Sample Hashing Randomness
    chr_begin = conf['chr_range'].start
    chr_end = conf['chr_range'].stop
    max_segs = conf['max_segs']

    print("=" * 50)
    print("=" * 20  + "Shared seed: ", seed)
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
                        h_keys[tid, trial_id] = ([(gen.integers(0, kmer_enc_len, dtype=np.int32),) for _ in range(conf['l'])],
                                                 [gen.integers(1, coef_uHash_range, dtype=np.uint32)
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
            merge_order[i] = gen.permutation(max_num_trial)
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

    # =============================================================================
    # Step 0b: Sample a subset of SNPs

    ratio = args_global['s']
    save_path = args_local['sketched_snps_dir']
    M = len(snps)
    print("=" * 51)
    print("*" * 2 + " " * 3 +  "SF-Relate Step 0b: Sample a subset of SNPs" + " " * 3 + "*" * 2)

    snp_len = int(ratio * M)
    print(f"  Sub-sampling rate = {ratio}\n  # of SNPs in subset = {snp_len}")
    snp_range = gen.choice(M, (snp_len), replace=False)
    snp_range = np.sort(snp_range).astype(np.int32)
    np.savez(save_path + f"/SNPs.npz", snp_range = snp_range)
    with open(save_path + f"/num_SNP.txt", "w") as f:
        f.write(f"{snp_len}")
    print(f"** Output saved in {save_path}" + ("/" if save_path[-1] != "/" else ""))