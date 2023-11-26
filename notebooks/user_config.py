import numpy as np

conf = {
    "test_dir": "data/2party_{}/",
    "ground_truth": 'ground_truths/KING.dat',
    "snp_threshold": 0.01,
    "cap_merged": 1, "cap_each": 1,
    "n": None,
    "w": 8,
    "len_seg": 600,
    "k": 8,
    "ver": 38, 
    "target_len": 80,
    "mode": 4, 
    "LSH_mode": "minHash",
    # "LSH_mode": "Hamming",
    "l": 3,
    "L": 6,
    "len_seg_cM": 8,
    "key_space_ratio": 64,
    # only mode 1 is supported in release
    "hash_mode": 1, 
    # 0 = FNV-1, 1 = FNV-1a, 2 = uhash from L^2LSH

    # max_cMs found in UKB dataset
    "max_cMs": np.asarray([286.2789055501355, 268.8318954936778, 223.25963499194336, 214.54387131588663, 204.08614941635688, 192.032233, 187.15981416488336, 168.00070572595456, 166.353007141962, 180.93048929289142, 158.21844015613382, 174.679006, 125.704447, 120.202182, 141.345171, 134.036757, 128.48575436470588, 117.70657414812156, 107.73365132875236, 108.215478, 62.785274288801574, 74.1058950938716]),
    "bucket_size": None,
    "chr_range": range(1, 23),
    'maf_LB': 0.01,
    "verbose": False,
    "n_party": 2,
    "bgen_name_template": 'chr{}.bgen',
    "cache_name_template": "chr{}_{}.npz",
    "n_job": 16,
    "out_dir": '{}/table/{}/',
    "out_dir_tmp": '{}/table/{}/',
}

