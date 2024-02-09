import numpy as np
import argparse

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=
                                   '''Python script to sketch the SNPs to get a subset for MHE

                                   Note: the boolean mask for the subset is stored as SNPs.npz file with the key 'snp_range' and the number of SNPs in num_SNP.txt.

                                   Please provide the same path to the MHE step and make sure both parties uses the same files.
                                   ''')
    parser.add_argument('-M', type=int, help='number of SNPs on which KING is computed', required=True)
    parser.add_argument('-s', type=float, help='subsampling rate (i.e. size(outputSNPs)/size(total))', required=True)
    parser.add_argument('-out', type=str, help='output directory', required=True)
    namespace = parser.parse_args()
    ratio = namespace.s
    save_path = namespace.out
    M = namespace.M
    print("=" * 51)
    print("*" * 2 + " " * 3 +  "SF-Relate Step 0b: Sample a subset of SNPs" + " " * 3 + "*" * 2)

    snp_len = int(ratio * M)
    print(f"  Sub-sampling rate = {ratio}\n  # of SNPs in subset = {snp_len}")
    snp_range = np.random.choice(M, (snp_len), replace=False)
    snp_range = np.sort(snp_range).astype(np.int32)
    np.savez(save_path + f"/SNPs.npz", snp_range = snp_range)
    with open(save_path + f"/num_SNP.txt", "w") as f:
        f.write(f"{snp_len}")
    print(f"** Output saved in {save_path}" + ("/" if save_path[-1] != "/" else ""))