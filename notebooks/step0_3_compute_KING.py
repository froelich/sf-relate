#%%
from step0_1_split_data import *
from my_util import write_king_file

if __name__ == "__main__":
    # Set test names and params here
    exec(param_from_bash())
    args = eval(ARG_STR)
    rid, rep = read_cmd_args(conf, args=args)
    k = 2
    ver = 0
    start = start_timer()
    out_dir, ratio = f'data/2party_{n}/' + (f'v{ver}' if ver != 0 else ''), 0.01
    chr_range = range(22, 0, -1)
    ratio = 0.01
    batch_splits = [0] + [4] * 17 + [2] * 3 + [1] * 2
    read_suff = ''
    out_suff = ''
    permuted = False 
    batch_size = min(25000, n)
    batch_cnt = int(np.ceil(n / batch_size))
            
    # split all chromo into batches so that we can compute matrix product in memory
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
        kings = np.zeros((batch_size, n), dtype=np.float32)

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
        
    Parallel(n_jobs=1)(delayed(compute_king)(p1, None) for p1 in range(0, batch_cnt))
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

# %%
