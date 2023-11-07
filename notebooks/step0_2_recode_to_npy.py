#%%
from step0_1_split_data import *
from user_config import conf

%cd notebooks

if __name__ == "__main__":
    # Set test names and params here
    exec(param_from_bash())
    args = eval(ARG_STR)
    rid, rep = read_cmd_args(conf, args=args)
    k = 2
    ver = 0
    test_dir = conf_subset['test_dir'] = f"data/{k}party_{n}{'' if ver == 0 else 'v'+str(ver)}/"
    list_names = [test_dir + f'party{i}/list.txt' for i in range(k)]
    chr_range = range(1, 23)
    maf_LB = 0.01
    batch_splits = [0] + [4] * 17 + [2] * 3 + [1] * 2
    haps_dir = 'haps/'
    geno_dir = 'geno/'
    geno_flag = True
    haps_flag = True
    ref = "data/maf"
    suff = ''
    permuted = False
    chromo = 20
    i = 0

    # compute maf and export to a test case directory
    maf_fname = f"data/maf/chr{chromo}.txt"
    maf_dir = test_dir + 'maf/'
    run_command(['mkdir', f'{maf_dir}'])
    run_command(['mkdir', f'{maf_dir}{maf_LB}/'])
    run_command(['mkdir', f'data/pos/'])
    for i in range(k):
        run_command(['mkdir', f"{test_dir}party{i}/{haps_dir}/"])
        run_command(['mkdir', f"{test_dir}party{i}/{geno_dir}/"])
    start = start_timer()
    def save_maf(chromo):
        maf_fname = f"data/maf/chr{chromo}.txt"

        maf_data = pd.read_csv(maf_fname, sep="\s+", header=None)
        party_dir = f"{test_dir}party0/"
        bgen_name = f"{party_dir}bgen/chr{chromo}.bgen"
        bgen = open_bgen(bgen_name, verbose=False)
        # save list of positions
        np.savetxt(f"data/pos/chr{chromo}.txt", np.asarray(bgen.positions.copy()))

        maf_rsids = maf_data[0]
        # find the intersection of rsids
        locs = np.isin(maf_rsids, bgen.rsids)
        locs_bgen = np.isin(bgen.rsids, maf_rsids)
        # verify the order is the same
        assert np.all(maf_rsids[locs] == bgen.rsids[locs_bgen])
        maf_data_bgen = maf_data[locs]
        # get a boolean array of the snps in the bgen file,
        # such that the corresponding maf is between maf_LB and 1 - maf_LB
        maf_filt = (maf_data_bgen[5] >= maf_LB) & (maf_data_bgen[5] <= (1 - maf_LB))
        locs_bgen[locs_bgen] = maf_filt
    
        new_maf_fname = f"{maf_dir}/{maf_LB:.2f}/chr{chromo}.txt"

        np.save(f"{test_dir}/maf/{maf_LB:.2f}/chr{chromo}", locs_bgen)
        mafs = maf_data_bgen[maf_filt][5]
        np.savetxt(new_maf_fname, mafs)

        numsnps_fname = f"{maf_dir}/{maf_LB:.2f}/chr{chromo}_numsnps.txt"
        write_number(numsnps_fname, f"{len(mafs)}")
        
        bgen.close()
        return len(mafs)
    
    nums = Parallel(n_jobs=16)(delayed(save_maf)(chromo) for chromo in chr_range)
    M = sum(nums)
    print(f"using {M} snps for hashing")
    fname = f"{test_dir}/maf/{maf_LB:.2f}/all_chrs_numsnps.txt"
    write_number(fname, M)

    stop_timer(start)
    #%%
    # save all samples that are randomly chosen in .txt
    samples = load_samples(test_dir, k)
    for i in range(k):
        np.savetxt(f"{test_dir}/party{i}/samples.txt", samples[i], fmt='%d')

    if permuted:
        # randomly permute the samples if needed
        # not recommended in tests
        perm = []
        for i in range(k):
            perm.append(np.random.permutation(n))
            samples[i] = samples[i][perm[i]]
            np.savetxt(f"{test_dir}party{i}/perm.txt", perm[i], fmt='%d')

    # save geno and hap files to npy for faster reading in scripts
    start = start_timer()
    def recode_genomes(chromo, i):
        if not geno_flag and not haps_flag:
            return
        if geno_flag:
            filt_king_fname = f"{ref}/chr{chromo}.npy"
            filt_king = np.load(filt_king_fname)

        if haps_flag:
            filt_fname = f"{test_dir}/maf/{maf_LB:.2f}/chr{chromo}.npy"
            filt = np.load(filt_fname)

        def apply_filter(data, filt):
            ret = data[:, filt]
            assert ret.shape[1] == np.count_nonzero(filt)
            return ret

        print(f"Recoding {chromo, i}")
        if chromo > 20:
            batch_size = n
        elif chromo > 17:
            batch_size = n // 2
        else:
            batch_size = n // 4
        rep = int(np.ceil(n/batch_size))
        assert batch_splits[chromo] == rep
        batch_size = n // batch_splits[chromo]

        for j in range(rep):
            print(f"chromo {chromo}, batch {j} party {i}")
            batch_begin = j * batch_size
            batch_end = min((j + 1) * batch_size, n)
            party_dir = f"{test_dir}party{i}/"
            bgen_name = f"{party_dir}bgen/chr{chromo}.bgen"
            bgen = open_bgen(bgen_name, verbose=False)

            def alleles2id(alleles, m):
                # Map genotype probs into ACGT encoding
                return np.unique(np.stack(np.char.split(alleles, ",")), return_inverse=True)[1] \
                    .astype(np.int8).reshape((m, 2))

            tmp = alleles2id(bgen.allele_ids, bgen.nvariants)
            allele_ids = np.hstack((tmp, tmp))
            genos_raw = bgen.read((range(batch_begin, batch_end), None)).astype(np.bool_)

            if rep == 1:
                fname_hap = f"{party_dir}{haps_dir}chr{chromo}_maf{maf_LB:.2f}"
                fname_gen = f"{party_dir}{geno_dir}chr{chromo}_maf{maf_LB:.2f}"
            else:
                fname_hap = f"{party_dir}{haps_dir}chr{chromo}_maf{maf_LB:.2f}_{j}"
                fname_gen = f"{party_dir}{geno_dir}chr{chromo}_maf{maf_LB:.2f}_{j}"


            # parse haplotypes data
            def parse_haps(genos_raw, allele_ids):
                n, m, p2 = genos_raw.shape
                tmp = np.broadcast_to(allele_ids, (n, m, p2))
                a = np.transpose(tmp[genos_raw[:, :, 0:4]].reshape((n, m, 2)), (0, 2, 1))
                return a.reshape((2 * n, m))

            if haps_flag:
                haps = parse_haps(genos_raw, allele_ids)
                haps = apply_filter(haps, filt)
                np.save(fname_hap, haps)
            
            # parse genotypes data
            def parse_geno(genos_raw):
                return genos_raw[:, :, 1] + genos_raw[:, :, 3]
            
            if geno_flag:
                genotype_scores = parse_geno(genos_raw.astype(np.int8))
                genotype_scores = apply_filter(genotype_scores, filt_king)
                np.save(fname_gen, genotype_scores)

            bgen.close()

    # mapping are split into two parts because of memory usage
    for i in range(k):
        Parallel(n_jobs=3)(delayed(recode_genomes)(chromo, i) for chromo in (range(22, 3, -1)))
    for i in range(k):
        Parallel(n_jobs=3)(delayed(recode_genomes)(chromo, i) for chromo in (range(3, 0, -1)))

    stop_timer(start)

    # Stack all the batches together
    start = start_timer()
    import numpy as np
    def decode_batch(name, b):
        fname = f"{name}_{b}.npy"
        batch_data = np.load(fname)
        run_command(["rm", fname])
        return batch_data

    def piece_one_chromo(name, i, chromo, is_phased):
        def apply_perm(data):
            ret = np.zeros(data.shape, dtype=data.dtype)
            s = np.arange(n)
            if is_phased:
                ret[2 * s] = data[perm[i] * 2]
                ret[2 * s + 1] = data[perm[i] * 2 + 1]
            else:
                ret = data[perm[i]]
            return ret

        print(f"working on {i, chromo}")
        if batch_splits[chromo] > 1:
            whole = np.vstack([decode_batch(name.format(i, chromo, maf_LB), b) for b in range(batch_splits[chromo])])
        else:
            whole = np.load(name.format(i, chromo, maf_LB) + ".npy")
        
        if permuted:
            whole = apply_perm(whole)
        np.save(name.format(i, chromo, maf_LB) + suff, whole)

    def piece_together(name, is_phased):
        Parallel(n_jobs=8)(delayed(piece_one_chromo)(name, i, chromo, is_phased) for i in range(k) for chromo in chr_range)

    if haps_flag:
        fname_hap = test_dir + "party{}/" + haps_dir + "chr{}_maf{:.2f}"
        piece_together(fname_hap, True)

    if geno_flag:
        fname_gen = test_dir + "party{}/" + geno_dir + "chr{}_maf{:.2f}"
        piece_together(fname_gen, False)

    stop_timer(start)

    def piece_all_chromo(i):
        all_geno = [] 
        fname_all_geno = f"{test_dir}party{i}/{geno_dir}all_chrs_maf{maf_LB:.2f}{suff}"
        for chromo in chr_range:
            print(f"reading {chromo}", end='...')
            all_geno.append(np.load(test_dir + f"party{i}/{geno_dir}chr{chromo}_maf{maf_LB:.2f}{suff}.npy"))
            print("done")
        whole = np.hstack(all_geno)
        print("HI")
        np.save(fname_all_geno, whole)
        run_command(['mkdir', f"{test_dir}{geno_dir}/party{i + 1}", '-p'])
        # convert the npy files to binary files for MHE
        whole.tofile(f"{test_dir}{geno_dir}/party{i + 1}/all_chrs_maf{maf_LB:.2f}.bin")
        print("BYE")

    for i in range(k):
        piece_all_chromo(i)