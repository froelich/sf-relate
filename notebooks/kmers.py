from my_util import *

class HaplotypeArray:
    # An object representing $n$ split haplotype/subchromosome/segment
    # storing kmer/kSNP/k-shingles in vectors
    def __init__(self, conf, chr_id, seg_id=None, cM=None, permuted=False, beg=None, end=None):
        if seg_id is not None:
            assert (beg is not None and end is not None) or cM is not None
        self.chr_id = chr_id
        self.conf = conf
        self.permuted=permuted
        
        num_segs = conf['max_segs'][chr_id]
        if seg_id is None:
            self.beg = 0
            self.end = -1
            return

        if seg_id == num_segs:
            raise ValueError(f"Segment {seg_id} does not exist for chromosome {chr_id}")
        if beg is not None and end is not None:
            self.beg = beg
            self.end = end
        else:
            self.beg, self.end = HaplotypeArray.get_boundaries(conf, seg_id, cM) 

        if cM is not None:
            self.cM = self.get_seg(cM)

    def get_cM(self):
        return self.cM

    @staticmethod
    def get_boundaries(conf, seg_id, cM):
        step_cM = conf['step_cM']
        len_seg_cM = conf['len_seg_cM']
        cM_begin = seg_id * step_cM
        cM_end = cM_begin + len_seg_cM
        return HaplotypeArray.cM_to_snp_indices(conf, cM, cM_begin, cM_end)

    @staticmethod
    def read_cM(conf, chr_range, plot=False):
        from matplotlib import pyplot as plt
        cMs = [None] * chr_range.start
        if plot:
            raise NotImplementedError("debug mode is disabled")
        for chrom in chr_range:
            pos = HaplotypeArray.read_pos(conf, chrom)
            maps = HaplotypeArray.read_gmap(conf, chrom)
            cM = maps['cM'].to_numpy()
            pos_map = maps['pos'].to_numpy()
            pos = pos[pos < pos_map[-1]]
            begins = np.searchsorted(maps['pos'], pos)
            bgen_cMs = cM[begins - 1] + (pos - pos_map[begins - 1]) * (cM[begins] - cM[begins - 1]) / (pos_map[begins] - pos_map[begins - 1])
            cMs.append(bgen_cMs)
        return  cMs

    @staticmethod
    def read_gmap(conf, chr_id):
        ver = conf['ver']
        gmap_fname = f'data/maps/chr{chr_id}.b{ver}.gmap.gz'
        return pd.read_csv(gmap_fname, sep="\t")

    @staticmethod
    def read_pos(conf, chr_id):
        pos = np.loadtxt(f"{conf['test_dir']}/pos/0.01/chr{chr_id}.txt")
        return pos

    @staticmethod
    def cM_to_snp_indices(conf, cM, cM_begin, cM_end):
        # get the min and max indices of the segment
        # which is the first index of the first SNP with cM >= cM_begin
        # and the last index of the last SNP with cM < cM_end
        idx_begin = np.searchsorted(cM, cM_begin)
        idx_end = np.searchsorted(cM, cM_end)
        if cM_end > cM[-1]:
            raise ValueError(f"cM_end {cM_end} is greater than the ending cM {cM[-1]}")

        if cM[idx_end] == cM_end:
            idx_end += 1
        
        SNP_CNT_MIN = conf['SNP_CNT_MIN']
        if idx_end - idx_begin < SNP_CNT_MIN:
            raise ValueError(f"Segment is too short. {idx_end - idx_begin} SNPs")

        return idx_begin, idx_end


    def get_seg(self, a):
        if len(a.shape) == 1:
            return a[self.beg:self.end]
        else:
            return a[:, self.beg:self.end]

    def read_mafs(self):
        maf_fname = f"{self.conf['test_dir']}../maf/{self.conf['maf_LB']:.2f}/chr{self.chr_id}.txt"
        maf = pd.read_csv(maf_fname, sep=" ", header=None).to_numpy().flatten()
        return self.get_seg(maf)

    def read_haps(self):
        npy_path = self.conf['test_dir']
        if self.permuted:
            npy_file = f"{npy_path}/haps/chr{self.chr_id}_maf0.01_p.npy" 
        else:
            npy_file = f"{npy_path}/haps/chr{self.chr_id}_maf0.01.npy" 
        x = np.load(npy_file, mmap_mode='r')
        return self.get_seg(x)
    

class KmerBuilder:
    def __init__(self, conf):
        self.conf = conf

    def sample_indices(self, mafs, cM=None, cM_begin=None):
        '''
        sample indices according to MAFs
        mode 1 - 3 are removed
        mode 4: sample according to MAFs from some varying-sized windows, such that the total number of SNPs is target_len
        mafs: a 1D array of MAFs
        cM: a 1D array of cM positions
        cM_begin: the beginning of the segment in cM
        '''
        m = len(mafs)
        mode = self.conf['mode']

        def enough_to_cover(m, w):
            return m // w + (m % w != 0)


        def sample_using_maf(m, w, target_len=None):
            # sample according to MAFs
            rnds = np.random.random_sample(enough_to_cover(m, w))
            if target_len is None:
                indices = np.zeros(len(rnds), int)
            else:
                indices = np.zeros(target_len, int)
            for i in range(len(indices)):
                if i == len(indices) - 1:
                    # use the rest of the SNPs
                    sect = mafs[i * w:]
                else:
                    sect = mafs[i * w: (i + 1) * w]
                prob_sum = sum(sect[sect >= self.conf['snp_threshold']])
                if prob_sum == 0:
                    indices[i] = i * w
                    continue
                dpf = 0
                # should increse w if it is the last window
                for j in range(len(sect)):
                    if i * w + j >= len(mafs):
                        break
                    if mafs[i * w + j] >= self.conf['snp_threshold']:
                        dpf += mafs[i * w + j] / prob_sum
                    if rnds[i] < dpf:
                        indices[i] = i * w + j
                        break
            return np.asarray(indices)

        if mode != 4:
            raise NotImplementedError("mode other than 4 are removed for simplicity")

        # for this mode
        # sample according to MAFs from some varying-sized windows, such that the total number of SNPs is target_len
        target_len = self.conf['target_len']
        indices = sample_using_maf(m, m // target_len, target_len=target_len)
        return np.asarray(indices)

    def window_project(self, hap_array, indices=None, subset=None):
        if subset is not None:
            # find the correct subsets 
            # the subset of rows are 2 * subset and 2 * subset + 1
            KmerBuilder.find_rows(subset)
        else:
            rows = np.arange(2 * hap_array.conf['n'])
        if indices is not None:
            if subset is None:
                return hap_array.read_haps()[:, indices]
            else:
                return hap_array.read_haps()[np.ix_(rows, indices)]
        else:
            if subset is None:
                return hap_array.read_haps()[:, self.sample_indices(hap_array.read_mafs())]
            else:
                return hap_array.read_haps()[np.ix_(rows, self.sample_indices(hap_array.read_mafs()))]

    @staticmethod
    def find_rows(subset):
        n = len(subset)
        rows = np.zeros(2 * n, dtype=np.int32)
        rows[2 * np.arange(n)] = 2 * subset
        rows[2 * np.arange(n) + 1] = 2 * subset + 1
        return rows


    def get_kmers(self, hap_array, rand_seed=0, indices=None, debug=False, h_key=None, subset=None):
        '''
        :return:
         set of kmers (stored by their 32-bit hash)
         Where it is stored as a1 a2 b1 b2 where a1 a2 are the two haplotypes of the same sample
        '''
        hash_mode = self.conf['hash_mode']
        data_set = self.window_project(hap_array, indices, subset=subset)
        inc = self.conf['inc']
        k = self.conf['k']
        assert (inc >= 1)
        if not debug:
            return apply_uhash(data_set, hash_mode, k, inc, rand_seed, h_key)
        else:
            return apply_uhash(data_set, hash_mode, k, inc, rand_seed, -1)

            
# the following are for unit-testing this class
def find_kmers(conf, chr_range, seg_range, modes, get_str=False):
    cMs = HaplotypeArray.read_cM(conf, chr_range)
    k = 2
    idx = [[] for _ in range(23)]
    for chrom in chr_range:
        for seg_id in seg_range:
            init_conf(conf, conf['n'], n)
            print(f'chrom {chrom} seg {seg_id}', end=' ')
            conf_party = conf.copy()
            conf_party['test_dir'] += f'party{0}/'
            builder = KmerBuilder(conf_party)
            haps = HaplotypeArray(conf_party, chrom, seg_id, cMs[chrom])
            maf = haps.read_mafs()
            cM_begin = seg_id * conf['step_cM'] 
            filt = builder.sample_indices(maf, haps.get_cM(), cM_begin)
            idx[chrom].append(filt)

    sums = [[] for seg_id in seg_range for chrom in range(23)]
    sums_hash = [[] for seg_id in seg_range for chrom in range(23)]
    for chrom in chr_range:
        for seg_id in seg_range:
            GG = []
            LL = []
            for i in range(k):
                conf_party = conf.copy()
                conf_party['test_dir'] += f'party{i}/'
                builder = KmerBuilder(conf_party)
                haps = HaplotypeArray(conf_party, chrom, seg_id, cMs[chrom])
                kmers_hash = builder.get_kmers(haps, indices=idx[chrom][seg_id])
                if get_str:
                    kmers = builder.get_kmers(haps, indices=idx[chrom][seg_id], debug=True)
                    GG.append(kmers)
                LL.append(kmers_hash)
            sums[chrom].append(GG)
            sums_hash[chrom].append(LL)
    return sums_hash, sums 

if __name__ == '__main__':
    # read in some haplotypes
    from user_config import conf
    n = 100
    seg_range = range(0, 2)
    chr_range = range(1, 2)

    conf['test_dir'] = 'data/2party_{}/'
    conf['mode'] = 4
    init_conf(conf, n)
    L0, L1 = find_kmers(conf, chr_range, seg_range, [1, 3, 4], get_str=False)