#%%
from bgen_reader import open_bgen
from my_util import *
import numpy as np
import pandas as pd
import sys
from joblib import Parallel, delayed
from user_config import conf

conf_raw = {
    "test_dir": "data/raw/",
    "bgen_name_template": 'chr{}.bgen',
    "n": None,
    "verbose": False,
}

conf_subset = {
    "bgen_name_template": 'chr{}.bgen',
    "n": None,
    "verbose": False,
}


def find_deg_stats(layers, bbank_n):
    num_in_class = np.zeros((4))
    densities = np.zeros((4))
    for (deg, data) in layers:
        num_in_class[deg] = int(data.shape[0])
        individuals = np.unique(np.concatenate((data['ID1'], data['ID2'])))
        densities[deg] = len(individuals) / bbank_n
    return num_in_class, densities


def read_related_info(rel_file):
    dataset = pd.read_csv(rel_file, sep=" ")
    dataset = dataset[dataset['Kinship'] != -1]

    # find samples that are not related with anyone
    all_samples = pd.read_csv(f"{conf_raw['test_dir']}chr1.sample", sep=' ', skiprows=2, header=None)[0]
    all_related_samples = np.unique(np.concatenate((dataset['ID1'], dataset['ID2'])))
    all_related_samples = all_related_samples[all_related_samples > 0]
    unrelated_samples = all_samples[~np.isin(all_samples, all_related_samples)].to_numpy()
    layers = dataset.groupby(lambda idx:  relativedegree(dataset.loc[idx]['Kinship']))

    return dataset, all_samples, all_related_samples, unrelated_samples, layers

if __name__ == '__main__':
    all_samples = pd.read_csv(f"{conf_raw['test_dir']}chr1.sample", sep=' ', skiprows=2, header=None)[0]
    all_samples = np.sort(all_samples.to_numpy())
    all_samples = all_samples[all_samples >= 0]
    bbank_n = len(all_samples)
    exec(param_from_bash())
    args = eval(ARG_STR)
    rid, rep = read_cmd_args(conf, args=args)
    k = 2
    ver = 0
    test_dir = f"data/{k}party_{n}{'' if ver == 0 else 'v'+str(ver)}/"
    list_names = [test_dir + f'party{i}/list' for i in range(k)]
    # only compute KING when we want to split by placing relatives on two sides first
    # recompute_KING(test_dir, k)
    def random_split_gen(test_dir, listnames, all_samples, n, k=2):
        run_command(['mkdir', f'{test_dir}'])
        print("Randomly generating test sets")
        list_ppl = []
        np.random.shuffle(all_samples)
        for i in range(k):
            list_ppl.append(all_samples[n * i: n * (i + 1)])

        np.random.shuffle(list_ppl[0])

        for i in range(len(listnames)):
            run_command(['mkdir', f'{test_dir}/party{i}'])
            with open(listnames[i] + '.txt', "w") as file:
                for ID in list_ppl[i]:
                    ID = int(ID)
                    file.write(str(ID) + ' ' + str(ID) + '\n')

        print(list_ppl[0])
        print(list_ppl[1])
        print(f"Sampled a total of {len(np.unique(np.concatenate((list_ppl[0], list_ppl[1]))))} individuals")

    def split_using_king_file(test_dir, list_names, n, k=2):
        rel = pd.read_csv(f"data/raw/KING.dat", sep=' ')
        rel['deg'] = relativedegree(rel['Kinship'])
        related_IDs = [[] for _ in range(2)]
        num_classes = 5
        to_sample = int(n / num_classes)
        unused = set(all_samples)
        total_n = 0
        sampled = []
        for d, data in rel.groupby('deg'):
            cnt = 0
            for row in data.sample(frac=1).itertuples():
                P0 = np.random.randint(2)
                P1 = 1 - P0
                i = row.ID1
                j = row.ID2
                if i in unused and j in unused:
                    related_IDs[P0].append(i)
                    related_IDs[P1].append(j)
                    unused.remove(i)
                    unused.remove(j)
                    cnt += 1
                    total_n += 1
                if cnt == to_sample:
                    break
            sampled.append(cnt)
        print(sampled)
        unused = list(unused)
        np.random.shuffle(unused)
        j = 0
        while total_n < n:
            total_n += 1
            j += 1
            for i in range(k):
                related_IDs[i].append(unused[2 * j + k])
        run_command(['mkdir', f'{test_dir}'])
        for i in range(k):
            run_command(['mkdir', f'{test_dir}/party{i}'])
            with open(list_names[i] + '.txt', "w") as file:
                for ID in related_IDs[i]:
                    file.write(f"{ID} {ID}\n")
        

    def split_bgen(test_dir, list_names, chr_range = (21, 0, -1), with_sample_file=True):
        # print(list_names)
        for i in range(len(list_names)):
            run_command(['mkdir', test_dir + f'party{i}/bgen'])

        command_recode = []
        for chromo in chr_range:
            for i in range(len(list_names)):
                list_name = list_names[i] + '.txt'
                bgen_name = test_dir + f"party{i}/bgen/chr{chromo}"

                command_with_file = ['./qctool',
                                    '-g', f'{conf_raw["test_dir"]}chr{chromo}.bgen']
                if with_sample_file:
                    command_with_file += ['-s', f'{conf_raw["test_dir"]}chr1.sample']

                command_recode.append(command_with_file.copy() + ['-incl-samples', list_name, '-og', f'{bgen_name}.bgen'])
        
        Parallel(n_jobs=16)(delayed(run_command)(cmd, False) for cmd in command_recode)


    def recodeIDs_chr22(test_dir, list_names):
        samples = pd.read_csv(f"{conf_raw['test_dir']}chr1.sample", sep=' ', skiprows=2, header=None)
        chr22_bgen_name = conf_raw['test_dir'] + conf_raw['bgen_name_template'].format(22)
        bgen = open_bgen(chr22_bgen_name)
        chr22_sample_IDs = np.asarray(bgen.samples.copy())
        for list_name in list_names:
            # print(pd.read_csv(list_name + '.txt', sep=' ', header=None))
            sample_subset = pd.read_csv(list_name + '.txt', sep=' ', header=None)[0]
            indices = np.isin(samples[0], sample_subset)
            special_ids = chr22_sample_IDs[indices]
            with open(f"{list_name}_chr22.txt", "w") as file:
                for ID in special_ids:
                    file.write(str(ID) + ' ' + str(ID) + '\n')


    # Now we sample two lists of IDs, containing related pairs
    if n >= 10000:
        # if n >= 10000, we randomly split
        random_split_gen(test_dir, list_names, all_samples, n)
    else:
        # Otherwise, we place relatives on the two sites and then split
        split_using_king_file(test_dir, list_names, n)

    start = start_timer()
    # call qctool to sub-sample genotypes according to the list
    chr_range = range(21, 0, -1)
    split_bgen(test_dir, list_names, chr_range)
    stop_timer(start)

    # CHROMO 22 seems to have a different set of samples
    # So have to call qctool on this and handle it specially
    recodeIDs_chr22(test_dir, list_names)
    split_bgen(test_dir, [f"{name}_chr22" for name in list_names], [22], False)
    stop_timer(start)