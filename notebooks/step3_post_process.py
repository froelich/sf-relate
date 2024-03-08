import numpy as np
import tomlkit
import pandas as pd
from my_util import *
import os

if __name__ == "__main__":
    # read environment variables in python
    parser = argparse.ArgumentParser(description='Post process outputs')
    parser.add_argument('-PARTY', type=int, help='party id', required=True)
    parser.add_argument('-FOLDER', type=str, help='path to the configuration folder', required=True)
    args = parser.parse_args()
    party_id = args.PARTY
    FOLDER = args.FOLDER
    args_local = tomlkit.load(open(f"{FOLDER}/configLocal.Party{party_id}.toml"))
    args_global = tomlkit.load(open(f"{FOLDER}/configGlobal.toml"))
    # find out which line in pvar_file start with #CHROM
    psam, snps, pvar = read_pgen_metadata(args_local)
    n = len(psam)
    # read in the results
    col_name = 'Result'
    reveal = args_global['reveal']
    if reveal == 0:
        name = ''
        nOutputs = 1
        out_name = 'boolean'
    elif reveal == 1:
        name = 'degree'
        nOutputs = 4
        out_name = 'degree'
    elif reveal == 2:
        name = 'kinship'
        nOutputs = len(args_global['discretized_thresh'])
        out_name = 'bin'
    elif reveal == 3:
        name = 'kinship_block'
        out_name = 'kinship'
    else:
        raise ValueError('reveal should be 0, 1, 2, or 3')

    out_dir = FOLDER + "out/"
    if reveal == 3:
        # change all computed block's kinship to the real kinship
        all_kinship = []
        ID_pos = np.load(args_local['hash_table_dir'] + '/ID_table.npz')['ID_table']
        rem = len(ID_pos) % 8192
        if rem == 0:
            rem = 8192
        ID_pos = np.hstack((ID_pos, [n] * (8192 - rem)))
        nB = len(ID_pos) // 8192
        for bb in range(nB):
            all_kinship.append(np.loadtxt(f"{out_dir}/raw/kinship_block_{bb}_party{party_id}.txt", dtype=float))
        all_kinship = np.concatenate(all_kinship, axis=0)
        IIDs = np.hstack((psam['#IID'], [0]))
        ID_pos[ID_pos == -1] = n
        # extend ID_pos by n until it's a multiple of 8192
        ls_psam = IIDs[ID_pos]
        computed = pd.DataFrame({"#IID" : ls_psam, "kinship" : all_kinship})
        ls_psam.to_csv(f"{out_dir}/kinship_block_party{party_id}.tsv", index=False, sep='\t')
    else:
        nB = int(np.ceil(n / args_global['batch_length']))
        stat_computed = pd.DataFrame({'ID' : np.arange(n), '#IID': psam['#IID'], col_name : nOutputs})
        for j in range(nOutputs):
            for bb in range(nB):
                IDs = []
                scores = []
                with open(f"{out_dir}/raw/{name}{j}_{bb}_party{party_id}.csv") as f:
                    q = 0
                    for line in f:
                        q += 1
                        if q == 1:
                            continue
                        a, b = line.split(',')
                        ID = int(a) + bb * args_global['batch_length']
                        IDs.append(ID)
                        b = b.replace('i', 'j')
                        gg = complex(b[1:-2])
                        scores.append(gg.real)
                    update_ID = np.array(IDs)
                    update_score = np.array(scores)
                    update_ID = update_ID[np.round(update_score) == 0]
                    # update 'deg to be j if those are larger
                    filt = stat_computed.loc[update_ID, col_name] > nOutputs - 1 - j
                stat_computed.loc[update_ID[filt], col_name] = nOutputs - 1 - j
        if reveal == 0:
            stat_computed[col_name] = 1 - stat_computed[col_name]
        
        stat_computed[['#IID', col_name]].to_csv(f"{out_dir}/{out_name}_party{party_id}.csv", index=False, sep='\t')