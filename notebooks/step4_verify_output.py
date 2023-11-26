from my_util import *
from user_config import conf

exec(param_from_bash())
args = eval(ARG_STR)

rid, rep = read_cmd_args(conf, args=args)
related = read_king_file_raw(f"{conf['test_dir']}/ground_truths/KING.dat")

max_kings, degrees = find_max_king(conf)

tests = ['demo']
names = ['70%']
data, caughts, missed, extras = [], [], [], []
for j in range(len(tests)):
    numBlock = 1
    kings = max_kings
    degs = degrees
    IDs_parties, scores_parties = get_data(conf, tests[j], numBlock)
    d, c, m, e= get_classes(IDs_parties, scores_parties, kings, degs)
    d['test'] = names[j]
    data.append(d)
    caughts.append(c)
    missed.append(m)
    extras.append(e)


