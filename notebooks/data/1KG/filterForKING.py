import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
data = pd.read_csv(f'pcavar.eigenvec.var', sep='\t')

wts = np.max(np.asarray([abs(data['PC1']), abs(data['PC2']), abs(data['PC2'])]), axis=0)
# wts = np.max(np.asarray([abs(okg_wts['PC1_loading']), abs(okg_wts['PC2_loading']), abs(okg_wts['PC2_loading'])]), axis=0)
plt.hist(wts, bins=100)
thresh = 75
cutoff = np.round(np.percentile(wts, thresh), 3)
percentile_95 = np.percentile(wts, 99)
plt.xlim(0,percentile_95)
plt.axvline(x=cutoff, color='red')
plt.show()
snps_to_keep = data.loc[wts < cutoff, 'ID']
print(len(snps_to_keep), ' vs ', len(data))
np.savetxt(f'snpsForKING.txt', snps_to_keep, fmt='%s')