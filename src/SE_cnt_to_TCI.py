### TPM CNT to COVs and TCI

import pandas as pd
import os,sys

TPM, CNT, COVs, TCI = sys.argv[1:]

dat_in = pd.read_table(CNT, sep = '\t', names = ['TSS'], usecols = [3])
ID = os.path.basename(CNT).split('.')[0]
dat_in[ID] = 1
dat_in = dat_in.groupby('TSS').sum()
covs = pd.DataFrame(dat_in[ID]*10**6 / dat_in[ID].sum())

covs.to_csv(COVs, sep = '\t', index_label = 'TSS')

dep = pd.read_table(COVs, sep = '\t', index_col = [0]).transpose()
tpm = pd.read_table(TPM, sep = ',')
tci = dep[[]]

for tissue_, sub_dat in tpm.groupby('Specific high expression tissue'):
    tissue = tissue_.split(' -')[0]
    TSSs = list(tpm[tpm['Specific high expression tissue'] == tissue_]['TSS ID'])
    if tissue == 'Whole Blood':
        tci['Blood cell'] = 50 - (dep[list(set(TSSs).intersection(set(dep.columns)))].sum(axis=1) / len(TSSs))
    else:
        tci[tissue] = 50 - (dep[list(set(TSSs).intersection(set(dep.columns)))].sum(axis=1) / len(TSSs))

tci.to_csv(TCI, index_label = 'sample')
