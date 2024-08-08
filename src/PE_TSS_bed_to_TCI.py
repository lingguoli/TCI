### TPM BED to COVs and TCI
import pandas as pd
import os, sys

TPM, BED, COVs, TCI = sys.argv[1:]
ID = os.path.basename(BED).split('.')[0]
TSS_dic = {}
with open(BED) as f:
    for line in f:
        TSS, cov = line.strip().split('\t')[3], line.strip().split('\t')[-1]
        #_,_,_,TSS,_,_,_,_,_,_,cov = line.strip().split('\t')
        try:
            TSS_dic[TSS] = TSS_dic[TSS] + int(cov)
        except:
            try:
                TSS_dic[TSS] = int(cov)
            except:
                continue
dat_in = pd.DataFrame(TSS_dic, index = [ID]).transpose()
covs = pd.DataFrame(dat_in[ID] * 10**6 / dat_in[ID].sum())
covs.to_csv(COVs, sep='\t', index_label = 'TSS')

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
