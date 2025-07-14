import pandas as pd
import time
import os
import sys

file_bed, out_dir = sys.argv[1:]

sample = os.path.basename(file_bed).split('.')[0]
file_bps = f'{out_dir}/{sample}.bps.corrected.csv'
src_path=os.path.dirname(sys.argv[0])

df_tpm  = pd.read_csv(f'{src_path}/../Ref/TPM_atlas.csv')
df_TSSs = pd.read_table(f'{src_path}/../Ref/ucsc.refgene.hg38.TSS.UD1000.bed',names=['chr','start','end','TSS'])
df_covs_raw = pd.read_csv(file_bps)
df_covs_sum = df_covs_raw[['bps', 'bps_WGS_corrected']].sum()


TSSs = df_tpm[df_tpm['Specific high expression tissue']=='Whole Blood']['TSS ID']
df_TSSs_target = df_TSSs[df_TSSs['TSS'].isin(TSSs)].set_index('TSS')
df_covs = df_covs_raw[df_covs_raw['TSS'].isin(TSSs)]


index1 = pd.MultiIndex(levels=[df_TSSs_target.index, ['cov_raw','cov_WGS_corrected']], 
                     codes=[list(range(0,len(df_TSSs_target.index)))*2, [0, 1]*len(df_TSSs_target.index)],
                     names=['TSS', 'pos'])
index2 = pd.MultiIndex(levels=[df_TSSs_target.index, ['cnt_raw','cnt_WGS_corrected']], 
                     codes=[list(range(0,len(df_TSSs_target.index)))*2, [0, 1]*len(df_TSSs_target.index)],
                     names=['TSS', 'pos'])
df_cov = pd.DataFrame(index=index1)
df_cnt = pd.DataFrame(index=index2)

UD_1000 = list(range(-1000,1001))
UD_1000_0_list = [0]*2001
start_time = time.time()
for TSS in TSSs:
    TSS_chr, TSS_start, TSS_end = df_TSSs_target.loc[TSS][['chr','start','end']].values
    sub_bed = df_covs[(df_covs['TSS']==TSS)&(df_covs['start'] < TSS_end)&(df_covs['end'] > TSS_start)]
    if len(sub_bed)!=0:
        sub_bed['central'] = sub_bed.apply(lambda x: int((x['end']+x['start'])/2), axis=1)
        sub_bed['cov_raw'] = 1
        sub_bed['cov_WGS_corrected'] = sub_bed.apply(lambda x: x['cov_raw']/x['WGS_GC_bias'] if x['WGS_GC_bias'] > 0.05 else 0, axis=1)
        cov_raw_list=[]
        cov_wgs_list=[]
        cnt_raw_list=[]
        cnt_wgs_list=[]
        for loc in range(TSS_start, TSS_end):
            covs = sub_bed[(sub_bed['start'] <= loc )&(sub_bed['end'] > loc)][['cov_raw','cov_WGS_corrected']].sum()
            cov_raw_list.append(covs['cov_raw'])
            cov_wgs_list.append(covs['cov_WGS_corrected'])
            cnts = sub_bed[(sub_bed['central'] == loc)][['cov_raw','cov_WGS_corrected']].sum()
            cnt_raw_list.append(cnts['cov_raw'])
            cnt_wgs_list.append(cnts['cov_WGS_corrected'])
        print(TSS,time.time()-start_time)
        df_cov.loc[(TSS,'cov_raw'),UD_1000] = cov_raw_list
        df_cov.loc[(TSS,'cov_WGS_corrected'),UD_1000] = cov_wgs_list
        df_cnt.loc[(TSS,'cnt_raw'),UD_1000] = cnt_raw_list
        df_cnt.loc[(TSS,'cnt_WGS_corrected'),UD_1000] = cnt_wgs_list
    else:
        df_cov.loc[(TSS,'cov_raw'),UD_1000] = UD_1000_0_list
        df_cov.loc[(TSS,'cov_WGS_corrected'),UD_1000] = UD_1000_0_list
        df_cnt.loc[(TSS,'cnt_raw'),UD_1000] = UD_1000_0_list
        df_cnt.loc[(TSS,'cnt_WGS_corrected'),UD_1000] = UD_1000_0_list
    print(TSS,time.time()-start_time)


df_cov_mean = df_cov.groupby('pos').mean().T
df_cov_mean['cov_raw'] = df_cov_mean['cov_raw']*1e6*2001/df_covs_sum['bps']
df_cov_mean['cov_WGS_corrected'] = df_cov_mean['cov_WGS_corrected']*1e6*2001/df_covs_sum['bps_WGS_corrected']
df_cov_mean.to_csv(f'{out_dir}/{sample}.cov.mean.csv')


df_cnt_mean = df_cnt.groupby('pos').mean().T
df_cnt_mean['cnt_raw'] = df_cnt_mean['cnt_raw']*1e6*2001/df_covs_raw.apply(lambda x: int(2*x['bps']>x['end']-x['start']),axis=1).sum()
df_cnt_mean['cnt_WGS_corrected'] = df_cnt_mean['cnt_WGS_corrected']*1e6*2001/df_covs_raw.apply(lambda x: 1/x['WGS_GC_bias'] if x['WGS_GC_bias'] > 0.05 else 0,axis=1).sum()
df_cnt_mean.to_csv(f'{out_dir}/{sample}.cnt.mean.csv')
