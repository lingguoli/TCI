import os
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.stats import norm
import math
from math import sqrt
import sklearn
from sklearn.metrics import roc_curve
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests


# Calculate the 95% confidence interval for the Pearson correlation coefficient
def bootstrap_ci(data_x, data_y, n_bootstrap=1000):
    bootstrapped_corrs = []
    np.random.seed(42)
    for _ in range(n_bootstrap):
        indices = np.random.randint(0, len(data_x), len(data_x))
        sampled_x = data_x[indices]
        sampled_y = data_y[indices]
        boot_corr, _ = stats.pearsonr(sampled_x, sampled_y)
        bootstrapped_corrs.append(boot_corr)
    lower = np.percentile(bootstrapped_corrs, 2.5)
    upper = np.percentile(bootstrapped_corrs, 97.5)
    return lower, upper

# Calculate the 95% confidence interval for the Spearman correlation coefficient
def bootstrap_ci2(data_x, data_y, n_bootstrap=1000):
    bootstrapped_corrs = []
    np.random.seed(42)
    for _ in range(n_bootstrap):
        indices = np.random.randint(0, len(data_x), len(data_x))
        sampled_x = data_x[indices]
        sampled_y = data_y[indices]
        boot_corr, _ = stats.spearmanr(sampled_x, sampled_y)
        bootstrapped_corrs.append(boot_corr)
    lower = np.percentile(bootstrapped_corrs, 2.5)
    upper = np.percentile(bootstrapped_corrs, 97.5)
    return lower, upper

# Calculate the 95% confidence interval for AUC
def roc_auc_ci(y_true, y_score, positive=1):
    AUC = roc_auc_score(y_true, y_score)
    N1 = sum(y_true == positive)  
    N2 = sum(y_true != positive)  
    Q1 = AUC / (2 - AUC)
    Q2 = 2 * AUC**2 / (1 + AUC)
    SE_AUC = sqrt((AUC * (1 - AUC) + (N1 - 1) * (Q1 - AUC**2) + (N2 - 1) * (Q2 - AUC**2)) / (N1 * N2))
    z = norm.ppf(1 - (1 - 0.95)/2)
    lower = AUC - z * SE_AUC
    upper = AUC + z * SE_AUC
    return (AUC, max(0, lower), min(1, upper))


# p value correction
def p_adjust(p_list):
    q_list = multipletests(p_list, alpha=0.05, method='bonferroni')[1]
	return q_list


# Single ROC plot
def AUC_plot(data, group, pos_label, pre_value, title):
    #AUC
    fpr,tpr,thresholds=roc_curve((data[group]==pos_label).astype('int'), data[pre_value])  #fpr,tpr,thresholds
    auc, lower_ci, upper_ci=roc_auc_ci((data[group]==pos_label).astype('int'), data[pre_value])
    auc, lower_ci, upper_ci = "{:.2f}".format(auc), "{:.2f}".format(lower_ci), "{:.2f}".format(upper_ci)
    AUC=f'{auc}[{lower_ci},{upper_ci}]'
    #plot
    fig1 = plt.figure()
    ax3 = plt.axes([0,0,0.5,0.65])
    ax3.plot(fpr*100,tpr*100,color='black')
    ax3.text(15,0,f'AUC:{AUC}',size=16,color='black',weight='bold')
    ax3.set_title(title,size=20)
    ax3.set_ylabel('Sensitivity(%)', size=18)
    ax3.set_xlabel('1-Specificity(%)', size=18)
	return fig1, ax1

# Multiple ROC plots
def AUC_plots(data_list, group_list, pos_label_list, pre_value_list, legend_list, title):
    fig1 = plt.figure()
    ax3 = plt.axes([0,0,0.5,0.65])
    i = 0
    platte = ['#d62728','#1f77b4','#2ca02c','#ff7f0f']
    line_style = ['--','-']
    for data, group, pos_label, pre_value in zip(data_list, group_list, pos_label_list, pre_value_list):
        fpr,tpr,thresholds=roc_curve((data[group]==pos_label).astype('int'), data[pre_value]) 
        auc, lower_ci, upper_ci=roc_auc_ci((data[group]==pos_label).astype('int'), data[pre_value])
        auc, lower_ci, upper_ci = "{:.2f}".format(auc), "{:.2f}".format(lower_ci), "{:.2f}".format(upper_ci)
        AUC=f'{auc}[{lower_ci},{upper_ci}]'
        ax3.plot(fpr*100,tpr*100,color=platte[i], linestyle = line_style[i])
        ax3.text(10,0+8*i,f'{legend_list[i]}:{AUC}',size=14,color=platte[i],weight='bold')
        i = i + 1
    ax3.set_title(title,size=20)
    ax3.text(10,0+8*i,f'AUC:',size=18,color='black',weight='bold')
    ax3.set_ylabel('Sensitivity(%)', size=18)
    ax3.set_xlabel('1-Specificity(%)', size=18)
    return fig1, ax1
	
# tSNE plot
def tSNE_plot(df_TPM, tissue_group):
	perplexity = min(30, df_TPM.shape[0] - 1)
	tsne = TSNE(n_components=2, perplexity=perplexity, n_iter=300)
	tsne_results = tsne.fit_transform(df_TPM)
	df_tsne_results = pd.DataFrame(tsne_results)
	df_tsne_results.columns = ['t-SNE Component 1','t-SNE Component 2']
	df_tsne_results['Tissue'] = tissue_group['group'].values
	fig1=plt.figure()
	ax1=plt.axes([0,0,2,1])
	sns.scatterplot(data = df_tsne_results, x = 't-SNE Component 1', y = 't-SNE Component 2',hue = 'Tissue', style = 'Tissue', palette = 'tab10_r')
	ax1.set_title(f'',size=25)
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles=handles, labels=labels, bbox_to_anchor=(0, 1.04), loc=2, frameon=True, borderaxespad=1.2,title='Tissue')
	ax1.yaxis.tick_right()
	#ax1.set_ylim(top=24)
	ax1.set_xlim(left=-22)
	ax1.yaxis.set_label_position("right")
	ax1.set_ylabel('t-SNE Component 2',size=20)
	ax1.set_xlabel('t-SNE Component 1',size=20)
	return fig1, ax1
