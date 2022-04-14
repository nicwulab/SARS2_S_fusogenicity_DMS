#!/usr/bin/python
import sys
import pandas as pd
import numpy as np

def count_to_freq(df, colname):
    df[colname+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return (df)

def exp_score_calculate(df, rep, freq_cutoff):
    print (rep)
    exp_weight = df['bin0_'+rep+'_freq']*0.25 + df['bin1_'+rep+'_freq']*0.5 + \
                 df['bin2_'+rep+'_freq']*0.75  + df['bin3_'+rep+'_freq']*1
    exp_norm_factor = df['bin0_'+rep+'_freq'] + df['bin1_'+rep+'_freq'] + \
                      df['bin2_'+rep+'_freq'] + df['bin3_'+rep+'_freq']
    df['exp_weight_'+rep] = exp_weight/exp_norm_factor
    df_high_freq = df[df['input_freq'] >= freq_cutoff]
    w_summary = df_high_freq.groupby('mut_class')['exp_weight_'+rep].mean()
    w_summary = w_summary.reset_index()
    print (w_summary)
    w_silent   = (float(w_summary.loc[w_summary['mut_class']=='silent']['exp_weight_'+rep]))
    w_nonsense = (float(w_summary.loc[w_summary['mut_class']=='nonsense']['exp_weight_'+rep]))
    df['exp_score_'+rep] = (df['exp_weight_'+rep]-w_nonsense)/(w_silent-w_nonsense)
    print ('w_nonsense', w_nonsense)
    print ('w_silent', w_silent)
    return (df)

def fusion_score_calculate(df, rep, freq_cutoff):
    df['fus_ratio_'+rep] = np.log10(df['fus_pos_'+rep+'_freq']/df['fus_neg_'+rep+'_freq'])
    df_high_freq = df[df['input_freq'] >= freq_cutoff]
    fus_summary = df_high_freq.groupby('mut_class')['fus_ratio_'+rep].mean()
    fus_summary = fus_summary.reset_index()
    fus_silent   = (float(fus_summary.loc[fus_summary['mut_class']=='silent']['fus_ratio_'+rep]))
    fus_nonsense = (float(fus_summary.loc[fus_summary['mut_class']=='nonsense']['fus_ratio_'+rep]))
    df['fus_score_'+rep] = (df['fus_ratio_'+rep]-fus_nonsense)/(fus_silent-fus_nonsense)
    return (df)

def read_count_data(count_file, freq_cutoff):
    print ('reading: %s' % count_file)
    df = pd.read_csv(count_file, sep='\t')
    colnames = [colname for colname in df]
    for colname in colnames:
        if 'mut' not in colname and 'frag' not in colname:
            df = count_to_freq(df, colname)
    df = exp_score_calculate(df, 'rep1', freq_cutoff)
    df = exp_score_calculate(df, 'rep2', freq_cutoff)
    df = exp_score_calculate(df, 'rep3', freq_cutoff)
    df = fusion_score_calculate(df, 'rep1', freq_cutoff)
    df = fusion_score_calculate(df, 'rep2', freq_cutoff)
    df['exp_score'] = (df['exp_score_rep1'] + df['exp_score_rep2'] + df['exp_score_rep3'])/3
    df['fus_score'] = (df['fus_score_rep1'] + df['fus_score_rep2'])/2
    return (df)

def main():
    #freq_cutoff = 0.0001
    freq_cutoff = 0
    outfile_1 = "result/S2_HR1_DMS_scores.tsv"
    outfile_2 = "result/S2_HR1_DMS_scores_by_resi.tsv"
    df = read_count_data("result/S2_HR1_DMS_count_aa.tsv", freq_cutoff)
    df = df[df['input_freq'] > freq_cutoff]
    print ('writing: %s' % outfile_1)
    df.to_csv(outfile_1, sep="\t", index=False)
    df['resi'] = df['mut'].str[0:-1]
    all_resi   = df[df['mut_class'] != 'WT'].groupby('resi').size().reset_index()
    df_by_resi = df
    df_by_resi = df_by_resi[df_by_resi['input_freq'] >= freq_cutoff]
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'WT']
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'silent']
    df_by_resi = df_by_resi[df_by_resi['mut_class'] != 'nonsense']
    df_by_resi_mean  = df_by_resi.groupby('resi')['exp_score'].mean().reset_index(name='mean_exp_score')
    df_by_resi_count = df_by_resi.groupby('resi').size().reset_index(name='count')
    df_by_resi = pd.merge(df_by_resi_mean, df_by_resi_count, on='resi', how='outer')
    df_by_resi = pd.merge(all_resi, df_by_resi, on='resi', how='outer')
    df_by_resi = df_by_resi.sort_values(by='resi', key=lambda x:x.str[1::].astype(int))
    df_by_resi['pos'] = df_by_resi['resi'].str[1::].astype(int)
    df_by_resi = df_by_resi[['resi','pos','count','mean_exp_score']]
    print ('writing: %s' % outfile_2)
    df_by_resi.to_csv(outfile_2, sep="\t", index=False)

if __name__ == "__main__":
    main()
