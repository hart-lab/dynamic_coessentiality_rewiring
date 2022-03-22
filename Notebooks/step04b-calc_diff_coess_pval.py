#!/usr/bin/python

import pandas as pd
import numpy as np
import scipy.stats as stats
import sys as sys
import random as random

def ut_as_list( dframe, diag=1, cols=['Row','Column','Value'] ):
  """
  for a symmetric dataframe, where cols=rows, get the upper triangle as a list of row/column pairs
  diag = 1 (default): ignore diagonal
  diag = 0: include diagonal
  """
  #if (dframe.index.name == dframe.columns.name):
  #dframe.index.name = cols[0]
  #dframe.columns.name = cols[1]
  #		dframe.index.name = dframe.index.name + '.1'
  #		dframe.index.name = dframe.index.name + '.2'
  d = dframe.where( np.triu( np.ones( dframe.shape ), k=diag).astype(np.bool))
  d.index.rename(cols[0], inplace=True)
  d.columns.rename(cols[1], inplace=True)
  d = d.stack().reset_index()
  d.columns=cols
  return d

def qnorm_dataframe( data ):
    """
    quantile normalize a dataframe with numeric values only! does not account for ties.
    """
    rank_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
    qnormed_data    = data.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return qnormed_data


features = pd.read_csv('./Data/features-bool-deDuped-2918feats-808cells.txt', sep='\t', index_col=0)
features.fillna(0, inplace=True)

bf = pd.read_table('./Data/table_Avana2020Q4_CRISPRcleanR_corrected_all', sep='\t', index_col=0)
genes_with_nans_idx = np.where( bf.isna().sum(1) > 0)[0]
bf.drop( bf.index.values[genes_with_nans_idx], axis=0, inplace=True)

qbf = qnorm_dataframe(bf)

cc_all = pd.DataFrame( index=qbf.index.values, columns=qbf.index.values, data=np.corrcoef( qbf.values))


strong_features_list = []
fin = open('./Data/unique_strong_features.txt', 'r')
for line in fin:
    this_feature = line.rstrip()
    if ( this_feature[0] != '#' ):
        strong_features_list.append( line.rstrip() )
fin.close()

numGenes, totNumCells = qbf.shape
cells = qbf.columns.values


for my_feature in strong_features_list:
    #
    # calculate baseline dPCC for my_feature
    #
    #my_feature = 'BRAF_GOF'

    cc_leave_out_feature = pd.DataFrame( index=qbf.index.values, 
            columns=qbf.index.values, 
            data=np.corrcoef( qbf.drop( qbf.columns.values[ features[ my_feature] ], axis=1).values))

    delta_cc = cc_all - cc_leave_out_feature

    #
    # shuffle cell lines, randomly remove n=features, calc cc, calc delta cc, compare to delta_cc
    #

    numToDrop = features[ my_feature ].sum()
    totSampledCells = totNumCells - numToDrop

    deltacc_pval = pd.DataFrame( index=qbf.index.values, 
            columns=qbf.index.values, 
            data=np.zeros( [numGenes, numGenes]) )

    for i in range(1000):
        sampled_cells = random.sample( set(cells), totSampledCells )
        sampled_cc = pd.DataFrame( index=qbf.index.values, 
                columns=qbf.index.values, 
                data=np.corrcoef( qbf[sampled_cells].values) )
        sample_deltacc = cc_all - sampled_cc
        deltacc_pval = deltacc_pval + (sample_deltacc > delta_cc) # boolean converts to integer

    cc_all_list = ut_as_list(cc_all, cols=['Gene1','Gene2','PCC_all'] )
    
    cc_leaveOneOut_list = ut_as_list( cc_leave_out_feature, cols=['Gene1','Gene2','PCC_L1out'])
    cc_all_list['PCC_L1out'] = cc_leaveOneOut_list.PCC_L1out
    del(cc_leaveOneOut_list)

    deltacc_list = ut_as_list(delta_cc, cols=['Gene1','Gene2','dPCC'] )
    cc_all_list['dPCC'] = deltacc_list.dPCC
    del(deltacc_list)

    deltacc_pval_list = ut_as_list( deltacc_pval, cols=['Gene1','Gene2','nRandom_dPCC_GT_Obs'])
    cc_all_list['nRandom_dPCC_GT_Obs'] = deltacc_pval_list.nRandom_dPCC_GT_Obs
    del(deltacc_pval_list)

    filename = 'dPCC-' + my_feature + '.txt'
    cc_all_list.to_csv('./Data/' + filename, sep='\t', index=False, float_format='%4.3f')
    del(cc_all_list)

