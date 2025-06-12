import os
import sys
import time
import warnings
import random
import gc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from scipy.sparse import *
import sklearn
from scBiMapping import *

color = {'L2/3':[205,205,0],'L3/4':[192,255,62],'L4/5':[72,118,255],'L4/5/6':[255,187,255],'L2':[255,130,71],'L4':[152,245,255],'L5/6':[171,130,255],'L6':[138,43,226],'L2/3/4':[255,255,0],'L3/4/5':[127,255,212],'PVALB':[79,148,205],'SST':[82,139,139],'RELN':[205,133,0],'LAMP5':[255,20,147],'VIP_RELN':[205,155,155],'VIP':[219,112,147],'PV_CHC':[176,226,255],'ASC':[255,106,106],'OLG':[205,205,180],'OPC':[255,182,193],'MG':[255,218,185],'EC':[255,211,155],'VLMC':[238,216,174]}
for i in color:
    color[i] = matplotlib.colors.to_hex(np.array(color[i])/255)

dir_in_scRNA = 'scRNA'        
dir_out = 'Final'          

# load StereoSeq data
dir_in_stereoSeq = 'StereoSeq'
adata_query = anndata.read(dir_in_stereoSeq+'/h5ad/monkey1.h5ad') # monkey1, monkey2, monkey3
 
# load scRNA data
cell_annotation = pd.read_csv(dir_in_scRNA + '/snRNA.metadata.2monkeys.csv',header = 0, index_col = 0)
adata_ref = anndata.read_h5ad('monkey1and2_reference.h5ad')
print(adata_ref.X.max())
 
sc.pp.normalize_total(adata_ref, target_sum=1e4)
sc.pp.log1p(adata_ref)
sc.pp.highly_variable_genes(adata_ref, n_top_genes=6000,subset=True,flavor= 'cell_ranger') # 

# pre-processing Query 
# Note: here adata_query.X is the normalized expression calculated by "SCTnormalize"

# get intersection genes
intersection_feature = list(set(adata_ref.var_names) & set(adata_query.var_names))
adata_ref = adata_ref[:,intersection_feature]
adata_query = adata_query[:,intersection_feature]
print(adata_ref)
print(adata_query)

# get class types for reference data
cell_names = adata_ref.obs_names.tolist() 
new_cell_names = [name.rsplit('-', 1)[0] for name in cell_names] 
adata_ref.obs_names = pd.Index(new_cell_names)
cell_annotation = cell_annotation.loc[adata_ref.obs_names]
adata_ref.obs['subclass'] = cell_annotation['SubClass'].values  # 23 types

# Main code: obtain nearest reference cells for each query cell
n_embedding = 30; 
K = 30; # this parameter will heavily influence the speed. Don't be too large
K_majority = 30; # parameter for marjority voting
knnMethod = 'HNSW';

adata_ref.obs['cell_annotation'] = adata_ref.obs['subclass'].values 

adata_query.obs['cell_type'] = scBiMapping_annotation(adata_ref, adata_query,n_embedding = n_embedding,K = K,K_majority = K_majority,knnMethod = knnMethod)

adata_query.obs['cell_type']=adata_query.obs['cell_type'].astype("category")
adata_query.obs['cell_type_colors'] = [color[i] for i in adata_query.obs['cell_type'].values.tolist()]
adata_query.obs.to_csv(dir_out + '/' + 'all_2.csv')
 
