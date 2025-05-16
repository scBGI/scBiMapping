import anndata 
import warnings  
from scipy import sparse
from scipy.sparse import *   
import numpy as np
import pandas as pd
import sys
warnings.filterwarnings('ignore')
from scBiMapping import *

#%% scRNA
dir_in = 'allen-brain-cell-atlas/'
dir_out = 'MerFish_and_stereSeq_annotation_results' 
 
scRNA_10xV2V3_version = '20230630' # 10x v2/v3 scRNA counts: '20230630'; 
Multiome_version = '20230830'; # '20230830'; 
ref_meta_version = '20230830'; #'20230630', '20230830', '20231215'
Zhuang_MerFish_count_version = '20230830'; # '20230830'
Zhuang_MerFish_meta_version = '20230830'; # '20230830', '20231215';
Zhuang_MerFish_ccf_version = '20230830'; # '20230830'
 
adata_ref = anndata.read_h5ad('scRNA_Zeng_all_genes_10Xv2_and_v3.h5ad') # # 10x v2/v3 scRNA counts: '20230630'; 

# preprocessing
import scanpy as sc
sc.pp.normalize_total(adata_ref, target_sum=1e4)
sc.pp.log1p(adata_ref)
 
import pandas as pd
ref_meta = pd.read_csv(dir_in + '/metadata/WMB-10X/' + ref_meta_version + '/views/cell_metadata_with_cluster_annotation.csv',index_col = 'cell_label',dtype = {'cell_label':'str'})
common_cells = list(set(ref_meta.index.values).intersection(set(adata_ref.obs_names.values))); 

print(len(common_cells)); print(adata_ref.shape); print(ref_meta.shape)


adata_ref = adata_ref[common_cells,:] 
ref_meta_sub = ref_meta.loc[common_cells]
adata_ref.obs = pd.concat([adata_ref.obs,ref_meta_sub], axis=1)    

#%% MerFish-1122
monkey_indexes = [1,2,3] # 1,2,3 represents three mices
xy = {1:{'x':'ccf_y','y':'ccf_z'},2:{'x':'ccf_y','y':'ccf_z'},3:{'x':'ccf_x','y':'ccf_y'}}
for monkey_index in monkey_indexes: 
    name = 'Zhuang-ABCA-'+str(monkey_index)
    print(name)
    adata_query = anndata.read_h5ad(dir_in + '/expression_matrices/'+name+ '/' + Zhuang_MerFish_count_version + '/' +name+ '-log2.h5ad') 

    MerFish_meta = pd.read_csv(dir_in + '/metadata/' + name + '/' + Zhuang_MerFish_meta_version + '/views/cell_metadata_with_cluster_annotation.csv',index_col = 'cell_label',dtype = {'cell_label':'str'}) # equvalent to  MerFish_meta = Get_MerFish_MetaData(dir_in,Zhuang_MerFish_meta_version,name)
    if not adata_query.obs_names.values.tolist() == MerFish_meta.index.tolist():
        print(MerFish_meta.shape);   print(len(set(adata_query.obs_names.values.tolist()).intersection(set(MerFish_meta.index.tolist()))))
        common_cells = list(set(adata_query.obs_names.values.tolist()).intersection(set(MerFish_meta.index.tolist())))
        adata_query = adata_query[common_cells,:] 
        MerFish_meta = MerFish_meta.loc[common_cells]
    adata_query.obs = pd.concat([adata_query.obs,MerFish_meta], axis=1)
           
    ccf = pd.read_csv(dir_in + '/metadata/'+name+'-CCF/' + Zhuang_MerFish_ccf_version + '/ccf_coordinates.csv',index_col = 'cell_label',dtype = {'cell_label':'str'})
    adata_query = adata_query[ccf.index,:] 
    adata_query.obs[['ccf_x','ccf_y','ccf_z']] = ccf[['x','y','z']] 
    
    # subset ref_data
    adata_ref = adata_ref[:,adata_query.var_names]
    print(adata_ref); print(adata_query)
     
    n_embedding = 100; K = 20; K_majority = 5；knnMethod = 'HNSW';  
    CellType_Key_for_ref = ['class','subclass'] 
    
    scBiMapping_annotation(adata_ref, adata_query,n_embedding = n_embedding,K = K,K_majority = K_majority,knnMethod = knnMethod,CellType_Key_for_ref = CellType_Key_for_ref)
    print(adata_query) # note that adata_query is a anndata; the prediction results are also stored in it.

    #%% evaluate accuracy 
    # Case 1 
    y_true = adata_query.obs['class']
    y_pred = adata_query.obs['cell_type_predicted1']
    print('case 1: class prediction...')
    accuracy,precision,recall,f1_score = evaluate_classification(y_true, y_pred)

    # Case 2  
    y_true = adata_query.obs['subclass']
    y_pred = adata_query.obs['cell_type_predicted2']
    print('case 2: subclass prediction...')
    accuracy,precision,recall,f1_score = evaluate_classification(y_true, y_pred)
  