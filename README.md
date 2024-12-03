# scBiMapping
Fast and Accurate Non-linear Dimensionality Reduction and Cell Annotation for Large Single-Cell Datasets

# install 
pip install scBiMapping

*note: the source code has not yet been uploaded (we will upload it once our paper is published); the currently uploaded codes have been compiled in python 3.11, and thus to run the package in python 3.11 is necessary)*

# how to use 
*there are two major functions in scBiMapping: scBiMapping_DR and scBiMapping_annotation*
 
## Task 1: Dimension reduction

**scBiMapping.scBiMapping_DR(adata,n_embedding,normalization = True):**

*adata: anndata;*    

*n_embedding: number of embeddings (default to 30);*  

*whether to normalize each embedded vector to be of norm one (default to True);*


## Task 2: cell type annotation
**scBiMapping.scBiMapping_annotation(adata_ref,adata_query,n_embedding = 30,normalization = True, K = 5, K_majority = 5, knnMethod = 'HNSW',reduction_method_on_cells_only = 'BiMapping',metric = 'euclidean',n_embedding_2nd = None, CellType_Key_for_ref = 'cell_annotation')**  

*adata_ref: referenc dataset (anndata format);* 

*adata_query: query dataset (anndata format);*  

*n_embedding: number of embedding (default to 30);*  

*normalization: whether to normalize each embedded vector to be of norm one (default to True);*  

*K: how many features are used as the new vector representation of each cell in the embedding (default to 5);*  

*K_majority: how many reference cells are used for majority voting (default to 5);*

*knnMethod: fast k-nearest neighbor searching method: 'HNSW' (default) or 'NNDescent';*

*reduction_method_on_cells_only: dimension reduction on the new representation in the embedded space: 'BiMapping' (default) or 'None';*

*metric: metric in the embedded space: 'euclidean' (default),'cosine', or, 'ip';*

*n_embedding_2nd: numbe of embeddings in the 2nd time dimension reduction: None (n_embedding will be used) or a value specfied by users;*

*CellType_Key_for_ref*: key in adata_ref.obs that stores the cell type labels of the reference cells* (**IMPORTANT!!!**);

# Tutorials




