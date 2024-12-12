# scBiMapping
Fast and Accurate Non-linear Dimensionality Reduction and Cell Annotation for Large Single-Cell Datasets

# install 
pip install scBiMapping

*note1: the source code has not yet been uploaded (we will upload it once our paper is published); the currently uploaded codes have been compiled in python 3.11, and thus to run the package in python 3.11 is necessary).*

*note 2: if you are a BGIer, you can directly use the public image (named scBiMapping) on the cloud platform.*

# how to use 
There are two major functions in scBiMapping, **scBiMapping_DR** and **scBiMapping_annotation**, corresponding to the following two tasks.
 
## Task 1: Dimension reduction

**scBiMapping_DR(adata,n_embedding = 30, normalization = True):**

* Input: 
  * adata: [anndata format](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) (cell-feature sparse matrix is stored in adata.X);   

  * n_embedding: an integer, denoting the number of embeddings (default to 30; slight adjustment may lead to better performance in practice);

  * normalization: whether to normalize each embedded vector to norm one (default to True);

* Output:
  * Embedded matrix is stored in adata.obsm['U'], where each row denotes the embedded vector of one cell;

## Task 2: reference-based cell type annotation
**scBiMapping_annotation(adata_ref,adata_query,n_embedding = 30,normalization = True, K = 5, K_majority = 5, knnMethod = 'HNSW',reduction_method_on_cells_only = 'BiMapping',metric = 'euclidean',n_embedding_2nd = None, CellType_Key_for_ref = 'cell_annotation')**  

* Input: 
  * adata_ref: referenc dataset (anndata format);

  * adata_query: query dataset (anndata format); **Note: the feature set of reference and query datasets should be the same, by using the following setttings for instance**
 
    * intersection_feature = list(set(adata_ref.var_names) & set(adata_query.var_names))
    * adata_ref = adata_ref[:,intersection_feature]
    * adata_query = adata_query[:,intersection_feature]

  * n_embedding: an integer, denoting the number of embeddings (default to 30; slight adjustment may lead to better performance in practice);  

  * normalization: whether to normalize each embedded vector to norm one (default to True);  

  * K: an integer, denoting how many features are used as the new vector representation of each cell in the embedding (default to 5; suggested to be no larger than 100); 

  * K_majority: an integer, denoting how many reference cells are used for majority voting (default to 5; suggested to be no larger than 50);

  * knnMethod: fast k-nearest neighbor searching method: 'HNSW' (default) or 'NNDescent' (recommended as well);

  * reduction_method_on_cells_only: dimension reduction on the new representation in the embedded space: 'BiMapping' (default) or 'None';

  * metric: metric in the embedded space: 'euclidean' (default),'cosine', or, 'ip';

  * n_embedding_2nd: numbe of embeddings in the 2nd time dimension reduction: None (n_embedding will be used) or a value specfied by users;

  * CellType_Key_for_ref: key in adata_ref.obs that stores the cell type labels of the reference cells (**IMPORTANT!!!**);

 * Output:
   * the predicted cell types for all query cells are stored in adata_query.obs['cell_type_predicted']

## Tutorials for tasks 1 and 2

We provide 8 tutorials to further demonstrate how to conduct dimension reduction and reference-based cell type annotation using scBiMapping; see details at https://cloud.stomics.tech/library/#/tool/detail/workspace_notebook/NB0120241204Yng3Pn/--?zone=sz




