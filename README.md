# scBiMapping
Fast and Accurate Non-linear Dimensionality Reduction and Cell Annotation for Large and High-dimensional Single-Cell Datasets

# Install 
pip install scBiMapping

*note1: the current version (v0.0.8) on PyPI has been compiled in python 3.11, and python 3.11 is necessary for using this package.*

*note 2: if you are a BGIer, you can directly use the public image (named scBiMapping) on the cloud platform.*

# How to use 
There are two major functions in scBiMapping, **scBiMapping_DR** and **scBiMapping_annotation**, corresponding to the following two tasks.
 
## Task 1: Dimension reduction

**scBiMapping_DR(adata,n_embedding):**

* Input: 
  * adata: [anndata format](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) (cell-feature sparse matrix is stored in adata.X);   

  * n_embedding: an integer, denoting embedding dimensions (default to 30; slight adjustment may lead to better performance in practice);

* Output:
  * Embedded matrix is stored in adata.obsm['U'], where each row denotes the embedded vector of one cell;

## Task 2: reference-based cell type annotation
**scBiMapping_annotation(adata_ref,adata_query,n_embedding, K, K_majority, CellType_Key_for_ref)**  

* Key Inputs: 
  * **adata_ref**: referenc dataset (anndata format);

  * **adata_query**: query dataset (anndata format); **Note: the feature set of reference and query datasets should be the same, by using the following setttings for instance**
 
    * intersection_feature = list(set(adata_ref.var_names) & set(adata_query.var_names))
    * adata_ref = adata_ref[:,intersection_feature]
    * adata_query = adata_query[:,intersection_feature]

  * **n_embedding**: an integer, denoting the number of embeddings (default to 30; slight adjustment may lead to better performance in practice);  

  * **K**: an integer, denoting how many features are used as the new vector representation of each cell in the embedding (default to 30; adjustment may be needed in practice); 

  * **K_majority**: an integer, denoting how many reference cells are used for majority voting (default to 10; adjustment may be needed in practice);
 
  * **CellType_Key_for_ref**: key in adata_ref.obs that stores the cell type labels of the reference cells (**IMPORTANT!!!**);
 
 * Output:
   * the predicted cell types for all query cells are stored in adata_query.obs['cell_type_predicted']

## Tutorials for tasks 1 and 2

We provide several demos to further demonstrate how to conduct dimension reduction and reference-based cell type annotation using scBiMapping; Details can be found in https://github.com/scBGI/scBiMapping/tree/main/Turtorials.

## Scripts to reproduce primary experimental results

See the corresponding files in this github. See also reproducible program in codeOcean: https://codeocean.com/capsule/3904732/tree.

