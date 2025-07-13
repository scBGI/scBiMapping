"""
scBiMapping: Accurate, fast, and robust dimension reduction of single-cell data with scBiMapping
License: GPL-3.0
Tuturial: https://github.com/scBGI/scBiMapping/tree/main/Turtorials

A comprehensive toolkit for single-cell RNA-seq data analysis featuring:
- Dimensionality reduction via spectral embedding
- Reference-based cell type annotation 

an open source under the MIT license.

"""

import numpy as np 
from scipy.sparse import diags, csr_matrix, coo_matrix, issparse
from random import randint  
from sklearn.preprocessing import normalize 
from sklearn.neighbors import NearestNeighbors
from pynndescent import NNDescent
from scipy import sparse
import anndata as ad 
import pandas as pd  
import warnings,random,hnswlib
from typing import Union, List, Optional
from numpy import exp, sqrt, tile, arange
rd = np.random.RandomState(888) 
 

###########################################################################
#                    dimension reduction (scBiMapping_DR)                 #
###########################################################################

def scBiMapping_DR(
    adata,
    n_embedding: int = 30
) -> ad.AnnData:
    """
    Perform dimensionality reduction using scBiMapping's spectral embedding approach.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix with cells (rows) and features (columns)
    n_embedding : int, optional
        Dimension of the embedding space (default: 30)
    normalization : bool, optional
        Whether to L2-normalize the embeddings (default: True)
        
    Returns:
    --------
    AnnData
        Modified AnnData object with embeddings added to:
        - obsm['U']: Cell embeddings (n_cells x n_embedding)
        - uns['V']: Feature embeddings (n_features x n_embedding)
    """
    
    # Ensure sparse matrix format   
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
        
    eps = 1e-10 # Numerical stability constant
    
    Dx = diags(np.ravel(1/(adata.X.sum(axis = 1)+eps)))  # actually Dx^(-1): Inverse of row sums (cells)
    Dy_sqrt = diags(sqrt(np.ravel(1/(adata.X.sum(axis = 0)+eps))))  # actually Dy^(-1): inverse column sums (features)
    C = sqrt(Dx)@adata.X@Dy_sqrt
    _,singularVals,evec = sparse.linalg.svds(C, k=n_embedding,return_singular_vectors = "vh",random_state =0) # evec (k*Ny): Unitary matrix having right singular vectors as rows.
    adata.uns['V'] = Dy_sqrt@evec.T; # Feature embeddings
    adata.obsm['U'] = Dx@(adata.X@(adata.uns['V']@diags(1/singularVals))); # Cell embeddings 
         
    adata.obsm['U'] = normalize(adata.obsm['U'], axis=1, norm='l2');
    adata.uns['V'] = normalize(adata.uns['V'], axis=1, norm='l2')
    return adata

###########################################################################
#    reference-dataset based cell annotation(scBiMapping_annotation)      #
###########################################################################

def scBiMapping_annotation(
    adata_ref,
    adata_query,
    n_embedding: int = 30,
    K: int = 30,
    K_majority: int = 10,
    knnMethod: str = 'HNSW',
    normalization: bool = True,
    reduction_method_on_cells_only: str = 'BiMapping',
    metric: str = 'euclidean',
    n_embedding_2nd: Optional[int] = None,
    CellType_Key_for_ref: Union[str, List[str]] = 'cell_annotation'
) -> List[str]:
    """
    Perform cell type annotation using scBiMapping approach.
    
    Parameters:
    -----------
    adata_ref : AnnData
        Reference dataset with known cell types
    adata_query : AnnData
        Query dataset to be annotated
    n_embedding : int
        Dimension of the embedding space (default: 30)
    K : int
        Number of nearest neighbors to find (default: 30)
    K_majority : int
        Number of neighbors for majority voting (default: 10)
    knnMethod : str
        KNN search method ('HNSW' or 'NNDescent') (default: 'HNSW')
    normalization : bool
        Whether to normalize embeddings (default: True)
    reduction_method_on_cells_only : str
        Dimensionality reduction method ('BiMapping' or 'None') (default: 'BiMapping')
    metric : str
        Distance metric for KNN search (default: 'euclidean')
    n_embedding_2nd : Optional[int]
        Alternative embedding dimension if needed (default: None)
    CellType_Key_for_ref : Union[str, List[str]]
        Key(s) for cell type annotations in adata_ref.obs (default: 'cell_annotation')
        
    Returns:
    --------
    List[str]
        Predicted cell types for query cells
    """
    
    # Print parameters
    print(f'n_embedding: {n_embedding}')
    print(f'K: {K}')
    print(f'K_majority: {K_majority}')
 
    # Ensure integer parameters
    K = int(K)
    K_majority = int(K_majority)
    
    if adata_ref.n_vars != adata_query.n_vars:
        intersection_feature = list(set(adata_ref.var_names) & set(adata_query.var_names))
        adata_ref = adata_ref[:,intersection_feature]
        adata_query = adata_query[:,intersection_feature]
    
    # Step 1: Compute embeddings for reference and query
    scBiMapping_DR(adata_ref,n_embedding = n_embedding)
    scBiMapping_DR(adata_query,n_embedding = n_embedding)
    
    # Step 2: Compute similarity matrices and merge
    knn_based_Sim_ref = knn_search_weighted(adata_ref, K, knnMethod,metric)
    knn_based_Sim_query = knn_search_weighted(adata_query, K, knnMethod,metric) 
    
    print('\nDirect merge softmax-weight coded reference and query dataset......')
    knn_based_Sim_both_ref_query = sparse.vstack((knn_based_Sim_ref, knn_based_Sim_query)) 
    
    # Step 3: Dimensionality reduction
    print(f'\nPerforming dimensionality reduction using {reduction_method_on_cells_only}...')
    if reduction_method_on_cells_only == 'BiMapping':
        n_embedding_final = n_embedding_2nd if n_embedding_2nd is not None else n_embedding
        U = BiTcut_embedding_v5_matrix_form(
            knn_based_Sim_both_ref_query,
            n_embedding=n_embedding_final,
            normalization=normalization
        )
    elif reduction_method_on_cells_only == 'None':
        U = knn_based_Sim_both_ref_query
    else:
        raise ValueError(f"Unsupported reduction method: {reduction_method_on_cells_only}")
                                    
    adata_ref.obsm['knn_based_Sim_ref'] = U[:knn_based_Sim_ref.shape[0],:]
    adata_query.obsm['knn_based_Sim_query'] = U[knn_based_Sim_ref.shape[0]:,:]                                    
    print(f'Reference shape: {adata_ref.obsm["knn_based_Sim_ref"].shape}')
    print(f'Query shape: {adata_query.obsm["knn_based_Sim_query"].shape}')
    
     
    # Step 5: Find nearest neighbors and predict cell types by majority voting
    if reduction_method_on_cells_only == 'None':
        print('when reduction_method_on_cells_only is None, NNDescent is recommended, which supports sparse input')   
    nnIndex, dis = knn_search(query_dataset = adata_query.obsm['knn_based_Sim_query'], index_dataset = adata_ref.obsm['knn_based_Sim_ref'], K = K_majority, knn_method = knnMethod)
    
    print('Performing majority voting...')
    if isinstance(CellType_Key_for_ref,str):  
        celltype_predicted = pd.DataFrame(np.array(adata_ref.obs[CellType_Key_for_ref])[nnIndex]).apply(find_most_frequent_np, axis=1).tolist()
        adata_query.obs['cell_type_predicted'] = celltype_predicted
    elif isinstance(CellType_Key_for_ref,list):
        for i in range(len(CellType_Key_for_ref)):
            celltype_predicted = pd.DataFrame(np.array(adata_ref.obs[CellType_Key_for_ref[i]])[nnIndex]).apply(find_most_frequent_np, axis=1).tolist()
            adata_query.obs['cell_type_predicted' + str(i+1)] = celltype_predicted
    else:
        raise TypeError("Unsupported type for count_elements function")
        
    return celltype_predicted

def scBiMapping_annotation_v2(
    adata_ref,
    adata_query,
    n_embedding: int = 30,
    K: int = 30,
    K_majority: int = 10,
    knnMethod: str = 'HNSW',
    normalization: bool = True,
    reduction_method_on_cells_only: str = 'BiMapping',
    metric: str = 'euclidean',
    n_embedding_2nd: Optional[int] = None,
    CellType_Key_for_ref: Union[str, List[str]] = 'cell_annotation'
) -> List[str]:
    """
    Perform cell type annotation using scBiMapping approach.
    
    Parameters:
    -----------
    adata_ref : AnnData
        Reference dataset with known cell types
    adata_query : AnnData
        Query dataset to be annotated
    n_embedding : int
        Dimension of the embedding space (default: 30)
    K : int
        Number of nearest neighbors to find (default: 30)
    K_majority : int
        Number of neighbors for majority voting (default: 10)
    knnMethod : str
        KNN search method ('HNSW' or 'NNDescent') (default: 'HNSW')
    normalization : bool
        Whether to normalize embeddings (default: True)
    reduction_method_on_cells_only : str
        Dimensionality reduction method ('BiMapping' or 'None') (default: 'BiMapping')
    metric : str
        Distance metric for KNN search (default: 'euclidean')
    n_embedding_2nd : Optional[int]
        Alternative embedding dimension if needed (default: None)
    CellType_Key_for_ref : Union[str, List[str]]
        Key(s) for cell type annotations in adata_ref.obs (default: 'cell_annotation')
        
    Returns:
    --------
    List[str]
        Predicted cell types for query cells
    """
    
    # Print parameters
    print(f'n_embedding: {n_embedding}')
    print(f'K: {K}')
    print(f'K_majority: {K_majority}')
 
    # Ensure integer parameters
    K = int(K)
    K_majority = int(K_majority)
    
    if adata_ref.n_vars != adata_query.n_vars:
        intersection_feature = list(set(adata_ref.var_names) & set(adata_query.var_names))
        adata_ref = adata_ref[:,intersection_feature]
        adata_query = adata_query[:,intersection_feature]
        
    # Step 1: Compute embeddings for reference and query
    scBiMapping_DR(adata_ref,n_embedding = n_embedding)
    scBiMapping_DR(adata_query,n_embedding = n_embedding)
    
    # Step 2: Compute similarity matrices and merge
    knn_based_Sim_ref = knn_search_weighted(adata_ref, K, knnMethod,metric)
    knn_based_Sim_query = knn_search_weighted(adata_query, K, knnMethod,metric) 
    

    # Step 3: Dimensionality reduction
    print(f'\nPerforming dimensionality reduction using {reduction_method_on_cells_only}...')
    if reduction_method_on_cells_only == 'BiMapping':
        print('\nDirect merge softmax-weight coded reference and query dataset......')
        knn_based_Sim_both_ref_query = sparse.vstack((knn_based_Sim_ref, knn_based_Sim_query)) 
        n_embedding_final = n_embedding_2nd if n_embedding_2nd is not None else n_embedding
        U = BiTcut_embedding_v5_matrix_form(
            knn_based_Sim_both_ref_query,
            n_embedding=n_embedding_final,
            normalization=normalization
        )
        adata_ref.obsm['knn_based_Sim_ref'],adata_query.obsm['knn_based_Sim_query'] = U[:knn_based_Sim_ref.shape[0],:],U[knn_based_Sim_ref.shape[0]:,:]  
        nnIndex, dis = knn_search(query_dataset = adata_query.obsm['knn_based_Sim_query'], index_dataset = adata_ref.obsm['knn_based_Sim_ref'], K = K_majority, knn_method = knnMethod)
    elif reduction_method_on_cells_only == 'None': 
        print('when reduction_method_on_cells_only is None, NNDescent is recommended, which supports sparse input')
        knnMethod = 'NNDescent'
        nnIndex, dis = knn_search(query_dataset = knn_based_Sim_query, index_dataset = knn_based_Sim_ref, K = K_majority, knn_method = knnMethod)
    else:
        raise ValueError(f"Unsupported reduction method: {reduction_method_on_cells_only}")
                                    
                                 
    print(f'Reference shape: {adata_ref.obsm["knn_based_Sim_ref"].shape}')
    print(f'Query shape: {adata_query.obsm["knn_based_Sim_query"].shape}')
    
     
    # Step 5: Find nearest neighbors and predict cell types by majority voting 
    print('Performing majority voting...')
    if isinstance(CellType_Key_for_ref,str):  
        celltype_predicted = pd.DataFrame(np.array(adata_ref.obs[CellType_Key_for_ref])[nnIndex]).apply(find_most_frequent_np, axis=1).tolist()
        adata_query.obs['cell_type_predicted'] = celltype_predicted
    elif isinstance(CellType_Key_for_ref,list):
        for i in range(len(CellType_Key_for_ref)):
            celltype_predicted = pd.DataFrame(np.array(adata_ref.obs[CellType_Key_for_ref[i]])[nnIndex]).apply(find_most_frequent_np, axis=1).tolist()
            adata_query.obs['cell_type_predicted' + str(i+1)] = celltype_predicted
    else:
        raise TypeError("Unsupported type for count_elements function")
        
    return celltype_predicted

def scBiMapping_batch_correction(
    adata: ad.AnnData,
    batch_key: str = 'batch',
    n_embedding: int = 30,
    K: int = 30, 
    knn_method: str = 'HNSW',
    normalization: bool = True,
    reduction_method: str = 'BiMapping',
    metric: str = 'euclidean',
    n_embedding_2nd: Optional[int] = None
) -> ad.AnnData:
    """
    Perform batch effect correction on multi-batch single-cell data using scBiMapping.
    
    This function integrates multiple batches by learning a joint embedding space while 
    preserving biological variations and removing technical batch effects.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix containing cells from multiple batches
    batch_key : str
        Key in adata.obs that identifies different batches (default: 'batch')
    n_embedding : int
        Dimension of the embedding space (default: 30)
    K : int
        Number of nearest neighbors to use for graph construction (default: 30)
    knn_method : str
        KNN search method, either 'HNSW' or 'NNDescent' (default: 'HNSW')
    normalization : bool
        Whether to L2-normalize the embeddings (default: True)
    reduction_method : str
        Dimensionality reduction method ('BiMapping' or 'None') (default: 'BiMapping')
    metric : str
        Distance metric for KNN search (default: 'euclidean')
    n_embedding_2nd : Optional[int]
        Alternative embedding dimension for second reduction step (default: None)

    Returns:
    --------
    AnnData
        Returns the modified AnnData object with batch-corrected embeddings stored in:
        - obsm['batch_corrected_embeddings']: Joint embeddings after batch correction
    """
    
    # Validate input parameters
    if batch_key not in adata.obs:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")
    
    print(f"\nStarting batch correction with parameters:")
    print(f"- Embedding dimensions: {n_embedding} (primary), {n_embedding_2nd} (secondary)")
    print(f"- Nearest neighbors: K={K}")
    print(f"- KNN method: {knn_method}")
    print(f"- Reduction method: {reduction_method}\n")
    
    # Initialize container for merged similarity matrices
    knn_based_sim_list = []
    batch_labels = adata.obs[batch_key].unique()
    
    # Process each batch separately
    for i, batch in enumerate(batch_labels):
        print(f"Processing batch {i+1}/{len(batch_labels)}: {batch}")
        
        # Subset the current batch
        batch_mask = adata.obs[batch_key] == batch
        adata_batch = adata[batch_mask, :].copy()
        
        # Compute batch-specific embeddings and similarity matrix
        scBiMapping_DR(adata_batch, n_embedding=n_embedding)
        sim_matrix = knn_search_weighted(adata_batch, K=K, knn_method=knn_method, metric=metric)
        knn_based_sim_list.append(sim_matrix)
        
        # Clean up memory
        del adata_batch, sim_matrix
    
    # Combine similarity matrices from all batches
    print("\nMerging similarity matrices from all batches...")
    knn_based_sim_combined = sparse.vstack(knn_based_sim_list)
    
    # Perform final dimensionality reduction
    print(f"\nPerforming final dimensionality reduction using {reduction_method}...")
    if reduction_method == 'BiMapping':
        n_embedding_final = n_embedding_2nd if n_embedding_2nd is not None else n_embedding
        corrected_embeddings = BiTcut_embedding_v5_matrix_form(
            knn_based_sim_combined,
            n_embedding=n_embedding_final,
            normalization=normalization
        )
    elif reduction_method == 'None':
        corrected_embeddings = knn_based_sim_combined
    else:
        raise ValueError(f"Unsupported reduction method: {reduction_method}")
    
    # Store results in the original AnnData object
    adata.obsm['batch_corrected_embeddings'] = corrected_embeddings
    print("\nBatch correction completed successfully.")
    
    return adata
    
def knn_search_weighted(adata, K, knn_method='HNSW',metric = 'euclidean'):
    """
    Find K nearest genes for each cell in co-embedded space with weighted similarity
    
    Parameters:
    -----------
    adata : AnnData
        Contains 'U' (cell embeddings) and 'V' (gene embeddings)
    K : int
        Number of nearest neighbors
    knn_method : str
        KNN search method
    metric : str
        Distance metric
        
    Returns:
    --------
    csr_matrix
        Weighted similarity matrix between cells and genes
    """

    print('for each cell, find K nearest genes in the co-embedded space......')
 
    # Get nearest neighbors
    nnIndex, dis = knn_search(adata.obsm['U'], adata.uns['V'], K, knn_method, metric)
    
    n = adata.obsm['U'].shape[0]
    m = adata.uns['V'].shape[0]
    
    # Convert distances to normalized similarities
    if metric == 'cosine':
        dis = 1 - dis  # Convert cosine distance to similarity
    else:
        dis = exp(-dis)  # Convert Euclidean distance to similarity
    
    # Normalize similarities per cell
    row_sums = dis.sum(axis=1, keepdims=True)
    dis = dis / (row_sums + 1e-15)
    
    # Create sparse similarity matrix
    I1 = tile(arange(n), K)
    J1 = nnIndex.flatten('F')
    Sim = dis.flatten('F')
    
    return coo_matrix((Sim, (I1, J1)), shape=(n, m)).tocsr()

def knn_search(query_dataset, index_dataset, K, knn_method='HNSW', metric='euclidean'):
    """
    Perform K-nearest neighbors search using specified method
    
    Parameters:
    -----------
    query_data : ndarray
        Query data points (n_samples x n_features)
    index_data : ndarray
        Index data points (n_samples x n_features) 
    K : int
        Number of neighbors to find
    knn_method : str (default: 'HNSW')
        Search method ('HNSW', 'NNDescent', or 'bruteforce')
    metric : str (default: 'euclidean')
        Distance metric ('euclidean', 'cosine', 'ip')
        
    Returns:
    --------
    nnIndex : ndarray
        Indices of nearest neighbors
    distances : ndarray
        Distances to nearest neighbors
    """

    print('find knn...')
    N, Dim = index_dataset.shape
    if knn_method == 'HNSW':
        # Map metric to HNSW space parameter
        space_map = {
            'euclidean': 'l2',
            'l2': 'l2',
            'cosine': 'cosine',
            'ip': 'ip'
        }
        
        if metric not in space_map:
            raise ValueError(f"Unsupported metric '{metric}' for HNSW")
    
        # Declaring index
        random.seed(1)
        # possible options are l2, cosine or ip
        p = hnswlib.Index(space=space_map[metric], dim=Dim)
        # Initializing index - the maximum number of elements should be known beforehand
        random.seed(1)
        p.init_index(max_elements=N, ef_construction=200,
                     M=16, random_seed=100)
        # Element insertion (can be called several times):
        random.seed(1)
        p.add_items(index_dataset, arange(N))
        # Controlling the recall by setting ef:
        p.set_ef(50)  # ef should always be > k
        # Query dataset, k - number of closest elements (returns 2 numpy arrays)
        random.seed(1)
        nnIndex, dis = p.knn_query(query_dataset, k=K)
        dis = sqrt(dis) if space_map[metric] == 'l2' else dis
    elif knn_method == 'NNDescent': # Find k nearest neighbors based on NNDescent
        # initialized by random projection trees; n_neighbors is a parameter in indexing for constructing the searching graph
        random.seed(1)
        index = NNDescent(index_dataset, n_neighbors=K, random_state=1) # n_jobs is avaible for this function
        nnIndex, dis = index.query(query_dataset, k=K)
    elif knn_method == 'bruteforce':
        print('knnMethod: ' + knn_method)
        
        neigh = NearestNeighbors(n_neighbors=2,algorithm='brute')
        neigh.fit(index_dataset)
        dis,nnIndex= neigh.kneighbors(query_dataset, 2, return_distance=True)
    else:
        print('No such knn method............................................')

    return nnIndex, dis

def BiTcut_embedding_v5_matrix_form(X,n_embedding,normalization = True): 

    """
    Compute embeddings for matrix X
    
    Parameters:
    -----------
    X : array or sparse matrix
        Input data matrix (cells x genes)
    n_embedding : int
        Number of embedding dimensions
    normalization : bool
        Whether to L2-normalize embeddings
        
    Returns:
    --------
    ndarray
        Cell embeddings (n_cells x n_embedding)
    """
   
    if not issparse(X):
        X = csr_matrix(X)
        
    eps = 1e-10
    Dx = diags(np.ravel(1/(X.sum(axis = 1)+eps)))  
    Dy_sqrt = diags(sqrt(np.ravel(1/(X.sum(axis = 0)+eps)))) 
    C = sqrt(Dx)@X@Dy_sqrt
 
    _,vals,evec = sparse.linalg.svds(C, k=n_embedding,return_singular_vectors = "vh",random_state =0) # evec (k*Ny): Unitary matrix having right singular vectors as rows.
 
    U = Dx@(X@(Dy_sqrt@evec.T@diags(1/vals)));
    if normalization:
        U = normalize(U, axis=1, norm='l2');   # normalize each row to unit norm
     
    return U


def find_most_frequent_np(row):
    """
    Function: majority voting implementation 
    Parameters:
      row : array
    Returns:
      Most frequent element (returns first occurrence in case of ties)
    """
    
    uni, counts = np.unique(row, return_counts=True)
    return uni[np.argmax(counts)]

 


