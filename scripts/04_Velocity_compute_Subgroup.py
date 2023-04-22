import scvelo as scv
import anndata
import pandas as pd
import numpy as np
import matplotlib as plt
from matplotlib import pyplot as plt
import sys
args = sys.argv

INPUT=args[1]
Status=args[2]
name_status=args[3]
n_genes = int(sys.argv[4])
n_neigh = int(sys.argv[5])
OUTPUT=args[6]

# import data
adata = scv.read(INPUT)

# Making categorical data
adata.obs["Status"] = adata.obs["Status"].astype('category')
adata.obs["Severity"] = adata.obs["Severity"].astype('category')
adata.obs["l1"] = adata.obs["l1"].astype('category')
adata.obs["l2"] = adata.obs["l2"].astype('category')
adata.obs["l3"] = adata.obs["l3"].astype('category')
# print(adata.obs["Status"].dtype)   # object to category

adata = adata[adata.obs[Status].isin([name_status]),:]          # pickup subgroup

print("data shape")
print(adata.shape)

#--------------------------------#
#      Preprocessing(PP)         #
#--------------------------------#
# Filtering, normalization and log-transformation
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=n_genes)

# Computes moments for velocity estimation: 
scv.pp.moments(adata, n_pcs=30, n_neighbors=n_neigh)             # number of neighbors to use : default = 30


#--------------------------------#
#          Tools(tl)             #
#--------------------------------#
# Estimate full splicing kinetics of specified genes
# skipping this process in state_state mode 
scv.tl.recover_dynamics(adata, var_names = 'all', n_jobs=8)

# Estimates velocities in a gene-specfic manner
scv.tl.velocity(adata, mode='dynamical')

print("Number os Velocity_genes")
print(adata.var.velocity_genes.sum())


# Computes velocity graph based on cosine similarities
scv.tl.velocity_graph(adata)

#----------------------------------#
#      Project the velocities      #
#----------------------------------#
#--- stream plot of velocities on the embedding ---#
# min_mass: 
scv.pl.velocity_embedding_stream(
	adata, basis="umap", color="l2", min_mass=1, 
	figsize=(10,10), dpi=300, save=f'Velocity_stream_1_default_{n_genes}.png')

scv.pl.velocity_embedding_stream(
	adata, basis="umap", color="l2", min_mass=3, 
	figsize=(10,10), dpi=300, save=f'Velocity_stream_3_{n_genes}.png')

scv.pl.velocity_embedding_stream(
	adata, basis="umap", color="l2", min_mass=3.5, 
	figsize=(10,10), dpi=300, save=f'Velocity_stream_3.5_{n_genes}.png')

# Scatter plot of velocities on the embedding
scv.pl.velocity_embedding(
	adata, basis="umap", color="l2", arrow_length=5, arrow_size=2,
	figsize=(10,10), dpi=300, save=f'Velocity_vector_2_{n_genes}.png')



#------------------------------------------#
#            save the data                 #
#------------------------------------------#
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata.write_h5ad(filename = f'{OUTPUT}_{n_genes}.h5ad', compression ='gzip')

