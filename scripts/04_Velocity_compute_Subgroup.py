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
r2_min = float(sys.argv[6])
likelihood_min = float(sys.argv[7])
OUTPUT=args[8]

# env: Seurat4.1 (Seurat4.2 does not work due to pandas version)

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
# adata = adata[adata.obs["status"].isin(["case"]),:] 

print("data shape")
print(adata.shape)

#--------------------------------#
#      Preprocessing(PP)         #
#--------------------------------#

# Filtering, normalization and log-transformation
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=n_genes)   # SeuratでNormalization行う必要なし

# Computes moments for velocity estimation: 
scv.pp.moments(adata, n_pcs=30, n_neighbors=n_neigh)             # number of neighbors to use : default = 30


#--------------------------------#
#          Tools(tl)             #
#--------------------------------#

# Estimate full splicing kinetics of specified genes
# skipping this process in state_state mode 
scv.tl.recover_dynamics(adata, var_names = 'all', n_jobs=8)                    # var_names = "velocity_genes" is default (https://github.com/theislab/scvelo/issues/129)

# Estimates velocities in a gene-specfic manner: too few genes wre selected as velocity genes　→ setting a lower threshold of min_r2 as suggested
# Define Velocity genes: taken into account for all downstream analysis
scv.tl.velocity(adata, mode='dynamical', min_r2=r2_min, min_likelihood=likelihood_min)     # min_r2=0.01, min_likelihood=0.001 

print("Number os Velocity_genes")
print(adata.var.velocity_genes.sum())
#adata.var['velocity_genes'] = True      # enforce all genes to velocity_gene


# Computes velocity graph based on cosine similarities: nobs*nobs matrix 
scv.tl.velocity_graph(adata)


#----------------------------------#
#      Project the velocities      #
#----------------------------------#

#--- stream plot of velocities on the embedding ---#
# min_mass: 
scv.pl.velocity_embedding_stream(
	adata, basis="umap", color="l2", min_mass=1, 
	figsize=(10,10), dpi=300, save=f'Velocity_stream_1_default_{n_genes}_{r2_min}_{likelihood_min}.png')

scv.pl.velocity_embedding_stream(
	adata, basis="umap", color="l2", min_mass=3, 
	figsize=(10,10), dpi=300, save=f'Velocity_stream_3_{n_genes}_{r2_min}_{likelihood_min}.png')

scv.pl.velocity_embedding_stream(
	adata, basis="umap", color="l2", min_mass=3.5, 
	figsize=(10,10), dpi=300, save=f'Velocity_stream_3.5_{n_genes}_{r2_min}_{likelihood_min}.png')

# Scatter plot of velocities on the embedding
scv.pl.velocity_embedding(
	adata, basis="umap", color="l2", arrow_length=5, arrow_size=2,
	figsize=(10,10), dpi=300, save=f'Velocity_vector_2_{n_genes}_{r2_min}_{likelihood_min}.png')



#------------------------------------#
#          Speed & coherence         #
#------------------------------------#
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence' 

scv.pl.scatter(adata, c=keys, cmap='coolwarm', size=7.5, perc=[5,95], figsize=(5,5), dpi=300, save='speed_coherence.png')



#------------------------------------------#
#            save the data                 #
#------------------------------------------#

# https://github.com/theislab/scvelo/issues/255 を参考にして以下のコマンドで adataの修正した上で、h5adでsaveする
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata.write_h5ad(filename = f'{OUTPUT}_{n_genes}_{r2_min}_{likelihood_min}.h5ad', compression ='gzip')
#adata.write_h5ad(filename = 'Innate_with_velocity_dynamical.h5ad', compression ='gzip')

# https://github.com/theislab/scvelo/issues/200 !!!!!
# https://github.com/theislab/scvelo/issues/158 !!!!!
