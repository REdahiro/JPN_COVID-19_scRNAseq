import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt
import sys
args = sys.argv

thread = int(sys.argv[1])
min_size = int(sys.argv[2])

#-------------------------------------------------------------------------------------#
# dir: /work22/home/redahiro/analysis/COVID-19_scRNAseq/TCR_BCR/BCR
# env: scirpy
#-------------------------------------------------------------------------------------#

# Import h5ad file
adata = sc.read("BCR.h5ad")               # 123,728
print(adata.shape)


# metadata edit
adata.obs.dtypes
print(adata.obs["l3"].dtype)   # Pandasのdata型 check    #int32 → categoryに変更
adata.obs["l3"] = adata.obs["l3"].astype('category')
print(adata.obs["orig.ident"].dtype)
adata.obs["l3"].value_counts()

#----- Rename categories -----#
adata.obs['l2'] = adata.obs['l3'].apply(lambda x : 'B_Naive' if x == 0.0 else ('B_Memory' if x == 1.0 else('B_Activated' if x == 2.0 else 'Plasma')))

#----- Categorical order -----#
l2_order = ('B_Naive','B_Memory','B_Activated','Plasma')
adata.obs['l2'] = adata.obs['l2'].cat.reorder_categories(list(l2_order),ordered=True) 

print(adata.obs['l2'].dtype) 
adata.obs['l2'].value_counts()


#-----------------------------#
#    TCR Quality Control      #
#-----------------------------#

ir.tl.chain_qc(adata)

# Remove multi_chain Barcode
adata = adata[adata.obs["multi_chain"] != "True", :].copy()    
print("Number of Cell after remove multi_chain:{:.2f}".format(adata.n_obs))
# Remove all cells that don’t have at least one full pair of receptor sequences
adata = adata[~adata.obs["chain_pairing"].isin(["orphan VDJ", "orphan VJ", "ambiguous","no IR"]), :].copy()
print("Number of Cell after remove orpha_chain:{:.2f}".format(adata.n_obs))                    

print("data shape after filtering")
print(adata.shape)

#------------------------------------------------#
#    Define clonotypes and clonotype clusters    #
#------------------------------------------------#

#----- Construct a neighborhood graph based on CDR3 nucleotide sequency similariy ----#

ir.pp.ir_dist(adata)     # computes sequence-distance metric btw all unique CDR3
                         # default: using identity (not alighment), sequecen: nucleotide

# detect connected modules in the graph and annotate them as clonotypes
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only", n_jobs=thread)

ir.tl.clonotype_network(adata, min_cells = 2)

#----- Recompute CDR3 neighborhood graph and difine clonotype clusters -----#
ir.pp.ir_dist(
	    adata, 
	    metric = "alignment",
	    sequence = "aa",
	    cutoff = 10,
	    n_jobs=thread
        )

ir.tl.define_clonotype_clusters(
	   adata, sequence = "aa", metric = "alignment", receptor_arms="all", dual_ir="any", n_jobs=thread
	   )

ir.tl.clonotype_network(adata, min_cells = min_size , sequence = "aa", metric = "alignment")


##-----------------------##
##        Export         ##
##-----------------------## 
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata.write_h5ad(filename = f'BCR_ClonotypeSize_{min_size}.h5ad', compression ='gzip')


