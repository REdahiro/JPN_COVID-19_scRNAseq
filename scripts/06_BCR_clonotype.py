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
# 個数の確認
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
# receptor_arms = "all": both TRA and TRB need to match
# dual_ir = "primary_only": only consider most abundant pair of TRA/TRB chains

ir.pp.ir_dist(adata)     # computes sequence-distance metric btw all unique CDR3
                         # default: using identity (not alighment), sequecen: nucleotide

# detect connected modules in the graph and annotate them as clonotypes
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only", n_jobs=thread)

ir.tl.clonotype_network(adata, min_cells = 2)        # default: sequence="nt", default="identity"

#ax = ir.pl.clonotype_network(adata, color="orig.ident", base_size=20, label_fontsize=9, panel_size=(7, 7))
#fig = ax.get_figure()
#fig.savefig("figures/clonotype_network_perSample_nt.png", dpi = 300, bbox_inches="tight")
#AttributeError: module 'matplotlib.cbook' has no attribute 'iterable' のため
 


#----- Recompute CDR3 neighborhood graph and difine clonotype clusters -----#
ir.pp.ir_dist(
	    adata, 
	    metric = "alignment",
	    sequence = "aa",
	    cutoff = 10,                     # 距離の値、大きいほど大きなclusterが出来上がる
	    n_jobs=thread
        )


ir.tl.define_clonotype_clusters(
	   adata, sequence = "aa", metric = "alignment", receptor_arms="all", dual_ir="any", n_jobs=thread
	   )


ir.tl.clonotype_network(adata, min_cells = min_size , sequence = "aa", metric = "alignment")


# ir.pl.clonotype_network(adata, color="orig.ident", label_fontsize=9, panel_size=(7, 7), base_size=20)
# AttributeError: module 'matplotlib.cbook' has no attribute 'iterable'


##-----------------------##
##        Export         ##
##-----------------------## 
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata.write_h5ad(filename = f'BCR_ClonotypeSize_{min_size}.h5ad', compression ='gzip')