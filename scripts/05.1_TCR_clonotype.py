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
# dir: /work22/home/redahiro/analysis/COVID-19_scRNAseq/TCR_BCR/TCR
# env: scirpy
#-------------------------------------------------------------------------------------#

#Import h5ad file
adata = sc.read("TCR.h5ad")                               # 628,715 cells
print(adata.shape)


# metadata edit
adata.obs.dtypes
print(adata.obs["l3"].dtype)
adata.obs["l3"] = adata.obs["l3"].astype('category')
print(adata.obs["orig.ident"].dtype)
adata.obs["l3"].value_counts()

#----- Rename categories -----#
adata.obs['l3'] = adata.obs['l3'].apply(lambda x : 'CD4_Naive' if x == 0.0 
	                                             else ('CD4_Memory' if x == 1.0 
	                                             else('CD4_Ef' if x == 2.0
	                                             else('Treg' if x == 3.0 
	                                             else('CD8_Naive' if x == 4.0 
	                                             else('CD8_Memory' if x == 5.0 
	                                             else('CD8_Ef' if x == 6.0 
	                                             else('NKT' if x == 7.0 
	                                             else('MAIT' if x == 8.0 
	                                             else('gdT' if x == 9.0 
	                                             else('Pro_T' if x == 10.0 
	                                             else('NK' if x == 11.0  
	                                             else 'NK_CD56bright'))))))))))))

#----- Categorical order -----#
l3_order = ('CD4_Naive','CD4_Memory','CD4_Ef','Treg','CD8_Naive','CD8_Memory','CD8_Ef',
            'NKT','MAIT','gdT','Pro_T','NK','NK_CD56bright')
adata.obs['l3'] = adata.obs['l3'].cat.reorder_categories(list(l3_order),ordered=True) 

print(adata.obs['l3'].dtype) 
adata.obs['l3'].value_counts()

# Status
adata.obs["Status"] = adata.obs["Status"].astype('category')
adata.obs['Status'] = adata.obs['Status'].apply(lambda x : 'HC' if x == 0.0 else 'COVID19')

adata.obs["Severity"] = adata.obs["Severity"].astype('category')
adata.obs['Severity'] = adata.obs['Severity'].apply(lambda x : 'HC' if x == 0.0 else ('Moderate' if x == 1.0 else 'Severe'))
Severity_order = ('HC','Moderate','Severe')
adata.obs['Severity'] = adata.obs['Severity'].cat.reorder_categories(list(Severity_order),ordered=True) 


##---------------------------------------------------##
##           select analysis T clusters              ##
##---------------------------------------------------##

# remove NK, NK_CD56bright, NKT, gdT
adata = adata[adata.obs["l3"].isin(['CD4_Naive','CD4_Memory','CD4_Ef','Treg','CD8_Naive','CD8_Memory','CD8_Ef','MAIT','Pro_T']),:]
print("data shape after pick up analysis clusters")
print(adata.shape)


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
adata.write_h5ad(filename = f'TCR_ClonotypeSize_{min_size}.h5ad', compression ='gzip')


