import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt
from cycler import cycler

#-------------------------------------------------------------------------------------#
# dir: /work22/home/redahiro/analysis/COVID-19_scRNAseq/TCR_BCR/TCR
# env: scirpy
#-------------------------------------------------------------------------------------#

plt.rcParams['figure.figsize'] = 4,4
plt.rcParams['figure.dpi'] = 300
plt.xlabel("xlabel",fontsize=8)

sc.set_figure_params(figsize=(8, 8))
sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)

# VDJ
vdjdb = sc.read("/work22/home/redahiro/software/scRNAseq_script/TCR_BCR/vdjdb.h5ad")
#vdjdb = ir.datasets.vdjdb()        # 55,586 

# TCR data
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

# TCR quality check

ir.tl.chain_qc(adata)

# Remove multi_chain Barcode
adata = adata[adata.obs["multi_chain"] != "True", :].copy()    
print("Number of Cell after remove multi_chain:{:.2f}".format(adata.n_obs))
# Remove all cells that don’t have at least one full pair of receptor sequences
adata = adata[~adata.obs["chain_pairing"].isin(["orphan VDJ", "orphan VJ", "ambiguous","no IR"]), :].copy()
print("Number of Cell after remove orpha_chain:{:.2f}".format(adata.n_obs))                    

print("data shape after filtering")
print(adata.shape)


#---------------------------------------------------------------#
#                     Epitope analysis                          #
#---------------------------------------------------------------#

# Compute sequence-based distance matrices
ir.pp.ir_dist(
    adata, vdjdb, 
    metric="alignment",     
    sequence="aa", 
    cutoff = 10,
    n_jobs=4
)            


# Identify matching entries in a reference database for each cell
ir.tl.ir_query(
    adata, vdjdb, 
    metric="alignment", sequence="aa", 
    receptor_arms="any", dual_ir="any",
    n_jobs=4
)

# Return a dataframe with all matching annotations
ir.tl.ir_query_annotate_df(
    adata,
    vdjdb,
    metric="alignment",
    sequence="aa",
    include_ref_cols=["antigen.species", "antigen.gene"]
)

ir.tl.ir_query_annotate(
    adata,
    vdjdb,
    metric="alignment",
    sequence="aa",
    include_ref_cols=["antigen.species", "antigen.gene"],
    strategy="most-frequent"
)

#-----------------------#
#         UAMP          #
#-----------------------# 
sc.pl.umap(adata, color="antigen.species", save='clonal_expansion_antigen.png')

sc.pl.umap(
    adata, 
    color = "antigen.species", 
    #groups = ["1","2",">= 3"],
    size = [
          50 if x in ["SARS-CoV-2"] else 0.5
          for x in adata.obs["antigen.species"]
          ],
    save='clonal_expansion_antigen.png')

adata.obs['SARS-CoV-2'] = adata.obs['antigen.species'].apply(lambda x : "SARS-CoV-2" if x == 'SARS-CoV-2' else "Non_SARS-CoV-2")
adata.obs['SARS-CoV-2'] = adata.obs['SARS-CoV-2'].astype('category') 
adata.obs['SARS-CoV-2'] = adata.obs['SARS-CoV-2'].fillna("Non_SARS-CoV-2")

sc.pl.umap(
    adata, 
    color = "SARS-CoV-2", 
    #groups = ["1","2",">= 3"],
    size = [
          15 if x in ["SARS-CoV-2"] else 0.5
          for x in adata.obs["SARS-CoV-2"]
          ],
    save='clonal_expansion_SARS-CoV-2.png')


#-----------------------------#
#          Export data        #
#-----------------------------#

# h5ad
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata.write_h5ad(filename = 'TCR_filtered_EpitopeInfo.h5ad', compression ='gzip')

# metadata
meta_data = adata.obs.loc[:,['antigen.species','antigen.gene','l3']] 
meta_data.to_csv('TCR_filtered_EpitopeInfo.txt', sep='\t', index=True, index_label='cellID')


