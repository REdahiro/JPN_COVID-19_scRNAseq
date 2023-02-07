library(Seurat)
library(tidyverse)

data <- readRDS("./Data/ALL_RNA.rds")

cluster_l2 <- c("CD4T","Treg","CD8T","MAIT","Pro_T","NK",
               "CD14p_Mono","CD16p_Mono","cDC","pDC","B","Plasma","Platelet")

data <- SetIdent(data, value = "l2") 
levels(data) <- cluster_l2

gene <- rownames(data)
module = "GOBP"

#--------------------------------#
#            Function            #
#--------------------------------#

for(list in c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON","GOBP_RESPONSE_TO_INTERFERON_GAMMA")){

   module_gene <- read.table(str_c("/work22/home/redahiro/scRNAseq/reference/GSEA_geneset/", module,"/",list,".txt", sep=""), sep="\t", h=F, skip=2) %>% 
                          as_tibble %>% .$V1 %>% intersect(gene) 

   module_gene <- list(c(module_gene))
   
   # calculation of module score
   data <- AddModuleScore(
            object = data,
            features = module_gene,           
            ctrl = 100,                   
            name = 'module_score')
    
    percentile_99 <- data$module_score1 %>% quantile(0.99); percentile_01 <- data$module_score1 %>% quantile(0.01)
    
    # edit module score for umap visualization
    module_score.editted <- data@meta.data %>% as_tibble %>% select(module_score1) %>%
                               mutate(module_score.editted = if_else(module_score1>percentile_99,percentile_99,
                                                             if_else(module_score1<percentile_01,percentile_01, module_score1)))


    data[["module_score.editted"]] <- module_score.editted$module_score.editted

    
    fp <- FeaturePlot(data,  reduction = "umap", c("module_score.editted")) + 
           scale_colour_viridis_c(option="plasma", limits=c(percentile_01,percentile_99)) +        
           theme(plot.title=element_blank(), title=element_blank(),axis.text=element_blank(),axis.ticks.length = unit(0, "cm"),
                 legend.title = element_blank(),legend.text = element_blank()) + NoLegend()
    ggsave(file = paste0("./Module_score/FP_", list,"_ALL_modified.png"), plot = fp, width = 4.8, height = 4.8, dpi = 300)
    
  }




 