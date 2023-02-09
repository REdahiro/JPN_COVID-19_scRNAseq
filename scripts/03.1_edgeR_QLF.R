library(edgeR)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
print(args)
clustering = args[1]                       # l1,l2,l3
gene_filter = as.numeric(args[2])          # 1,5,10%


pwd="/work22/home/redahiro/analysis/COVID-19_scRNAseq/"

#---------------- cluster setup -----------------#
cluster_list <- list()
cluster_1 <- read_csv(paste0(pwd,"CellCount/l1.csv")) %>% .$cluster %>% unique
cluster_2 <- read_csv(paste0(pwd,"CellCount/l2.csv")) %>% .$cluster %>% unique
cluster_3 <- read_csv(paste0(pwd,"CellCount/l3.csv")) %>% .$cluster %>% unique

cluster_list[["l1"]] <- cluster_1
cluster_list[["l2"]] <- cluster_2
cluster_list[["l3"]] <- cluster_3

# cluster setting 
cluster_list <- cluster_list[[clustering]]


#---------------- Phenotype_data -----------------#
pheno_base <- read_csv(paste0(pwd,"Pheno_data/Phenotype_analysis_220702.csv"))

for(i in cluster_list){

    #------------------------------#
    #        Data import           #
    #------------------------------#

     # Pseudo_Bulk raw count data
     count_data <- read_delim(paste0(pwd,"edgeR/",clustering,"/Count/PseudoBulk_perSample_counts_",i,".txt"), "\t", col_names=T)
     
     # sample pick up: cell count>5 
     QCed.sample <- read_csv(paste0(pwd,"sceQTL/CellCount/",clustering,".csv")) %>%
                                filter(cluster==i) %>% filter(ncells>5) %>%.$sample 
     
     # filtering count data 
     count_data <- count_data %>% select(Gene,QCed.sample)

     analysis_ID <- colnames(count_data) %>% as_tibble %>% setNames("ID")
     
     pheno <- analysis_ID %>% inner_join(pheno_base, by="ID")

    # Gene filtering
     Case_gene <- read_delim(paste0(pwd,"edgeR/",clustering,"/Gene_Info/Case_", i, "_",gene_filter,".txt"), "\t", col_names=F) %>% .$X1
     Ctrl_gene <- read_delim(paste0(pwd,"edgeR/",clustering,"/Gene_Info/Ctrl_", i, "_",gene_filter,".txt"), "\t", col_names=F) %>% .$X1
     Severe_gene <- read_delim(paste0(pwd,"edgeR/",clustering,"/Gene_Info/Severe_", i, "_",gene_filter,".txt"), "\t", col_names=F) %>% .$X1
     Moderate_gene <- read_delim(paste0(pwd,"edgeR/",clustering,"/Gene_Info/Moderate_", i, "_",gene_filter,".txt"), "\t", col_names=F) %>% .$X1

     Case.Ctrl_gene <- c(Case_gene,Ctrl_gene) %>% unique; length(Case.Ctrl_gene)
     In.Case_gene <- c(Severe_gene,Moderate_gene) %>% unique; length(In.Case_gene)
     Severe.HC_gene <- c(Severe_gene,Ctrl_gene) %>% unique; length(Severe.HC_gene)
     Moderate.HC_gene <- c(Moderate_gene,Ctrl_gene) %>% unique; length(Moderate.HC_gene)
        
    #---------------------------------------------------#
    #                Case vs Control                    #
    #---------------------------------------------------#
    Case.Ctrl_count_data <- count_data %>% filter(Gene %in% Case.Ctrl_gene) %>%
                                  as.data.frame %>% column_to_rownames("Gene")
    
    Status <- pheno$Status 
    Age <- pheno$Age
    Sex <- pheno$Sex

    #----- edgeR -----# 

    dge <- DGEList(Case.Ctrl_count_data, group = Status)          
    dge <- calcNormFactors(dge)
    design <- model.matrix(~Status+Age+Sex)                  
     
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit,coef="Status")        # Case vs Control

    tt <- topTags(qlf, n = Inf)$table %>% rownames_to_column("Gene") %>% as_tibble
    write.csv(tt, paste0(pwd,"edgeR/",clustering,"/Result/Gene_",gene_filter,"/",i,"_Case.Control_adjusted.csv"), row.names=F)


    #-----------------------------------------------------#
    #                 Severe vs Moderate                  #
    #-----------------------------------------------------#

     # Pheno_data
     candidate_ID <- pheno %>% filter(Status ==1) %>% .$ID

     In.Case_count_data <- count_data %>% select(Gene, candidate_ID) %>% filter(Gene %in% all_of(In.Case_gene)) %>%
                                   as.data.frame %>% column_to_rownames("Gene")
     
     candidate_ID <- colnames(In.Case_count_data) %>% as_tibble %>% setNames("ID")

     pheno_SM <- candidate_ID %>% inner_join(pheno_base, by="ID")
    
     Severity <- pheno_SM$Severity
     Status <- pheno_SM$Status 
     Age <- pheno_SM$Age
     Sex <- pheno_SM$Sex
     days_after_onset <- pheno_SM$days_after_onset
     days_after_steroids <- pheno_SM$days_after_steroids 
    
    #--- edgeR ---#
    
    # adjust with age, sex, days since symptom onsets, steroids therapy
    dge <- DGEList(In.Case_count_data, group = Severity)          
    dge <- calcNormFactors(dge)
    design <- model.matrix(~Severity+Age+Sex+days_after_onset+days_after_steroids)
     
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit, coef="Severity")        # Severe vs Moderate
    
    tt <- topTags(qlf, n = Inf)$table %>% rownames_to_column("Gene") %>% as_tibble
    write.csv(tt, paste0(pwd,"edgeR/",clustering,"/Result/Gene_",gene_filter,"/",i,"_Severe.Moderate_adjusted.csv"), row.names=F)


    #-----------------------------------------------------#
    #                 Severe vs HC                        #
    #-----------------------------------------------------#

     # Pheno_data
     candidate_ID <- pheno %>% filter(Severity %in% c(0,2)) %>% .$ID

     S.H_count_data <- count_data %>% select(Gene, candidate_ID) %>% filter(Gene %in% all_of(Severe.HC_gene)) %>%
                                   as.data.frame %>% column_to_rownames("Gene")
     
     candidate_ID <- colnames(S.H_count_data) %>% as_tibble %>% setNames("ID")

     pheno_SH <- candidate_ID %>% inner_join(pheno_base, by="ID")
    
     Severity <- pheno_SH$Severity
     Age <- pheno_SH$Age
     Sex <- pheno_SH$Sex

    
    #--- edgeR ---#
    
    # adjust with age, sex
    dge <- DGEList(S.H_count_data, group = Severity)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~Severity+Age+Sex)
     
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit, coef="Severity")       # Severe vs HC
    
    tt <- topTags(qlf, n = Inf)$table %>% rownames_to_column("Gene") %>% as_tibble
    write.csv(tt, paste0(pwd,"edgeR/",clustering,"/Result/Gene_",gene_filter,"/",i,"_Severe.HC_adjusted.csv"), row.names=F)



    #-----------------------------------------------------#
    #                Moderate vs HC                       #
    #-----------------------------------------------------#

     # Pheno_data
     candidate_ID <- pheno %>% filter(Severity %in% c(0,1)) %>% .$ID

     M.H_count_data <- count_data %>% select(Gene, candidate_ID) %>% filter(Gene %in% all_of(Moderate.HC_gene)) %>%
                                   as.data.frame %>% column_to_rownames("Gene")
     
     candidate_ID <- colnames(M.H_count_data) %>% as_tibble %>% setNames("ID")

     pheno_MH <- candidate_ID %>% inner_join(pheno_base, by="ID") 
    
     Severity <- pheno_MH$Severity
     Age <- pheno_MH$Age
     Sex <- pheno_MH$Sex

    
    #--- edgeR ---#
    
    # adjust with age, sex
    dge <- DGEList(M.H_count_data, group = Severity)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~Severity+Age+Sex)
     
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit, coef="Severity")     # Moderate vs HC
    
    tt <- topTags(qlf, n = Inf)$table %>% rownames_to_column("Gene") %>% as_tibble
    write.csv(tt, paste0(pwd,"edgeR/",clustering,"/Result/Gene_",gene_filter,"/",i,"_Moderate.HC_adjusted.csv"), row.names=F)

   }





