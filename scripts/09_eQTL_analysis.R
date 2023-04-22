library(tidyverse)

PWD="/work22/home/redahiro/analysis/COVID-19_scRNAseq/sceQTL/"

args <- commandArgs(trailingOnly = TRUE)
print(args)
clustering = args[1]

# gene setup
gene_list <- read_csv(paste0(PWD,"/Genome/VCF/round5_CHR.POS_V2G.csv")) %>% .$V2G

# cluster setup
cluster_list <- list()
cluster_1 <- c("CD4T","CD8T","NK","Mono","DC","B")
cluster_2 <- c("CD14p_Mono","CD16p_Mono","cDC","pDC")
cluster_3 <- c("cMono_S100A","cMono_IL1B","cMono_CCL3","intMono","ncMono","cDC","pDC")

cluster_list[["l1"]] <- cluster_1
cluster_list[["l2"]] <- cluster_2
cluster_list[["l3"]] <- cluster_3

#---------------------------#
#           Genome          #
#---------------------------#
# Dosage
dosage <- read_csv(paste0(PWD,"/Genome/Dosage/round5_genotype.dosage.csv"))

# Genome_PCA
pc_dna <- read.table(paste0(PWD,"/Genome/PLINK/PCA_COV70_HC75.evec"), sep="", h=F, skip=1) %>% 
             as_tibble %>% select(V1:V6) %>% setNames(c("ID","PC1_d","PC2_d","PC3_d","PC4_d","PC5_d"))

genome <- dosage %>% left_join(pc_dna, by ="ID") %>% mutate(IID=ID) %>% select(-ID) %>% select(IID,everything())

# SGID.CTID
sample <- read_csv(paste0(PWD,"/Genome/ID/COVID_SG_TaskForce_ID_94samples.csv"))

genome <- genome %>% left_join(sample,by="IID") %>%
               mutate(Status = str_sub(IID, start = 1, end = 2),
                      ID=if_else(Status=="SG",IID,ID)) %>% 
               select(-IID,-Status) %>% select(ID,everything())


# RNA_PCA
pc_rna <- read_delim("./PCA/PCA_RNA_Pseudobulk_ALL.txt", "\t", col_names=T)

# Pheno
pheno <- read_csv("./Phenotype/Phenotype_data_220617.csv") %>% 
             mutate(Status = if_else(Severity == "HC", 0, 1),
                    Severity = if_else(Severity == "Severe", 2,
                               if_else(Severity == "Moderate", 1, 0)),
                    Sex = if_else(Sex == "M", 1, 2)) %>%
              select(ID,Age,Sex,Status,Severity,days_after_onset,days_after_steroids)

# ALL: 143 samples
pheno <- pheno %>% 
            inner_join(genome, by="ID") %>%
            inner_join(pc_rna, by="ID")

#-------------------------------------#
#            eQTL analysis            #
#-------------------------------------#

rm(result)
result <- list()

   count = read_csv(paste0("./CellCount/",clustering,".csv"))
   
   cluster_list <- cluster_list[[clustering]]
   #cluster_list <- count$cluster %>% unique
   
# eQTL analysis per cluster

for(i in cluster_list){

     print(paste0("START ",i))

     # PseudoBulk_Matrix import
     gene_matrix <- read_delim(paste0("./Matrix_",clustering,"/PseudoBulk_matrix_", i,".txt"))

     # Sample filtering
     sub_sample <- count %>% filter(cluster==i) %>% filter(ncells>=6) %>%.$sample
     
     # filtering of samples & genes
     gene_analysis <- intersect(gene_list,colnames(gene_matrix))
     gene_matrix <- gene_matrix %>% filter(ID %in% sub_sample) %>% select(ID,gene_analysis)
     
     ##--------------------------------------------##
     ##     merge expression and phenotype data    ##
     ##--------------------------------------------##

     analysis <- gene_matrix %>% inner_join(pheno,by="ID")


   # eQTL analysis per Gene
   for(g in gene_analysis){
         print(paste0("START ",g))

         gene.d <- paste0(g,".d")

         sub <- analysis %>% select(ID,g, gene.d, Age, Sex, Status, Severity, days_after_onset, days_after_steroids, PC1:PC5, PC1_d:PC5_d)
         col.names_base <- colnames(sub)[4:19]
         colnames(sub) <- c("ID","eGene","dosage",col.names_base)

         #----- making subset -----#
         sub_case <- sub %>% filter(Status=="1")
         sub_hc <- sub %>% filter(Status=="0")

        # Case: eQTL anlaysis
         nocov <- lm(eGene ~ dosage, data=sub_case)
         Beta_no = summary(nocov)$coef[2,1]
         SE_no = summary(nocov)$coef[2,2]
         P_no = summary(nocov)$coef[2,4]

         PC2 <- lm(eGene ~ dosage + Age + Sex + Severity + days_after_onset + days_after_steroids + PC1 + PC2 + PC1_d + PC2_d, data=sub_case)
         Beta_PC2 = summary(PC2)$coef[2,1]
         SE_PC2 = summary(PC2)$coef[2,2]
         P_PC2 = summary(PC2)$coef[2,4]

         sub_result = c(i,g,"Case_ALL",Beta_no,SE_no,P_no,Beta_PC2,SE_PC2,P_PC2) %>% as_tibble %>% t() %>% as.matrix
         colnames(sub_result) <- c("Cluster","Gene","Condition","B_no","SE_no","P_no","B_PC2","SE_PC2","P_PC2")
        
         result[[paste0("Case_ALL_",i,"_",g)]] <- sub_result
        
        # HC: eQTL anlaysis
         nocov <- lm(eGene ~ dosage, data=sub_hc)
         Beta_no = summary(nocov)$coef[2,1]
         SE_no = summary(nocov)$coef[2,2]
         P_no = summary(nocov)$coef[2,4]

         PC2 <- lm(eGene ~ dosage + Age + Sex + PC1 + PC2 + PC1_d + PC2_d, data=sub_hc)
         Beta_PC2 = summary(PC2)$coef[2,1]
         SE_PC2 = summary(PC2)$coef[2,2]
         P_PC2 = summary(PC2)$coef[2,4]
   
         sub_result = c(i,g,"HC_ALL",Beta_no,SE_no,P_no,Beta_PC2,SE_PC2,P_PC2,Beta_PC2,SE_PC2,P_PC2) %>% as_tibble %>% t() %>% as.matrix
         colnames(sub_result) <- c("Cluster","Gene","Condition","B_no","SE_no","P_no","B_PC2","SE_PC2","P_PC2")
        
         result[[paste0("HC_ALL_",i,"_",g)]] <- sub_result


        # Interaction (Status * dosage)
         nocov <- lm(eGene ~ dosage + Status + dosage*Status, data=sub)
         Beta_no = summary(nocov)$coef[4,1]
         SE_no = summary(nocov)$coef[4,2]
         P_no = summary(nocov)$coef[4,4]

         PC2 <- lm(eGene ~ dosage + Status + dosage*Status + Age + Sex + PC1 + PC2 + PC1_d + PC2_d, data=sub)
         Beta_PC2 = summary(PC2)$coef[10,1]
         SE_PC2 = summary(PC2)$coef[10,2]
         P_PC2 = summary(PC2)$coef[10,4]

         sub_result = c(i,g,"Interaction",Beta_no,SE_no,P_no,Beta_PC2,SE_PC2,P_PC2) %>% as_tibble %>% t() %>% as.matrix
         colnames(sub_result) <- c("Cluster","Gene","Condition","B_no","SE_no","P_no","B_PC2","SE_PC2","P_PC2")
        
         result[[paste0("Interaction_",i,"_",g)]] <- sub_result

         }
     
     }


   summary <- do.call(rbind, result) %>% as_tibble %>% mutate_at(vars("B_no","SE_no","P_no","B_PC2","SE_PC2","P_PC2"),as.numeric)

   write.csv(summary, paste0("eQTL_round5_variant_results_",clustering,"_PCrna_ALL.csv"), row.names=F)


