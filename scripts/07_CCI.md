

##-----------------------------------------------##
##      step1: Making count_data & meta_data     ##
##-----------------------------------------------##
#--- DownSampling dataset ---#
# 4645912 (67G)
qsub -q large@node101 -l m_mem_free=20G -l s_vmem=20G -pe OpenMP 8 -N CCI -cwd -o step1.o -e step1.e -V -b y \
"bash /work22/home/redahiro/analysis/COVID-19_scRNAseq/Script/CCI/01_Making.Matrix/CCI_prepare.sh"


#--- CXCL10 & CXCR3, TNFSF10(=TRAIL) & TNFRSF10B ---#
Rscript CXCL10_ds_rds.file.R
python3 /work22/home/redahiro/analysis/COVID-19_scRNAseq/Script/CCI/01_Making.Matrix/02_sparse_to_txt.py Ctrl
python3 /work22/home/redahiro/analysis/COVID-19_scRNAseq/Script/CCI/01_Making.Matrix/02_sparse_to_txt.py Severe
python3 /work22/home/redahiro/analysis/COVID-19_scRNAseq/Script/CCI/01_Making.Matrix/02_sparse_to_txt.py Moderate





##--------------------------------##
##            NATMI               ##
##--------------------------------##
# dir: /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/NATMI
## "lrdbs" has to be placed on the directory of the analysis

######------- STEP1: Extract ligand-receptor-mediated interactions -----#####

# Case, Ctrl: 4683408(620G), 4684206(642G)
for pheno in Case Ctrl; do
qsub -q large@node101 -l m_mem_free=100G -l s_vmem=100G -pe OpenMP 8 -N ${pheno}_NATMI -cwd -o ${pheno}.o -e ${pheno}.e -V -b y \
"python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/ExtractEdges.py \
         --species human \
         --emFile /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/${pheno}_counts.txt \
         --annFile /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/${pheno}_ALL_meta_l3.txt \
         --interDB lrc2p \
         --coreNum 8 \
         --out ${pheno}";
done


# Severe: 4683410
qsub -q large@node101 -l m_mem_free=60G -l s_vmem=60G -pe OpenMP 8 -N Severe_NATMI -cwd -o Severe.o -e Severe.e -V -b y \
"python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/ExtractEdges.py \
         --species human \
         --emFile /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Severe_counts.txt \
         --annFile /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Severe_ALL_meta_l3.txt \
         --interDB lrc2p \
         --coreNum 8 \
         --out Severe"

# Moderate: 4683411
qsub -q large@node101 -l m_mem_free=60G -l s_vmem=60G -pe OpenMP 1 -N Moderate_NATMI -cwd -o Moderate.o -e Moderate.e -V -b y \
"python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/ExtractEdges.py \
         --species human \
         --emFile /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Moderate_counts.txt \
         --annFile /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Moderate_ALL_meta_l3.txt \
         --interDB lrc2p \
         --coreNum 1 \
         --out Moderate"




######------- STEP2: Visualise ligand-receptor-mediated interaction network -----#####
# Filtering: --expressionThreshold  or --specificityThreshold 0.1 etc (Default: 0) 

# Default setting
python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/VisInteractions.py \
       --sourceFolder Ctrl \
       --interDB lrc2p \
       --expressionThreshold 0 \
       --specificityThreshold 0 \
       --weightType mean --detectionThreshold 0.2 --drawNetwork y \
       --plotWidth 12 --plotHeight 12 --layout circle --fontSize 12 \
       --edgeWidth 6 --maxClusterSize 0 --clusterDistance 1 --plotFormat png



##-----------------------------------------------##
##         Changes  between Condition            ##
##-----------------------------------------------##

##--- Identification of changes in L-R edge between in two-conditions ---##

# Case vs Ctrl
python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/DiffEdges.py \
       --refFolder Ctrl --targetFolder Case \
       --interDB lrc2p --out Status

python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/VisInteractions.py \
       --sourceFolder Status --interDB lrc2p \
       --specificityThreshold 0 \
       --weightType mean --detectionThreshold 0.2 --drawNetwork y \
       --plotWidth 12 --plotHeight 12 --layout circle --fontSize 12 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 1 \
       --plotFormat png

# Severe vs Moderate
python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/DiffEdges.py \
       --refFolder Moderate --targetFolder Severe \
       --interDB lrc2p --out Severity


python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/VisInteractions.py \
       --sourceFolder Severity --interDB lrc2p \
       --specificityThreshold 0 \
       --weightType mean --detectionThreshold 0.2 --drawNetwork y \
       --plotWidth 12 --plotHeight 12 --layout circle --fontSize 12 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 1 \
       --plotFormat png







##------ CXCL10 -----#
# dir: /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/NATMI/CXCL10

######------- STEP1: Extract ligand-receptor-mediated interactions -----#####

for Status in Ctrl Severe Moderate; do
qsub -q large@node101 -l m_mem_free=20G -l s_vmem=20G -pe OpenMP 8 -N ${Status}_NATMI -cwd -o ${Status}.o -e ${Status}.e -V -b y \
"python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/ExtractEdges.py \
         --species human \
         --emFile /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/CXCL10/${Status}_counts.txt \
         --annFile /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/CXCL10/${Status}_ALL_meta_l3.txt \
         --interDB lrc2p \
         --coreNum 8 \
         --out ${Status}";
done




######------- STEP2: Visualise ligand-receptor-mediated interaction network -----#####
# Filtering: --expressionThreshold   or --specificityThreshold 0.1 

python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/VisInteractions.py \
       --sourceFolder Severe \
       --interDB lrc2p \
       --expressionThreshold 0 \
       --specificityThreshold 0 \
       --weightType mean --detectionThreshold 0 --drawNetwork y \
       --plotWidth 12 --plotHeight 12 --layout circle --fontSize 12 \
       --edgeWidth 6 --maxClusterSize 0 --clusterDistance 1 --plotFormat png












##-----------------------------##
##        CellphoneDB          ##
##-----------------------------##
#--- statistical analysis (Directory automatically created)---#
# Case: 4737394          4343674(280G)
# Control: 4736873 (440G)       4343675(280G)

qsub -q large@node101 -l m_mem_free=140G -l s_vmem=140G -pe OpenMP 4 -N Case_CellPhone -cwd -o Case.o -e Case.e -V -b y \
"cellphonedb method statistical_analysis \
/work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Case_ALL_meta_l3.txt \
/work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Case_counts.txt \
--output-path=/work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/CellphoneDB/Case \
--counts-data=gene_name --iterations=1000 --result-precision 5 --threads=4 --verbose"


qsub -q large@node101 -l m_mem_free=120G -l s_vmem=120G -pe OpenMP 4 -N Ctrl_CellPhone -cwd -o Ctrl.o -e Ctrl.e -V -b y \
"cellphonedb method statistical_analysis \
/work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Ctrl_ALL_meta_l3.txt \
/work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Ctrl_counts.txt \
--output-path=/work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/CellphoneDB/Ctrl \
--counts-data=gene_name --iterations=1000 --result-precision 5 --threads=4 --verbose"

#--- heat_map ---#
cellphonedb plot heatmap_plot \
     --pvalues-path=./Case/pvalues.txt \
     --output-path=./Case \
     /work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/Data/Case_ALL_meta_l3.txt


#--- heat_map modified (range: manually edit)---#   
Rscript Heatmap_modified.R \
../Data/Case_ALL_meta_l3.txt \
./Case/pvalues.txt \
./annotation_list/l3.txt \
./Case/heatmap_for_paper.pdf



# args
meta_data = args[1]
pvalues_data = args[2]
annotation_list = args[3]
graph_name = args[4]


meta_data ="../Data/DS/Moderate_ALL_meta_l3_v2.txt"
pvalues_data ="./Moderate_DS_l3_v2/pvalues.txt" 
annotation_list ="./annotation_list/l3_v2.txt"