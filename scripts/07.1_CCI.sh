#!/bin/bash
pheno=$1

##--- NATMI ---#
source/work22/home/redahiro/conda_env/NATMI/bin/activate

python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/ExtractEdges.py \
         --species human \
         --emFile ./Data/${pheno}_counts.txt \
         --annFile ./Data/${pheno}_ALL_meta_l3.txt \
         --interDB lrc2p \
         --coreNum 8 \
         --out ${pheno}

python /work22/home/redahiro/software/scRNAseq_script/CCI/NATMI/VisInteractions.py \
       --sourceFolder ${pheno} \
       --interDB lrc2p \
       --expressionThreshold 0 \
       --specificityThreshold 0 \
       --weightType mean --detectionThreshold 0.2 --drawNetwork y \
       --plotWidth 12 --plotHeight 12 --layout circle --fontSize 12 \
       --edgeWidth 6 --maxClusterSize 0 --clusterDistance 1 --plotFormat png

source/work22/home/redahiro/conda_env/NATMI/bin/deactivate

##--- CellPhoneDB ---#
source/work22/home/redahiro/conda_env/CellPhoneDB/bin/activate

cellphonedb method statistical_analysis \
     ./Data/${pheno}_ALL_meta_l3.txt \
     ./Data/${pheno}_counts.txt \
     --output-path=/work22/home/redahiro/analysis/COVID-19_scRNAseq/CCI/CellphoneDB/${pheno} \
     --counts-data=gene_name --iterations=1000 --result-precision 5 --threads=4 --verbose

cellphonedb plot heatmap_plot \
     --pvalues-path=${pheno}/pvalues.txt \
     --output-path=${pheno} \
     ./Data/${pheno}_ALL_meta_l3.txt

