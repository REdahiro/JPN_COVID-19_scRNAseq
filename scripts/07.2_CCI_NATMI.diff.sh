#!/bin/bash

##--- NATMI ---#
source/work22/home/redahiro/conda_env/NATMI/bin/activate

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


