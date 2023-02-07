#!/bin/bash

scdrs munge-gs \
    --out-file munge.gs \
    --zscore-file zscore.txt \
    --weight zscore \
    --n-max 100

mkdir -p output

scdrs compute-score \
    --h5ad-file covid-19.h5ad\
    --h5ad-species human\
    --gs-file munge.gs\
    --gs-species human\
    --out-folder output\
    --cov-file covariate.txt\
    --flag-filter-data True\
    --flag-raw-count True\
    --n-ctrl 1000\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl-norm-score True

scdrs perform-downstream \
    --h5ad-file covid-19.h5ad\
    --score-file output/Trait.full_score.gz\
    --out-folder output\
    --group-analysis l1\
    --gene-analysis\
    --flag-filter-data True\
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file covid-19.h5ad\
    --score-file output/Trait.full_score.gz\
    --out-folder output\
    --group-analysis l2\
    --gene-analysis\
    --flag-filter-data True\
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file covid-19.h5ad\
    --score-file output/Trait.full_score.gz\
    --out-folder output\
    --group-analysis l3\
    --gene-analysis\
    --flag-filter-data True\
    --flag-raw-count True




