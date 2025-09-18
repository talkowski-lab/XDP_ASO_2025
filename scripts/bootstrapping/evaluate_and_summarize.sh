#!/bin/bash

# Author: Siddharth Reed, Talkowski Lab
# Wrapper around R script to DEG bootstrapping analysis
# Usage 
# $ ./evaluate_and_summarize.sh 

ROOT_DIR="/data/talkowski/Samples/XDP/ASO_dSVA/Results_draft1"
SCRIPT_DIR="${ROOT_DIR}/XDP_ASO_2025/scripts"
PERMUTATION_DIR="${ROOT_DIR}/BootstrapResults"
JOB_DIR="${ROOT_DIR}/Jobs"
METHOD="glm"
CHUNK_SIZE=200  # Number of permutatiosn to compute in parallel before caching

COMPARISONS=(NO.CON.NO.CombinedXDP)

for comparison_dir in "${COMPARISONS[@]}"; do
    comparison_dir="$PERMUTATION_DIR/$comparison_dir"
    comparison=$(basename $comparison_dir)
    for bootstrap_chunk in "$comparison_dir/bootstraps-$comparison".tsv.{000..100}; do
        chunk="$(echo $bootstrap_chunk | rev | cut -d'.' -f1 | rev)"
        finished_sub_chunks=$(find "$comparison_dir" -name "summary-$METHOD-original-$comparison.tsv.$chunk.*" | wc -l)
        # Check if chunk is finshed being computed or areadly submitted and skip
        if [[ -n "$(bjobs -l | grep "$comparison-$chunk")" ]]; then
            echo "$comparison-$chunk Job already submitted" 
        elif [[ $finished_sub_chunks -eq 10 ]]; then
            echo "$comparison-$chunk already finished" 
        else 
            # submit job to run bootstrapping analysis for this chunk of bootstraps
            WORKERS=16
            MEMORY=96000
            bsub -q bigmem \
                 -o "$JOB_DIR/summary-$METHOD-$comparison-$chunk".out \
                 -e "$JOB_DIR/summary-$METHOD-$comparison-$chunk".err \
                 -J "$comparison-$chunk" \
                 -n $WORKERS \
                 -M $MEMORY \
                 -R rusage[mem=$MEMORY] \
            "
            module load R/3.6.3
            Rscript $SCRIPT_DIR/evaluate_and_summarize.R $bootstrap_chunk $METHOD $WORKERS $CHUNK_SIZE
            "
        fi
    done
done
