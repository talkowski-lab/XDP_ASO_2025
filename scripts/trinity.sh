#!/bin/bash

module load /apps/modulefiles/lab/miket/trinity/2.2.0/trinity-2.2.0

bam=$1
name="${bam%%.bam}"
mkdir -p trinity_$name

Trinity \
--SS_lib_type RF \
--min_contig_length 100 \
--genome_guided_max_intron 100000 \
--max_memory 50G \
--CPU 8 \
--genome_guided_bam /data/talkowski/Samples/XDP/Brain_RNASeq/data/Trinity/bamfiles/$bam \
--output trinity_$name

