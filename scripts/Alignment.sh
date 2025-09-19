#!/bin/bash

# adapted from Rachita's Alignment script in scripts_RY, to accommodate directory structure in trimmed
# run in main directory

fdata=$1
prefix=$2
DIR=`pwd`
DATADIR=$DIR/$fdata
WORKDIR=$DIR/data/alignment_${prefix}
mkdir -p $WORKDIR
mkdir -p $WORKDIR/logs
module load star/2.5.3

for file in `ls -1 $DATADIR/*/*.R1{.,_}{fastq,fq.gz,fq,fastq.gz}`
do

ID=`basename "$file"`
ID=${ID/.R1.fq.gz/}
mkdir -p $WORKDIR/$ID

echo $DIR
echo $ID

bsub -sla miket_sc -q big-multi -n 8 -M 50000 -J ${ID} \
      -o $WORKDIR/logs/${ID}.out -e $WORKDIR/logs/${ID}.err \
      "STAR --runThreadN 8 \
        --genomeDir /data/talkowski/tools/ref/RNA-Seq/human/star_2.4.2a_ref/star_GRCh37.75_sva_corrected_ercc \
        --twopassMode Basic \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNoverLmax 0.05 \
        --outSAMtype BAM Unsorted \
        --outReadsUnmapped Fastx \
        --readFilesCommand zcat \
        --alignIntronMin 21 \
        --alignIntronMax 0 \
        --quantMode GeneCounts \
        --alignEndsType EndToEnd \
        --outFileNamePrefix $WORKDIR/${ID}/${ID}. \
        --readFilesIn $DATADIR/${ID}/${ID}.R1.fq.gz $DATADIR/${ID}/${ID}.R2.fq.gz; \
	samtools sort -o $WORKDIR/${ID}/${ID}.Aligned.sortedByCoord.out.bam $WORKDIR/${ID}/${ID}*.Aligned.out.bam; \
      samtools index $WORKDIR/${ID}/${ID}.Aligned.sortedByCoord.out.bam;"
done
