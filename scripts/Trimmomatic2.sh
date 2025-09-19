#!/bin/bash

# adapted from Rachita's trimmomatic script
# new parameters for trimming
# results would be created within the folder where the script is run

data=$1 ## Path to fastq files.

mkdir -p logs

for file in `ls ${data}/*{.,_}R1.{fastq,fq.gz,fq,fastq.gz}`
do
        filename="${file##*/}"      
        #echo $filename
		f2=${file/R1./R2.}
		SN=${filename/.R1.fastq.gz/}
		echo $SN, $f2
        mkdir -p ${SN}
        bsub -q big -sla miket_sc -o logs/report_$filename.out -e logs/report_$filename.err -J Trim_${SN} \
        	"java -jar /apps/lab/miket/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
        	-threads 8 \
    		-phred33 \
    		$file $f2 ${SN}/${SN}.R1.fq.gz ${SN}/${SN}.R1.unpaired.fq.gz ${SN}/${SN}.R2.fq.gz ${SN}/${SN}.R2.unpaired.fq.gz \
    		ILLUMINACLIP:/apps/lab/miket/Trimmomatic-0.36/adapters/TruSeq3-PE-miket.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50"
done
