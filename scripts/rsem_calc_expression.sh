# run in directory where fastqs as loop through fastq1s

module load bowtie2
module load /apps/modulefiles/lab/miket/rsem/1.2.31/rsem-1.2.31

WorkDir=`pwd`

for fastq1 in `ls *_trimmed_1.fq`; do
	
	fastq2=${fastq1/_1.fq/_2.fq};
	name="${fastq1%%_trimmed_1.fq}";
	echo $name;

	bsub -q normal -sla miket_sc -J rsem_$name -R 'hname!=cmu085' -R 'hname!=cn078' -R 'hname!=cn063' -o logs/$name.out -e logs/$name.err "rsem-calculate-expression \
  		--forward-prob=0 \
  		--bowtie2 \
  		--bowtie2-mismatch-rate 0.05 \
  		--estimate-rspd \
  		--paired-end $fastq1 $fastq2 $WorkDir/TAF1_contigs.fa $name";

done

