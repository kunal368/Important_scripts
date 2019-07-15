FASTQ1=$1
FASTQ2=${FASTQ1/_read1/_read2}
SAMPLE1=$(basename $FASTQ1 .fastq.gz)
SAMPLE2=$(basename $FASTQ2 .fastq.gz)
DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 20:00
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -M 20000
#BSUB -e /data/zorn_lab_praneet/Kunal_projects/Rafi_mutant_vs_epithelium_mice/Trimmed_fastq/log/%J.err
#BSUB -o /data/zorn_lab_praneet/Kunal_projects/Rafi_mutant_vs_epithelium_mice/Trimmed_fastq/log/%J.out
#BSUB -J $SAMPLE1

cd $DIR
module load bbmap/38.32.0
#module load fastqc
bbduk.sh in1=$FASTQ1 in2=$FASTQ2 out1=./Trimmed_fastq/$SAMPLE1.fastq out2=./Trimmed_fastq/$SAMPLE2.fastq ref=/usr/local/bbmap/38.32.0/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo 
#mlf=50 (minlengthfraction=50) would discard reads under 50% of their original length after trimming

#fastqc -o ./Qc_sample43 E01316_S43_L004_R1_001_clean.fq
EOF