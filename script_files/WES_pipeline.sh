
FASTQ1=$1
FASTQ2=${FASTQ1/_read1/_read2}
ref_fasta=$2
picard_path=/usr/local/picard/2.18.22/picard.jar
Mills_file=/data/zorn_lab_praneet/Kunal_projects/Akaljot_WGS_sequencing/Mills_and_1000G_gold_standard.indels.hg38.vcf
a1000g_file=/data/zorn_lab_praneet/Kunal_projects/Akaljot_WGS_sequencing/1000G_phase1.snps.high_confidence.hg38.vcf
dbSNP_path=/data/zorn_lab_praneet/Kunal_projects/Akaljot_WGS_sequencing/dbsnp_146.hg38.vcf
annovar_humandb_path=/usr/local/annovar/20170616/humandb/


SAMPLE1=$(basename $FASTQ1 .fastq.gz)
SAMPLE2=$(basename $FASTQ2 .fastq.gz)
DIR=$(pwd)
cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 100:00
#BSUB -n 6
#BSUB -R "span[ptile=4]"
#BSUB -M 30000
#BSUB -e /data/zorn_lab_praneet/Kunal_projects/Akaljot_WGS_sequencing/log/%J.err
#BSUB -o /data/zorn_lab_praneet/Kunal_projects/Akaljot_WGS_sequencing/log/%J.out
#BSUB -J $SAMPLE1

cd $DIR
module load bwa/0.7.17
module load gatk/4.1.2.0
module load picard/2.18.22
module load samtools/1.9.0

##Creating UBAM file from fastq files
#java -jar $picard_path FastqToSam F1=$FASTQ1 F2=$FASTQ2 O=unlaigned_read_pairs.bam SM=374.9

#######################################################################################################
			#SamToFastqAndBwaMemAndMba#
######################################################################################################


##Creating BWA Index
#bwa index $ref_fasta

#Overview of command structure
#[SamToFastq] | [BWA-MEM] | [MergeBamAlignment]


#Picard's SamToFastq takes read identifiers, read sequences, and base quality scores to write a Sanger FASTQ format file
#java -Xmx8G -jar $picard_path SamToFastq I=unlaigned_read_pairs.bam FASTQ=unaligned_interleaved.fastq INTERLEAVE=true NON_PF=true

###Align reads with BWA-MEM
#bwa mem -K 100000000 -p -v 3 -t 8 $ref_fasta unaligned_interleaved.fastq > aligned.bam

#MergeBamAlignment
#java -jar $picard_path MergeBamAlignment ALIGNED=aligned.bam UNMAPPED=unlaigned_read_pairs.bam O=merge_alignments.bam R=$ref_fasta SORT_ORDER=unsorted PAIRED_RUN=true PRIMARY_ALIGNMENT_STRATEGY=MostDistant MAX_INSERTIONS_OR_DELETIONS=-1 UNMAP_CONTAMINANT_READS=true VALIDATION_STRINGENCY=SILENT CLIP_ADAPTERS=false


########MarkDuplicates #ASSUME_SORT_ORDER=queryname
#java -jar $picard_path MarkDuplicates I=merge_alignments.bam O=marked_duplicates.bam M=marked_dup_metrics.txt VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname CREATE_MD5_FILE=true

#Sort sam
#java -jar $picard_path SortSam I=marked_duplicates.bam O=/dev/stdout SORT_ORDER=coordinate CREATE_INDEX=false CREATE_MD5_FILE=false | java -jar $picard_path SetNmAndUqTags INPUT=/dev/stdin OUTPUT=sorted.bam CREATE_INDEX=true CREATE_MD5_FILE=true REFERENCE_SEQUENCE=$ref_fasta

#Adding Platform information
#java -jar $picard_path AddOrReplaceReadGroups I=sorted.bam O=sorted_platform.bam SORT_ORDER=coordinate RGPL=illumina RGLB=lib1 RGPU=lane12 RGSM=374.9

##Creating index file for known sites
#gatk IndexFeatureFile -F $a1000g_file
#gatk IndexFeatureFile -F $Mills_file
#gatk IndexFeatureFile -F dbsnp_146.hg38.vcf

###Base Quality Recalibration
#gatk BaseRecalibrator -R $ref_fasta -I sorted_platform.bam -O recal_data.table --known-sites $Mills_file --known-sites $a1000g_file

#gatk ApplyBQSR -R $ref_fasta -I sorted_platform.bam --bqsr-recal-file recal_data.table -O output_bqsr.bam

###hAPLOTYPE CALLER
#gatk HaplotypeCaller -R $ref_fasta -I output_bqsr.bam --dbsnp $dbSNP_path -stand-call-conf 30.0 -O output.vcf


########################### VariantFiteration ##############################

#1. Subset to SNPs-only callset with SelectVariants
#gatk SelectVariants -V output.vcf -select-type SNP -O snps.vcf.gz

#2. Subset to indels-only callset with SelectVariants
#gatk SelectVariants -V output.vcf -select-type INDEL -O indels.vcf.gz

#3. Hard-filter SNPs on multiple expressions using VariantFiltration 
#gatk VariantFiltration \
#-V snps.vcf.gz \
#-filter "QD < 2.0" --filter-name "QD2" \
#-filter "QUAL < 30.0" --filter-name "QUAL30" \
#-filter "SOR > 3.0" --filter-name "SOR3" \
#-filter "FS > 60.0" --filter-name "FS60" \
#-filter "MQ < 40.0" --filter-name "MQ40" \
#-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#-O snps_filtered.vcf.gz

#4. Similarly, hard-filter indels on multiple expressions using VariantFiltration
#gatk VariantFiltration \
#-V indels.vcf.gz \
#-filter "QD < 2.0" --filter-name "QD2" \
#-filter "QUAL < 30.0" --filter-name "QUAL30" \
#-filter "FS > 200.0" --filter-name "FS200" \
#-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
#-O indels_filtered.vcf.gz

#5. Combining filtered variants
#java -jar $picard_path MergeVcfs \
#I=snps_filtered.vcf.gz \
#I=indels_filtered.vcf.gz \
#O=combined_Filtered_Variants.vcf.gz

#Removing reads that fail the filter
#gatk \
#SelectVariants \
#-R $ref_fasta \
#-V combined_Filtered_Variants.vcf.gz \
#-O combined_filtered_removed_variants.vcf.gz \
#--exclude-filtered

#####################Variant Annotation################
table_annovar.pl combined_filtered_removed_variants.vcf $annovar_humandb_path -buildver hg38 -out myanno_1 -remove -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp30a,clinvar_20190513 -operation g,r,f,f,f,f -nastring . -vcfinput


EOF
####