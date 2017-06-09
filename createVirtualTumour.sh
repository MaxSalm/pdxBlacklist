#!/bin/bash

## Script to create "synthetic" tumour data using high-coverage 1000 Genome Project Data

## Global variables
PARENT="/farm/babs/data/genomes/homo_sapiens/1000GenomeProject/PlatinumGenomes/Trio/NA12878.bam"
CHILD="/farm/babs/data/genomes/homo_sapiens/1000GenomeProject/PlatinumGenomes/Trio/NA12882.bam"
DIR_sh="/farm/scratch/rs-bio-lif/salm01/shell_script_graveyard/"
DIR_out="/farm/babs/data/genomes/homo_sapiens/1000GenomeProject/PlatinumGenomes/Trio/VirtualTumor/"


## Programs etc.
JAVA="/farm/babs/bin/java1.6"
Xmx=" -Xmx50g"
Xms=" -Xms8g"
PICARD=$JAVA$Xmx$Xms" -jar /farm/home/salm01/NGS_gen/picard-tools-1.59/"
DSF=$PICARD"DownsampleSam.jar"
MSF=$PICARD"MergeSamFiles.jar"
options=" SORT_ORDER=coordinate CREATE_INDEX=true MAX_RECORDS_IN_RAM=2000000 TMP_DIR=/farm/home/salm01/tmp VALIDATION_STRINGENCY=LENIENT "
options0=" CREATE_INDEX=true MAX_RECORDS_IN_RAM=2000000 TMP_DIR=/farm/home/salm01/tmp VALIDATION_STRINGENCY=LENIENT "

VARSCAN="/farm/babs/redhat6/bin/java -Xmx20g -jar /farm/home/salm01/bin/VarScan.v2.3.3.jar somatic"

MPILEUP="/farm/babs/redhat6/bin/samtools mpileup -A -B -C 50 -d 10000 -m 3 -F 0.0002 -q 20 -Q 20 -f /farm/babs/data/genomes/homo_sapiens/GATK_bundle_2.2/hg19/ucsc.hg19.fasta"


## Start
## Split files into 'tumor' and 'germline' (restrict to chr17 (NB UCSC encoding & contains BRCA1) to speed-up analyses)
## ~200x each, so should end up ~100x
## NB have to use PICARD, as Samtools won't do this properly
## downsample BAMs
cd $DIR_out 
samtools view -b $PARENT chr17 chr18 chr19 | $DSF INPUT=/dev/stdin OUTPUT=NA12878_fakeTumour.bam RANDOM_SEED=123 PROBABILITY=0.5$options0	
samtools view -b $PARENT chr17 chr18 chr19 | $DSF INPUT=/dev/stdin OUTPUT=NA12878_fakeNormal.bam RANDOM_SEED=456 PROBABILITY=0.5$options0
#samtools view -b $CHILD chr17 chr18 chr19 | $DSF INPUT=/dev/stdin OUTPUT=NA12882_fakeTumour.bam PROBABILITY=0.5$options0

## Merge BAMs, mix mother and daughter (50/50) to simulate tumour (won't work as they have different sequence dictionaries) 
# samtools merge -r -f synthetic.bam NA12878_fakeTumour.bam NA12882_fakeTumour.bam 	
	
## Run VarScan (note the removal of zero-coverage events, and relax calling params to get more false positives)
$MPILEUP NA12878_fakeTumour.bam NA12878_fakeNormal.bam | awk '{ if ($4 > 1 && $7 > 1) print }' > falsePositive.mpileup
$VARSCAN falsePositive.mpileup falsePositive --mpileup 1 --min-coverage 8 --min-coverage-normal 10 --min-coverage-tumor 6 --min-var-freq 0.01 --min-freq-for-hom 0.75 --normal-purity 0.9 --p-value 0.99 --somatic-p-value 0.1 --tumor-purity 0.5 --strand-filter 0

grep -i Somatic falsePositive.snp | awk '{if ($15 < 0.05) print}' > falsePositive_somatic2.snp
grep -i Somatic falsePositive.indel | awk '{if ($15 < 0.05) print}' > falsePositive_somatic2.indel

## True positives (mum ('tumor') specific SNVs, not shared with daughter ('germline') - models 100% tumour)
#$MPILEUP NA12882_fakeTumour.bam NA12878_fakeTumour.bam | awk '{ if ($4 > 9 && $7 > 5) print }' > TruePositive.mpileup
#$VARSCAN TruePositive.mpileup truePositive --mpileup 1 --min-coverage 8 --min-coverage-normal 10 --min-coverage-tumor 6 --min-var-freq 0.01 --min-freq-for-hom 0.75 --normal-purity 0.9 --p-value 0.99 --somatic-p-value 0.05 --tumor-purity 0.5 --strand-filter 0





## High-quality SNV/Indel for NA12878 in GenomeInABottle (http://www.nature.com/nbt/journal/v32/n3/full/nbt.2835.html)






