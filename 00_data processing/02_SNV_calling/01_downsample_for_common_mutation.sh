#
# downsample chrM bam file for pseudo-bulk bam:
# run the strandard pepiline for finding the positive control(ref: scATAC_mtSMC  pipeline)
# Filter mapping results with q 30, PE 
# 

# Step1: get the chrM bam file and downsampling;
AR3_C4_bam=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/outs/possorted_bam.bam
nohup samtools view -b $AR3_C4_bam chrM  > AR3_C4_100bp_possorted_chrM.bam &
nohup samtools index AR3_C4_100bp_possorted_chrM.bam &

AR3_C5_bam=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/outs/possorted_bam.bam
nohup samtools view -b $AR3_C5_bam chrM  > AR3_C5_100bp_possorted_chrM.bam &
nohup samtools index AR3_C5_100bp_possorted_chrM.bam &

AR3_C4_chrM_bam=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/common_mutation/AR3_C4_100bp_possorted_chrM.bam
AR3_C5_chrM_bam=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/common_mutation/AR3_C5_100bp_possorted_chrM.bam
nohup samtools view -@ 10 -s 0.3 -b $AR3_C4_chrM_bam > AR3_C4_chrM_downsampled.bam &
nohup samtools view -@ 10 -s 0.3 -b $AR3_C5_chrM_bam > AR3_C5_chrM_downsampled.bam &

# Step2: bwa realign the bam file：
AR3_C4_chrM_downsampled=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/common_mutation/AR3_C4_chrM_downsampled.bam
AR3_C5_chrM_downsampled=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/common_mutation/AR3_C5_chrM_downsampled.bam
bwa index assembly_hardmask.fa

cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/common_mutation
samtools collate -Oun128 $AR3_C4_chrM_downsampled | samtools fastq -OT RG,BC - \
  | bwa mem -pt8 -CH <(samtools view -H $AR3_C4_chrM_downsampled|grep ^@RG) /md01/nieyg/ref/hard-mask/bwa_index/assembly_hardmask.fa - \
  | samtools sort -@4 -m4g -o AR3_C4_chrM_downsampled_bwa.bam -

cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/common_mutation
samtools collate -Oun128 $AR3_C5_chrM_downsampled | samtools fastq -OT RG,BC - \
  | bwa mem -pt8 -CH <(samtools view -H $AR3_C5_chrM_downsampled|grep ^@RG) /md01/nieyg/ref/hard-mask/bwa_index/assembly_hardmask.fa - \
  | samtools sort -@4 -m4g -o AR3_C5_chrM_downsampled_bwa.bam -

# Step3: picard remove the replicates
conda activate MitoSort
picard MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
INPUT=AR3_C4_chrM_downsampled_bwa.bam \
OUTPUT=AR3_C4_chrM_downsampled_bwa_rmdup.bam  \
METRICS_FILE=AR3_C4_chrM_downsampled_bwa.rmdup.metrics

picard MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
INPUT=AR3_C5_chrM_downsampled_bwa.bam \
OUTPUT=AR3_C5_chrM_downsampled_bwa_rmdup.bam  \
METRICS_FILE=AR3_C5_chrM_downsampled_bwa.rmdup.metrics

conda deactivate
# Step4: Filter genome reads, q >30,PE for SNP calling  
samtools view -b -q  30 -f 0x2 AR3_C4_chrM_downsampled_bwa_rmdup.bam -o AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam chrM
samtools view -b -q  30 -f 0x2 AR3_C5_chrM_downsampled_bwa_rmdup.bam -o AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam chrM
samtools index AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam
# mtQ30PE_rate:

# Step5: make re-align target regions
java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T RealignerTargetCreator  -nt 32 \
-I AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam \
-o AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam.target.intervals

java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam \
-targetIntervals AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam.target.intervals \
-o AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.realign.bam

java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T RealignerTargetCreator  -nt 32 \
-I AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam \
-o AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam.target.intervals

java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam \
-targetIntervals AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.bam.target.intervals \
-o AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.realign.bam

# Step6: mpileup get file for snv calling

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
-f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa \
-q 10 -Q 10 AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.realign.bam > AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.mpileup

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
-f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa \
-q 10 -Q 10 AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.realign.bam > AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.mpileup

# make count table for each allele 
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.mpileup > AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.mpileup.count
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.mpileup > AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.mpileup.count

# Step5: varscan call somatic mutation in chrM
varscan  pileup2snp AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE.mpileup --min-coverage 3 --min-var-freq 0.0000001  --min-reads2 2 > AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE_filtered.snv
varscan  pileup2snp AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE.mpileup --min-coverage 3 --min-var-freq 0.0000001  --min-reads2 2 > AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE_filtered.snv

# filter somatic mutation  --step6
# STEP 1: filter mtSNV
# $17>=1 && $18>=1 ( make sure both positive and negative strand were called 
# Minor  Allele frequency more than 0.01 
# remove G->T and C->A, "N",based on previous empirical evidence. 
# sequence depth >= 20 for each cell, the sequence depth criteria should be adjust according to the data. 

awk '$17>=1 && $18>=1' AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE_filtered.snv | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'   | awk ' $6+$5>=20' >>AR3_C4_chrM_downsampled_bwa_rmdup_genomeQ30PE_filtered.mtDNA.snv.filter
awk '$17>=1 && $18>=1' AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE_filtered.snv | sed '1,1d' |awk '$6/($6+$5)>=0.10' |awk '!($3=="G" && $6=="T")' |awk '!($3=="C" && $6=="A")' |awk '$3!="N"'   | awk ' $6+$5>=20' >>AR3_C5_chrM_downsampled_bwa_rmdup_genomeQ30PE_filtered.mtDNA.snv.filter


