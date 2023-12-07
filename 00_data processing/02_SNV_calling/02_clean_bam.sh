# Step1: bwa realign the chrM bam file 
AR3_C4_bam=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/outs/possorted_bam.bam
nohup samtools view -b $AR3_C4_bam chrM  > AR3_C4_100bp_possorted_chrM.bam &
nohup samtools index AR3_C4_100bp_possorted_chrM.bam &

AR3_C5_bam=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/outs/possorted_bam.bam
nohup samtools view -b $AR3_C5_bam chrM  > AR3_C5_100bp_possorted_chrM.bam &
nohup samtools index AR3_C5_100bp_possorted_chrM.bam &


AR3_C4_chrM_bam=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/01_clean_bam/AR3_C4_100bp_possorted_chrM.bam
AR3_C5_chrM_bam=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/01_clean_bam/AR3_C5_100bp_possorted_chrM.bam

cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/01_clean_bam
nohup samtools collate -Oun128 $AR3_C4_chrM_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $AR3_C4_chrM_bam|grep ^@RG) /md01/nieyg/ref/hard-mask/bwa_index/assembly_hardmask.fa - \
  | samtools sort -@4 -m4g -o AR3_C4_chrM_bwa2.bam - &

cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/01_clean_bam
nohup samtools collate -Oun128 $AR3_C5_chrM_bam | samtools fastq -OT RG,CB,CR,CY,TR,TQ - \
  | bwa mem -pt8 -CH <(samtools view -H $AR3_C5_chrM_bam|grep ^@RG) /md01/nieyg/ref/hard-mask/bwa_index/assembly_hardmask.fa - \
  | samtools sort -@4 -m4g -o AR3_C5_chrM_bwa2.bam - &

# Step2: picard remove the replicates
conda activate MitoSort
nohup picard -Xmx100g MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
INPUT=AR3_C4_chrM_bwa.bam \
OUTPUT=AR3_C4_chrM_bwa_rmdup.bam  \
METRICS_FILE=AR3_C4_chrM_bwa.rmdup.metrics &

nohup picard -Xmx100g MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true \
INPUT=AR3_C5_chrM_bwa.bam \
OUTPUT=AR3_C5_chrM_bwa_rmdup.bam  \
METRICS_FILE=AR3_C5_chrM_bwa.rmdup.metrics &
conda deactivate

# Step3: Filter chrM reads, q >60,PE for SNP calling  
samtools view -b -q  60 -f 0x2 AR3_C4_chrM_bwa_rmdup.bam -o AR3_C4_chrM_bwa_rmdup_chrMQ60PE.bam chrM
samtools view -b -q  60 -f 0x2 AR3_C5_chrM_bwa_rmdup.bam -o AR3_C5_chrM_bwa_rmdup_chrMQ60PE.bam chrM
samtools index AR3_C4_chrM_bwa_rmdup_chrMQ60PE.bam
samtools index AR3_C5_chrM_bwa_rmdup_chrMQ60PE.bam

# remove reads mapped with different chrom
samtools view -H AR3_C4_chrM_bwa_rmdup_chrMQ60PE.bam > header.sam
samtools view  -@ 12 AR3_C4_chrM_bwa_rmdup_chrMQ60PE.bam |awk '{if ($7 == "=") print $0}' | grep -v "MM:i:1" | cat header.sam - | samtools view -Sb -@ 12  - > AR3_C4_chrM_bwa_rmdup_chrMQ60PE.unique.bam
samtools index AR3_C4_chrM_bwa_rmdup_chrMQ60PE.unique.bam
samtools sort -@ 12 AR3_C4_chrM_bwa_rmdup_chrMQ60PE.unique.bam -o AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.bam
samtools index AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.bam

# remove reads mapped with different chrom
samtools view -H AR3_C5_chrM_bwa_rmdup_chrMQ60PE.bam > header.sam
samtools view  -@ 12 AR3_C5_chrM_bwa_rmdup_chrMQ60PE.bam |awk '{if ($7 == "=") print $0}' | grep -v "MM:i:1" | cat header.sam - | samtools view -Sb -@ 12  - > AR3_C5_chrM_bwa_rmdup_chrMQ60PE.unique.bam
samtools index AR3_C5_chrM_bwa_rmdup_chrMQ60PE.unique.bam
samtools sort -@ 12 AR3_C5_chrM_bwa_rmdup_chrMQ60PE.unique.bam -o AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.bam
samtools index AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.bam


# Step4: make re-align target regions
java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T RealignerTargetCreator  -nt 32 \
-I AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.bam \
-o AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.bam.target.intervals

java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.bam \
-targetIntervals AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.bam.target.intervals \
-o AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.realign.bam


java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T RealignerTargetCreator  -nt 32 \
-I AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.bam \
-o AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.bam.target.intervals


java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.bam \
-targetIntervals AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.bam.target.intervals \
-o AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.realign.bam

samtools view -h AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.realign.bam | awk '{if ($0 !~ /NM:i:[2-9]/ && $0 !~ /NM:i:[0-9][0-9]/) print $0}' | samtools view -b -o AR3_C4_chrM_bwa_rmdup_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam -

samtools view -h AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.realign.bam | awk '{if ($0 !~ /NM:i:[2-9]/ && $0 !~ /NM:i:[0-9][0-9]/) print $0}' | samtools view -b -o AR3_C5_chrM_bwa_rmdup_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam -

###### not rmdup ######

# Step3: Filter chrM reads, q >60,PE for SNP calling  
# remove reads mapped with different chrom

samtools view -b -q  60 -f 0x2 ../AR3_C4_chrM_bwa.bam -o AR3_C4_chrM_bwa_chrMQ60PE.bam chrM
samtools index AR3_C4_chrM_bwa_chrMQ60PE.bam
samtools view -H AR3_C4_chrM_bwa_chrMQ60PE.bam > header.sam
samtools view  -@ 12 AR3_C4_chrM_bwa_chrMQ60PE.bam |awk '{if ($7 == "=") print $0}' | grep -v "MM:i:1" | cat header.sam - | samtools view -Sb -@ 12  - > AR3_C4_chrM_bwa_chrMQ60PE.unique.bam
samtools index AR3_C4_chrM_bwa_chrMQ60PE.unique.bam
samtools sort -@ 12 AR3_C4_chrM_bwa_chrMQ60PE.unique.bam -o AR3_C4_chrM_bwa_chrMQ60PE_filtered.bam
samtools index AR3_C4_chrM_bwa_chrMQ60PE_filtered.bam


samtools view -b -q  60 -f 0x2 ../AR3_C5_chrM_bwa.bam -o AR3_C5_chrM_bwa_chrMQ60PE.bam chrM
samtools index AR3_C5_chrM_bwa_chrMQ60PE.bam
samtools view -H AR3_C5_chrM_bwa_chrMQ60PE.bam > header.sam
samtools view  -@ 12 AR3_C5_chrM_bwa_chrMQ60PE.bam |awk '{if ($7 == "=") print $0}' | grep -v "MM:i:1" | cat header.sam - | samtools view -Sb -@ 12  - > AR3_C5_chrM_bwa_chrMQ60PE.unique.bam
samtools index AR3_C5_chrM_bwa_chrMQ60PE.unique.bam
samtools sort -@ 12 AR3_C5_chrM_bwa_chrMQ60PE.unique.bam -o AR3_C5_chrM_bwa_chrMQ60PE_filtered.bam
samtools index AR3_C5_chrM_bwa_chrMQ60PE_filtered.bam

# Step4: make re-align target regions
java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T RealignerTargetCreator  -nt 32 \
-I AR3_C4_chrM_bwa_chrMQ60PE_filtered.bam \
-o AR3_C4_chrM_bwa_chrMQ60PE_filtered.bam.target.intervals

nohup java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I AR3_C4_chrM_bwa_chrMQ60PE_filtered.bam \
-targetIntervals AR3_C4_chrM_bwa_chrMQ60PE_filtered.bam.target.intervals \
-o AR3_C4_chrM_bwa_chrMQ60PE_filtered.realign.bam &

samtools view -h AR3_C4_chrM_bwa_chrMQ60PE_filtered.realign.bam | awk '{if ($0 !~ /NM:i:[2-9]/ && $0 !~ /NM:i:[0-9][0-9]/) print $0}' | samtools view -b -o AR3_C4_chrM_bwa_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam -


java -Xmx8g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T RealignerTargetCreator  -nt 32 \
-I AR3_C5_chrM_bwa_chrMQ60PE_filtered.bam \
-o AR3_C5_chrM_bwa_chrMQ60PE_filtered.bam.target.intervals


nohup java -Xmx100g -jar /md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa  \
-T IndelRealigner -filterNoBases -maxReads 10000000  \
-I AR3_C5_chrM_bwa_chrMQ60PE_filtered.bam \
-targetIntervals AR3_C5_chrM_bwa_chrMQ60PE_filtered.bam.target.intervals \
-o AR3_C5_chrM_bwa_chrMQ60PE_filtered.realign.bam &

samtools view -h AR3_C5_chrM_bwa_chrMQ60PE_filtered.realign.bam | awk '{if ($0 !~ /NM:i:[2-9]/ && $0 !~ /NM:i:[0-9][0-9]/) print $0}' | samtools view -b -o AR3_C5_chrM_bwa_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam -
