# Step1: trim the reads form R1 and R2 

# trim fastq files 
## PBS configure
#PBS -N trim_fastq
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=6
#PBS -l mem=10G

cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/rawdata

cat AR3_C5.input|while read id;

do echo $id
arr=($id)
sample=${arr[0]}
input1=${arr[1]}
input2=${arr[2]}

java -jar /md01/nieyg/ori/biosoft/package/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
 -phred33 $input1 $input2 ../trim_data/$input1 ../trim_data/$sample_unpaired_R1.fq.gz ../trim_data/$input2 ../trim_data/$sample_unpaired_R3.fq.gz LEADING:35 TRAILING:35 CROP:100 MINLEN:35 AVGQUAL:30 -threads 6
done

# Step2: subseq the reads of R2 by header of R1 or R3

#!/bin/bash

# List of sample names to process
sample_names=("S6" "S7" "S8" "S9")

for sample_name in "${sample_names[@]}"; do

    # Decompress R1 and R2 FASTQ files
    gunzip -c "/data/R03/zhangwx/rawdata/mouse/Yage-mouse/AR3/AR3_C4/AR3_C4_${sample_name}_L001_R2_001.fastq.gz" > "AR3_C4_${sample_name}_L001_R2_001.fastq"
    gunzip -c "AR3_C4_${sample_name}_L001_R1_001.fastq.gz" > "AR3_C4_${sample_name}_L001_R1_001.fastq"

    # Use fastq_pair to pair the reads
    fastq_pair "AR3_C4_${sample_name}_L001_R1_001.fastq" "AR3_C4_${sample_name}_L001_R2_001.fastq"

    # Remove temporary files
    rm "AR3_C4_${sample_name}_L001_R2_001.fastq"
    rm "AR3_C4_${sample_name}_L001_R1_001.fastq"

    # Rename the paired file
    mv "AR3_C4_${sample_name}_L001_R2_001.fastq.paired.fq" "AR3_C4_${sample_name}_L001_R2_001.fastq"

    # Compress the paired file
    gzip "AR3_C4_${sample_name}_L001_R2_001.fastq"
done

sed -i 's/AR3_C4/AR3_C5/g' *.sh


# Step3: cellranger for trimed fastq 

#PBS -N masked_mouse
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC
cellranger-atac count --id=AR3_C4_add500G_100bp \
                        --reference=/md01/nieyg/ref/hard-mask/mm10_hard_masked \
                        --fastqs=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/trim_data/AR3_C4/ \
                        --sample=AR3_C4 \
                        --localcores=20 \
                        --localmem=64


#PBS -N masked_mouse
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC
cellranger-atac count --id=AR3_C5_add500G_100bp \
                        --reference=/md01/nieyg/ref/hard-mask/mm10_hard_masked \
                        --fastqs=/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/trim_data/AR3_C5/ \
                        --sample=AR3_C5 \
                        --localcores=20 \
                        --localmem=64

#Step4: mgatk for trimed fastq 

#PBS -N mgatk_AR3_C4
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/mgatk
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/outs/possorted_bam.bam \
-o  /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/mgatk \
-n AR3_C4_add500G_100bp -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 10 -bt CB -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/barcode.tsv

#PBS -N mgatk_AR3_C5
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/mgatk
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/outs/possorted_bam.bam \
-o  /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/mgatk \
-n AR3_C5_add500G_100bp -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 10 -bt CB -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk/barcode.tsv










