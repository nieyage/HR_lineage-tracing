# Step1: Extract chrM bam file 

#PBS -N scATAC
#PBS -j oe
#PBS -q batch
#PBS -o example.stdout
#PBS -e example2.stdout
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=2
#PBS -l mem=10000m

samtools view -h /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/mito_sort_out/AR3_C4_chrM_bwa.bam > 1.sam
sed -i "s/chrM/chrMT/g"  1.sam
sed -i "27,87d"  1.sam
samtools view -bS 1.sam > possorted_bam_chrM_new.bam
samtools index possorted_bam_chrM_new.bam
rm 1.sam

samtools index AR3_C4_chrM_bwa.bam

nohup java -Xmx8g -jar /public/home/chenbzh5/Tools/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa \
-T RealignerTargetCreator \
-nt 32 -I AR3_C4_chrM_bwa.bam -o AR3_C4_chrM_bwa.target.intervals &

java -Xmx10g -jar /public/home/chenbzh5/Tools/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa \
-T IndelRealigner -maxReads 10000000 \
-I AR3_C4_chrM_bwa.bam -targetIntervals AR3_C4_chrM_bwa.target.intervals \
-o AR3_C4_chrM_bwa_realign.bam
samtools sort -t CR AR3_C4_chrM_bwa_realign.bam -o AR3_C4_chrM_bwa_realign_sorted.bam

# For AR3_C5

samtools view -h /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/mito_sort_out/AR3_C5_chrM_bwa.bam > 2.sam
sed -i "s/chrM/chrMT/g"  2.sam
sed -i "27,87d"  2.sam
samtools view -bS 2.sam > AR3_C5_possorted_bam_chrM_new.bam
samtools index AR3_C5_possorted_bam_chrM_new.bam
rm 2.sam

samtools index AR3_C5_chrM_bwa.bam

nohup java -Xmx8g -jar /public/home/chenbzh5/Tools/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa \
-T RealignerTargetCreator \
-nt 32 -I AR3_C5_chrM_bwa.bam -o AR3_C5_chrM_bwa.target.intervals &

java -Xmx10g -jar /public/home/chenbzh5/Tools/GenomeAnalysisTK_3.5-0.jar \
-R /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa \
-T IndelRealigner -maxReads 10000000 \
-I AR3_C5_chrM_bwa.bam -targetIntervals AR3_C5_chrM_bwa.target.intervals \
-o AR3_C5_chrM_bwa_realign.bam

samtools sort -t CR AR3_C5_chrM_bwa_realign.bam -o AR3_C5_chrM_bwa_realign_sorted.bam

# Step2: 04_Split_bam.py; Split bam file by cell barcord(tag'CR')

#!/bin/python
import pysam
import os
import subprocess
# file to split on
unsplit_file = "AR3_C5_chrM_bwa_realign_sorted.bam"
# where to place output files
out_dir = "./AR3_C5_splitbam/"

# variable to hold barcode index
CB_hold = 'unset'
itr = 0
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile( unsplit_file, "rb")
lit = []
for read in samfile.fetch( until_eof=True):
    # barcode itr for current read
    try:
        CB_itr = read.get_tag('CB')
    except:
        continue
    CB_itr = read.get_tag('CB')
    # if change in barcode or first line; open new file  
    if( CB_itr!=CB_hold or itr==0):
        # close previous split file, only if not first read in file
        if( itr!=0):
            if len(lit)>1:
                split_file = pysam.AlignmentFile( out_dir + "{}.bam".format(CB_hold), "wb", template=samfile)
                for kk in lit:
                    split_file.write(kk)
                split_file.close()
                lit = []
            lit = []
        CB_hold = CB_itr
        itr+=1
       # print(CB_hold)

    lit.append(read)
   # print(len(lit))
        #split_file = pysam.AlignmentFile( out_dir + "CB_{}.bam".format(CB_hold), "wb", template=samfile)
    
split_file.close()

samfile.close()


# Step3: Single cell SNV 
#coding=utf-8
import os
import time

# Tools to identify mpileUP's reference genome (Human_G1K_V37), target region (chrM), and SNV Calling
reference_chrM = "/public/home/chenbzh5/DB/hg19/human_g1k_v37_chrM/chrM.len"
reference_gemone = "/public/home/tangzj/PB_ATAC/GSE142745/SRR10804559/refdata-cellranger-atac-hg19-1.2.0/fasta/genome.fa"
picard = "/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
GATK = "/public/home/chenbzh5/Tools/GenomeAnalysisTK_3.5-0.jar"
varscan = "/public/home/chenbzh5/Tools/VarScan.v2.3.7.jar"

# Set input and outout path
input_dir = "/md01/tangzj/project/scRNAseq/GSE159929/bladder/bladder/splitbam"
output_dir = "/public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/rna/bladder/allcells"
pool_output_dir = ""
print("mkdir " + output_dir+"\n")
print("mkdir " + pool_output_dir)

# Read bam files named by cellbarcords.
dirs = os.listdir(input_dir)
cellname = []
for j in dirs:
	cellname.append(j.split(".")[0].strip())       

# GATK was used to remove PCR duplicates and re-comparisons

#Pileup_inf_rj.pl generate count file.
#Varscan generate SNV file.
# for i in range(0, int(len(cellname)/40) + 1):
for number in range(0,len(cellname)):
    #if  os.path.exists("/public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/TZJ/cancer_cellrange_output/CCL2+_1cellranger_output/SRR10804671/outs/splitbam/"+cellname[number]+"/"+cellname[number]+".snv"):
       # print(str(cellname[number]))
     #   continue
    pbs ='''#PBS -N scATAC\n#PBS -j oe\n#PBS -q batch\n#PBS -o example.stdout\n#PBS -e example2.stdout\n#PBS -S /bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=6000m\nsource /public/home/jinxu/software/cellranger-atac-1.1.0/sourceme.bash\nmodule load samtools
'''
    f = open("./temp3.sh","w")
    f.write(pbs)
    f.write("mkdir " + os.path.join(output_dir, cellname[number])+"\n")
    f.write("cd "+output_dir+"/"+cellname[number]+"\n")
    f.write("cp "+input_dir+"/"+cellname[number]+".bam ./\n")

    input_bam = os.path.join(input_dir,cellname[number] + ".bam")
    # picard remove duplicates
    rmdup_bam = os.path.join(output_dir, cellname[number], cellname[number] + ".rmdup.bam")
    metric = os.path.join(output_dir, cellname[number], cellname[number] + ".metrics")
    f.write("java -Xmx4g -jar " + picard + \
                      " MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true INPUT=" + os.path.join(output_dir, cellname[number], cellname[number] + ".bam") + \
                      " OUTPUT=" + rmdup_bam + \
                      " METRICS_FILE=" + metric+"\n")
    f.write("samtools mpileup -l "+reference_chrM+" -q 30 -Q 30 -f "+reference_gemone+"  -x "+cellname[number]+".rmdup.bam > "+cellname[number]+".rmdup.mpileup\n")
    f.write("java -Xmx4g  -jar /public/home/chenbzh5/Tools/VarScan.v2.3.7.jar  pileup2snp "+cellname[number]+".rmdup.mpileup --min-var-freq 0.01  --min-reads2 2 >  "+cellname[number]+".snv\n")
    mpileup = os.path.join(output_dir, cellname[number], cellname[number] + ".rmdup.mpileup")
    count = os.path.join(output_dir, cellname[number], cellname[number] + ".counts")
    f.write("/md01/jinxu/bin/pileup_inf_rj.pl " + mpileup + " > " + count)
#if number%20==0 and number!=0:
    #	time.sleep(1500)
    f.close()
    os.system("qsub temp3.sh")








