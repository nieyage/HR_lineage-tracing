# Method1: mgatk 
#PBS -N mgatk_AR3_C4
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/02_Single_cell_SNV_calling
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/01_clean_bam/not_rmdup/AR3_C4_chrM_bwa_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam \
-o  /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/02_Single_cell_SNV_calling \
-n AR3_C4_100bp_bwa -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 10 -bt CB -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/barcode.tsv

#PBS -N mgatk_AR3_C5
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/02_Single_cell_SNV_calling
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/01_clean_bam/not_rmdup/AR3_C5_chrM_bwa_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam \
-o  /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/02_Single_cell_SNV_calling \
-n AR3_C5_100bp_bwa -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 10 -bt CB -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk/barcode.tsv

#PBS -N mgatk_AR3_C4_not_rmdup
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/02_Single_cell_SNV_calling
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/01_clean_bam/not_rmdup/AR3_C4_chrM_bwa_chrMQ60PE_filtered.realign_filtered_unique_mismatches.bam \
-o  /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/02_Single_cell_SNV_calling \
-n AR3_C4_100bp_bwa -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 10 -bt CB -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/barcode.tsv


# Method2: Single cell ATAC-seq were processed similarly to the bulk ATAC-seq, 
#          taking each individual cell as one sample.

# ref: https://github.com/ChangLab/ATAC_mito_sc

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
combined<- readRDS("./03_all_celltype/03_recall_peak/AR3_integrated_all_celltype_annotated_recall_peak.rds")
Idents(combined)<- combined$orig.ident
AR3_C4<- subset(combined,idents="AR3_C4")
AR3_C5<- subset(combined,idents="AR3_C5")

# for CM
AR3_C4_barcode<- gsub("AR3_C4_","",rownames(AR3_C4@meta.data[which(AR3_C4$Annotation=="CM"),]))
AR3_C5_barcode<- gsub("AR3_C5_","",rownames(AR3_C5@meta.data[which(AR3_C5$Annotation=="CM"),]))
dir.create(paste("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/","CM",sep=""))
write.table(AR3_C4_barcode,paste0("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/","CM","/AR4_barcode.txt",sep=""),row.name=F,col.names=F)
write.table(AR3_C5_barcode,paste0("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/","CM","/AR5_barcode.txt",sep=""),row.name=F,col.names=F)

# for every cell type 
Idents(combined)<- combined$Annotation
celltype<- levels(combined)
for (i in celltype){
	print(i)
	AR3_C4_barcode<- gsub("AR3_C4_","",rownames(AR3_C4@meta.data[which(AR3_C4$Annotation==i),]))
    AR3_C5_barcode<- gsub("AR3_C5_","",rownames(AR3_C5@meta.data[which(AR3_C5$Annotation==i),]))
    dir.create(paste("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/",i,sep=""))
    write.table(AR3_C4_barcode,paste0("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/",i,"/AR4_barcode.txt",sep=""),row.name=F,col.names=F)
    write.table(AR3_C5_barcode,paste0("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/",i,"/AR5_barcode.txt",sep=""),row.name=F,col.names=F)
}


# in shell
cd CM
sed -i 's/"//g' *.txt
nohup /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/bottleneck_ppl/AR3_C4_100bp_possorted_chrM_realign.bam --cell-barcodes AR4_barcode.txt --cores 20 --out-bam AR3_C4_barcode.bam &
nohup /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/bottleneck_ppl/AR3_C5_100bp_possorted_chrM_realign.bam --cell-barcodes AR5_barcode.txt --cores 20 --out-bam AR3_C5_barcode.bam &

mkdir splitbam
cp ../AR3_C4_split_bam.py .
cp ../AR3_C5_split_bam.py .
python AR3_C4_split_bam.py
python AR3_C5_split_bam.py


#!/bin/python
import pysam
import os
import subprocess
# file to split on
unsplit_file = "AR3_C4_barcode.bam"
# where to place output files
out_dir = "./splitbam/"

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

#PBS -N split_bam 
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam

celltypes=("CM"       "EC"       "FB"       "Epi"      "SMC"      "Pericyte" "MP_DC"    "T"        "B"        "Glial")
for celltype in "${celltypes[@]}"
do
    cd "$celltype" || exit 1
    sed -i 's/"//g' *.txt
    /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/bottleneck_ppl/AR3_C4_100bp_possorted_chrM_realign.bam --cell-barcodes AR4_barcode.txt --out-bam AR3_C4_barcode.bam
    /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/bottleneck_ppl/AR3_C5_100bp_possorted_chrM_realign.bam --cell-barcodes AR5_barcode.txt --out-bam AR3_C5_barcode.bam
    mkdir splitbam
    cp "../AR3_C4_split_bam.py" .
    cp "../AR3_C5_split_bam.py" .
    python AR3_C4_split_bam.py
    python AR3_C5_split_bam.py
    cd ..
done


# for every cell type 
Idents(combined)<- combined$detail_anno
celltype<- levels(combined)
for (i in celltype){
	print(i)
	AR3_C4_barcode<- gsub("AR3_C4_","",rownames(AR3_C4@meta.data[which(AR3_C4$detail_anno==i),]))
    AR3_C5_barcode<- gsub("AR3_C5_","",rownames(AR3_C5@meta.data[which(AR3_C5$detail_anno==i),]))
    dir.create(paste("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_cluster_bam/",i,sep=""))
    write.table(AR3_C4_barcode,paste0("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_cluster_bam/",i,"/AR4_barcode.txt",sep=""),row.name=F,col.names=F)
    write.table(AR3_C5_barcode,paste0("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_cluster_bam/",i,"/AR5_barcode.txt",sep=""),row.name=F,col.names=F)
}


# in shell

#PBS -N split_bam_cluster
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_cluster_bam

celltypes=( "27:EC"       "19:FB"       "4:EC"        "26:SMC"      "11:FB" "17:EC"       "3:MP_DC"     "29:Epi"      "0:MP_DC"     "7:FB" "2:FB"        "14:Pericyte" "5:CM"        "1:MP_DC"     "6:CM" "18:NA"       "10:CM"       "15:FB"       "8:EC"        "32:FB" "9:CM"        "12:FB"       "16:FB"       "30:CM"       "20:MP_DC" "33:EC"       "31:Epi"      "25:MP_DC"    "37:MP_DC"    "35:B" "21:Epi"      "13:FB"       "24:FB"       "34:SMC"      "23:T" "39:EC"       "22:EC"       "28:MP_DC"    "41:MP_DC"    "36:T" "38:Glial"    "40:MP_DC")
for celltype in "${celltypes[@]}"
do
    cd "$celltype" || exit 1
    sed -i 's/"//g' *.txt
    # 执行 subset-bam_linux 命令
    /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/bottleneck_ppl/AR3_C4_100bp_possorted_chrM_realign.bam --cell-barcodes AR4_barcode.txt --cores 20 --out-bam AR3_C4_barcode.bam
    /md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/bottleneck_ppl/AR3_C5_100bp_possorted_chrM_realign.bam --cell-barcodes AR5_barcode.txt --cores 20 --out-bam AR3_C5_barcode.bam
    # 创建目录并复制相关文件
    mkdir splitbam
    cp "../AR3_C4_split_bam.py" .
    cp "../AR3_C5_split_bam.py" .
    python AR3_C4_split_bam.py
    python AR3_C5_split_bam.py
    cd ..
done


# single cell SNV calling by bottleneck ppl 

#coding=utf-8
import os
import time

# Tools to identify mpileUP's reference genome (genome_modify mm10), target region (chrM), and SNV Calling
reference_chrM = "/md01/nieyg/ref/hard-mask/genome_modify/chrM.len"
reference_gemone = "/md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa"
picard = "/public/home/chenbzh5/Tools/picard-tools-2.4.1/picard.jar"
GATK = "/md01/nieyg/software/GenomeAnalysisTK_3.5-0.jar"

# Set input and outout path
input_dir = "/md01/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/CM/splitbam"
output_dir = "/md01/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/CM/allcells"
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
    #f.write(pbs)
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
    f.write("varscan  pileup2snp "+cellname[number]+".rmdup.mpileup --min-var-freq 0.0005  --min-reads2 2 >  "+cellname[number]+".snv\n")
    mpileup = os.path.join(output_dir, cellname[number], cellname[number] + ".rmdup.mpileup")
    count = os.path.join(output_dir, cellname[number], cellname[number] + ".counts")
    f.write("/md01/jinxu/bin/pileup_inf_rj.pl " + mpileup + " > " + count)
#if number%20==0 and number!=0:
    #	time.sleep(1500)
    f.close()
    os.system("bash temp3.sh")


#PBS -N snvcalling
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=64G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/02_AR3_add500G/09_celltype_bam/CM
python SNVcalling.py







 cp ../EC/step2_SNVcalling.pbs .
 cp ../EC/SNVcalling.py .
 mkdir allcells

 qsub step2_SNVcalling.pbs










