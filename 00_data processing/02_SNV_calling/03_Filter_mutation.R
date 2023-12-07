
# Step6: mpileup get file for snv calling

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
-f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa \
-q 10 -Q 10 AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.realign.bam > AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup

samtools mpileup  -l /md01/nieyg/ref/hard-mask/genome_modify/chrM.len \
-f /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa \
-q 10 -Q 10 AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.realign.bam > AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup

# make count table for each allele 
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup > AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup.count
/md01/nieyg/software/ATAC_mito_sc-master/src/pileup_inf_rj.pl AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup > AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup.count

# Step5: varscan call somatic mutation in chrM
varscan  pileup2snp AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup --min-coverage 3 --min-var-freq 0.0000001  --min-reads2 2 > AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.snv
varscan  pileup2snp AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup --min-coverage 3 --min-var-freq 0.0000001  --min-reads2 2 > AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.snv

varscan  mpileup2snp AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup --variants SNP --p-value 0.99 --output-vcf 1 --min-coverage 3 --min-var-freq 0.0001  --min-reads2 2 > AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.vcf
varscan  mpileup2snp AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup --variants SNP --p-value 0.99 --output-vcf 1 --min-coverage 3 --min-var-freq 0.0001  --min-reads2 2 > AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.vcf

# plot the mutation signature 
# Any neighbouring SNVs will be merged into DBS/MBS variants.
library(MutationalPatterns)
library(BSgenome)
head(available.genomes())
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
library(ref_genome, character.only = TRUE)
vcf_files <- list.files(path="/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/common_mutation/vcf/",pattern = ".vcf", full.names = TRUE)
sample_names <- c(  "AR3-C4", "AR3-C5")

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,group="circular", type = c("snv"))
time <- c( "AR3-C4", "AR3-C5")
#snv_grl <- get_mut_type(grl, type = "snv")
muts <- mutations_from_vcf(grl[[1]])
head(muts, 12)
library("gridExtra")
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)
pdf("./AR3_pseudobulk_96_profile.pdf",width=12,height=6)
plot_96_profile(mut_mat)
dev.off()

# plot the location signature 
# 1. coverage 
AR3_C4_data<- read.table("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mito_sort_out/AR3_C4_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup.count")
AR3_C5_data<- read.table("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mito_sort_out/AR3_C5_chrM_bwa_rmdup_chrMQ30PE_filtered.mpileup.count")
AR3_C4_data$sample<-"AR3_C4"
AR3_C5_data$sample<-"AR3_C5"

mdf<- rbind(AR3_C4_data,AR3_C5_data)
library(ggplot2)
pdf("./AR3_pseudobulk_location_coverage_profile.pdf",width=24,height=6)
ggplot(mdf, aes(x = V2, y = V4,color=sample)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "coverage") + 
  theme(legend.position = "none")

dev.off()

mdf$mutation<- (mdf$V9+mdf$V10)
pdf("./AR3_pseudobulk_location_alt_profile.pdf",width=24,height=6)
ggplot(mdf, aes(x = V2, y = V4,color=sample)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "coverage") + 
  theme(legend.position = "none")

ggplot(mdf, aes(x = V2, y = Freq,color=sample)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "# of alt") + 
  theme(legend.position = "none")


dev.off()


pdf("./AR3_pseudobulk_location_alt_profile16000.pdf",width=24,height=6)
ggplot(mdf, aes(x = V2, y = V4,color=sample)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "coverage") + xlim(16080,16150)+
  theme(legend.position = "none")

ggplot(mdf, aes(x = V2, y = Freq)) + 
  geom_bar(stat="identity") +  
  #pretty_plot(fontsize = 8)  + 
  #scale_color_fill(values = c("red", "blue")) +
  labs(x = "", y = "# of alt") + xlim(16080,16150)+
  theme(legend.position = "none")
dev.off()









