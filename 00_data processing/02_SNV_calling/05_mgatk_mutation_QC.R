# mgatk global 
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
combined<- readRDS("./03_all_celltype/03_recall_peak/AR3_integrated_all_celltype_annotated_recall_peak.rds")
# load mgatk output
AR3_C4_mito.data <- ReadMGATK(dir = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/02_Single_cell_SNV_calling/final")
AR3_C5_mito.data <- ReadMGATK(dir = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/02_Single_cell_SNV_calling/final")

# create an assay
AR3_C4_mito_counts<- AR3_C4_mito.data$counts
colnames(AR3_C4_mito_counts)<- paste("AR3_C4_",colnames(AR3_C4_mito_counts),sep="")
AR3_C5_mito_counts<- AR3_C5_mito.data$counts
colnames(AR3_C5_mito_counts)<- paste("AR3_C5_",colnames(AR3_C5_mito_counts),sep="")
all_counts<- cbind(AR3_C4_mito_counts,AR3_C5_mito_counts)

AR3_mito <- CreateAssayObject(counts = all_counts)
# Subset to cell present in the scATAC-seq assat
AR3_mito <- subset(AR3_mito, cells = colnames(combined))

# add assay and metadata to the seurat object
combined[["mito"]] <- AR3_mito

rownames(AR3_C4_mito.data$depth) <- paste("AR3_C4_",rownames(AR3_C4_mito.data$depth),sep="")
rownames(AR3_C5_mito.data$depth) <- paste("AR3_C5_",rownames(AR3_C5_mito.data$depth),sep="")
depth_all<- rbind(AR3_C4_mito.data$depth,AR3_C5_mito.data$depth)

combined <- AddMetaData(combined, metadata = depth_all, col.name = "mtDNA_depth")
Idents(combined)<- combined$detail_anno
pdf("./09_mgatk_for_clean_bam/all_cell_type_mtDNA_depth.pdf",width=16,height=8)
VlnPlot(combined, c("mtDNA_depth","nCount_mito","nFeature_mito"), pt.size = 0,ncol=1) + scale_y_log10()
dev.off()
saveRDS(combined,"./09_mgatk_for_clean_bam/all_cell_type_mgatk.rds")

# 1. mtDNA coverage among celltype 
combined<- readRDS("./09_mgatk_for_clean_bam/all_cell_type_mgatk.rds")
# load mgatk output
AR3_C4_cov.data <- read.table(gzfile("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/02_Single_cell_SNV_calling/final/AR3_C4_100bp_bwa.coverage.txt.gz"),sep=",")
AR3_C5_cov.data <- read.table(gzfile("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_add500G_100bp/02_Single_cell_SNV_calling/final/AR3_C5_100bp_bwa.coverage.txt.gz"),sep=",")

# add prefix in barcode 
AR3_C4_cov.data$V2<- paste("AR3_C4_",AR3_C4_cov.data$V2,sep="")
AR3_C5_cov.data$V2<- paste("AR3_C5_",AR3_C5_cov.data$V2,sep="")
cov.data<- rbind(AR3_C4_cov.data,AR3_C5_cov.data)
cov.data$celltype<- combined$Annotation[match(cov.data$V2,rownames(combined@meta.data))]
# Gental smooth of just 1 bp for plot aesthetics
library(roll)
smooth <- 5
cov.data$V3<- roll_mean(cov.data$V3, smooth)

mdf <- cov.data[,-2]
colnames(mdf)<- c("pos","coverage","celltype")
mdf[which(is.na(mdf$coverage)),]$coverage=0

df<- aggregate(mdf$coverage,by=list(mdf$pos,mdf$celltype),mean)
colnames(df)<- c("pos","celltype","coverage")

# Mean coverage X
aggregate(mdf$coverage,by=list(mdf$celltype),mean)
> aggregate(mdf$coverage,by=list(mdf$celltype),mean)
    Group.1         x
1        CM 149.47138
2        EC  76.12622
3        FB  80.12615
4       Epi  88.23092
5       SMC  54.85159
6  Pericyte  73.28716
7     MP_DC  37.91284
8         T  42.59819
9         B  54.14661
10    Glial  51.64424

# the coverage plot 
library(ComplexHeatmap)
library(GenomeInfoDb)
library(BuenColors) 
# Theme to remove any of the other plot riff raff
xxtheme <-   theme(
  axis.line = element_blank(),
  axis.ticks.y = element_blank(),   
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.margin = unit(c(0, 0, 0, 0), "cm"),
  plot.margin = unit(c(-0.35, -0.35, -0.35, -0.35), "cm")) 
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
# Visualize the rolled means
P1 <- ggplot(df, aes(x = pos, y = coverage, color = celltype)) + 
  geom_line() +  expand_limits(y = c(-5, 4)) +
  pretty_plot(fontsize = 8)  + scale_color_manual(values = myUmapcolors[1:11]) +
  coord_polar(direction = 1) + labs(x = "", y = "log2 Coverage") + scale_y_log10() +
  theme(legend.position = "none") 
cowplot::ggsave2(P1, file = "./09_mgatk_for_clean_bam/rollMean_coverage.pdf", width = 4, height = 4)

# IdentifyVariants
# filter cells based on mitochondrial depth
crc <- subset(combined, mtDNA_depth >= 10)

refAllele<- read.table("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_add500G_100bp/02_Single_cell_SNV_calling/final/chrM_refAllele.txt")
colnames(refAllele)<- c("pos","ref")
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = refAllele)

pdf("./09_mgatk_for_clean_bam/VariantPlot_global.pdf")
VariantPlot(variants = variable.sites)
dev.off()

# all celltype mutation 
# Establish a filtered data frame of variants based on this processing

# plot the frequency and sequening depth distribution:
pdf("./09_mgatk_for_clean_bam/QC_variable_sites_density.pdf",width=8,height=4)
plot(density(variable.sites$mean_coverage))
plot(density(variable.sites$mean*100),xlim=c(0,0.1))
plot(density(variable.sites$vmr),xlim=c(0,0.15))
plot(density(variable.sites$position),xlim=c(0,16299))
plot(density(variable.sites$n_cells_conf_detected),xlim=c(0,50))
plot(density(na.omit(variable.sites$strand_correlation)))
dev.off()

cutoff_n_cells_conf_detected<- as.matrix(summary(variable.sites$n_cells_conf_detected))[2]
cutoff_mean_coverage<- as.matrix(summary(variable.sites$mean_coverage))[2]
if(cutoff_mean_coverage>40){cutoff_mean_coverage=40}else{cutoff_mean_coverage=cutoff_mean_coverage}
if(cutoff_n_cells_conf_detected<5){cutoff_n_cells_conf_detected=5}else{cutoff_n_cells_conf_detected=cutoff_n_cells_conf_detected}

# step1: Consider strand balance(DNA damage casued C-->A and G-->T error) and the VMR(rm the germline mutation) 
# Step2: Restriction of minimum frequency and sequening depth:
high.conf <- subset(
  variable.sites, 
  subset = n_cells_conf_detected >= cutoff_n_cells_conf_detected &
    strand_correlation >= 0.65 &
    mean_coverage >= cutoff_mean_coverage &
    vmr > 0.01
)

high.conf<- high.conf[order(high.conf$mean,decreasing=T),]
high.conf[,c(1,2,5)]
crc <- AlleleFreq(
    object = crc,
    variants = high.conf$variant,
    assay = "mito"
  )

saveRDS(crc,"./09_mgatk_for_clean_bam/all_cell_type_mgatk_alleles.rds")

write.csv(high.conf,"./09_mgatk_for_clean_bam/high_conf_alleles_info.csv")
# After check the alleles 
combined<- readRDS("./09_mgatk_for_clean_bam/all_cell_type_mgatk_alleles.rds")

# Step1: QC metric for mutation 

Idents(combined)<- combined$detail_anno
combined$mito_reads_rate<- (combined$mitochondrial/combined$total)*100
library(scCustomize)
pdf("./09_mgatk_for_clean_bam/all_cell_type_mtDNA_depth.pdf",width=20,height=10)
VlnPlot(combined, c("mito_reads_rate","mtDNA_depth","nCount_alleles","nFeature_alleles"),col=myUmapcolors, pt.size = 0,ncol=1) + scale_y_log10()
VlnPlot(combined,log = TRUE, c("mito_reads_rate","mtDNA_depth","nCount_alleles","nFeature_alleles"), pt.size = 0,ncol=1) + scale_y_log10()
Stacked_VlnPlot(seurat_object = combined, features = c("mito_reads_rate","mtDNA_depth","nCount_mito","nFeature_mito","nCount_alleles","nFeature_alleles"), x_lab_rotate = TRUE,colors_use = myUmapcolors)
dev.off()

# Step2: VAF distribution
# global VAF distribution of mutation 
# the mutation freq barplot 
DefaultAssay(combined)<- "alleles"
data<- as.numeric(GetAssayData(combined))
# all mutation freq in single cell not all cells(remove the 0)
data<- data[data!=0]
data<- data.frame(variance=data)
data$variance<- data$variance;
pdf("./09_mgatk_for_clean_bam/All_celltype_mutation_freq.pdf",width=5,height=3)
  p=ggplot(data,aes(x=variance))+
  geom_histogram(
                 binwidth = 0.1,
                 fill="#69b3a2",##69b3a2
                 color="#e9ecef",##e9ecef
                 alpha=0.9,
                 breaks=seq(0,1,0.1))+ 
  theme_bw()+
  labs(x="Frequence",y="Count",title="All_celltype_mutation_freq")+
  scale_x_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.2))
  p
dev.off()

# VAF distribution of mutation in difference celltype 
# the mutation freq barplot 
Idents(combined)<- combined$Annotation
pdf("./09_mgatk_for_clean_bam/celltype_mutation_freq.pdf",width=10,height=6)
for (j in levels(combined)){
	print(j);
	obj<- subset(combined,idents=j)
	variable.sites <- IdentifyVariants(obj, assay = "mito", refallele = refAllele)
    # select the cutoff 
  cutoff_n_cells_conf_detected<- as.matrix(summary(variable.sites$n_cells_conf_detected))[2]
  cutoff_mean_coverage<- as.matrix(summary(variable.sites$mean_coverage))[2]
  if(cutoff_mean_coverage>40){cutoff_mean_coverage=40}else{cutoff_mean_coverage=cutoff_mean_coverage}
  if(cutoff_n_cells_conf_detected<5){cutoff_n_cells_conf_detected=5}else{cutoff_n_cells_conf_detected=cutoff_n_cells_conf_detected}
  
  # step1: Consider strand balance(DNA damage casued C-->A and G-->T error) and the VMR(rm the germline mutation) 
  # Step2: Restriction of minimum frequency and sequening depth:
  high.conf <- subset(
    variable.sites, 
    subset = n_cells_conf_detected >= cutoff_n_cells_conf_detected &
      strand_correlation >= 0.65 &
      mean_coverage >= cutoff_mean_coverage &
      vmr > 0.01
  )
  p1 <- VariantPlot(variants = variable.sites,concordance.threshold=0.65)
  print(p1)
  if(nrow(high.conf)>1){
  crc <- AlleleFreq(
    object = obj,
    variants = high.conf$variant,
    assay = "mito"
  )
  DefaultAssay(crc) <- "alleles"
  data<- GetAssayData(crc)
  d_varition1=as.vector(data)
  d_varition1=d_varition1[which(d_varition1!=0)]
  data<- data.frame(variance=d_varition1)
  data$variance<- data$variance;
  p=ggplot(data,aes(x=variance))+
  geom_histogram(
                 binwidth = 0.1,
                 fill="#69b3a2",##69b3a2
                 color="#e9ecef",##e9ecef
                 alpha=0.9,
                 breaks=seq(0,1,0.1))+ 
  theme_bw()+
  labs(x="Frequence",y="Count",title=j)+
  theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size=18,colour="black", face = "bold"),          
        axis.text.y = element_text(size=18,colour="black", face = "bold"),
        axis.title.x = element_text(size=14, face = "bold"), 
        axis.title.y = element_text(size=14, face = "bold"),
        panel.background = element_blank(),
        line = element_line(size=1),
        axis.line = element_line(size =1.0,colour = "black"),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank())+scale_x_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.2))
    print(p)
}}
dev.off()

# Step3: mutation lacation 

# plot the location signature 
library(ggplot2)

high.conf<- read.csv("./09_mgatk_for_clean_bam/high_conf_alleles_info.csv")
head(high.conf)
tmp_data<- as.data.frame(table(high.conf$position))

pdf("./09_mgatk_for_clean_bam/AR3_pseudobulk_location_alt_profile.pdf",width=24,height=6)

ggplot(tmp_data, aes(x = Var1, y = Freq)) + 
  geom_line() +  
  #pretty_plot(fontsize = 8)  + 
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "", y = "# of alt") + 
  theme(legend.position = "none")
dev.off()

# Step4: mutation signature 
# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}
library(data.table)
# Process 3 digit signature based on letters
ref_all <- fread("/data/R02/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/final/chrM_refAllele.txt")
colnames(ref_all) <- c("pos", "ref")
ref_all$ref <- toupper(ref_all$ref)
l <- as.character(ref_all$ref)

# Gs happen to be at the first and last position
ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

# Remove Ns
ref_all <- ref_all[!grepl("N", ref_all$three),]

# Make every possible mutation
ref_all_long <- rbind(ref_all,ref_all, ref_all,ref_all)
ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]

# add some meta data
ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

# A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
table(ref_all$ref) # so the reference strand is light (more C/T)
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")

# Change to C/T as ref allele
ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)

# Annotate with called variants
called_variants <- high.conf$variant
ref_all_long$called <- ref_all_long$variant %in% called_variants

# Compute changes in expected/observed
total <- dim(ref_all_long)[1]
total_called <- sum(ref_all_long$called)
library(dplyr)
library(tidyr)
prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
  summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
  mutate(fc_called = observed_prop_called/expected_prop)

prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)

# Visualize
library(tidyverse)
library(prettyGraphs)
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom",
    axis.text.x =element_text(angle=90,hjust=1,size=2)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "./09_mgatk_for_clean_bam/all_mito_signature.pdf", width = 4, height = 2.4)



# some high VAF mutation distribution in UMAP 
high.conf<- high.conf[order(high.conf$mean,decreasing=TRUE),]
mean_VAF_high_mutation<- high.conf$variant[1:9]
high.conf<- high.conf[order(high.conf$n_cells_conf_detected,decreasing=TRUE),]
cell_high_mutation<- high.conf$variant[1:9]


DefaultAssay(combined) <- "alleles"
pdf("./09_mgatk_for_clean_bam/all_cell_type_mutation_featureplot.pdf",width=9,height=9)
FeaturePlot(
  object = combined,
  features = mean_VAF_high_mutation,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 3
) & NoLegend()
FeaturePlot(
  object = combined,
  features = cell_high_mutation,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 3
) & NoLegend()
dev.off()

Idents(combined)<- combined$Annotation

pdf("./09_mgatk_for_clean_bam/all_cell_type_top100mutation_heatmap.pdf",width=9,height=12)
DoHeatmap(combined, features =  rownames(combined), slot = "data", disp.max = 1)
dev.off()






















library(SummarizedExperiment)
library(Matrix)

# Function that quickly computes the allele frequency matrix from a summarized experiment mgatk object
computeAFMutMatrix <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
}


library(data.table)
library(dplyr)
# Pretty simple plot... plotting distribution of heteroplasmy for top variants
pt1_vars <- readRDS("../output/PT1_specificVariants_forSupplement.rds")
pt2_vars <- readRDS("../output/PT2_specificVariants_forSupplement.rds")
plot_df <- data.frame(
  het = c(pt1_vars$X5140G.A, pt1_vars$X14858G.A,
          pt1_vars$X1872T.C, pt1_vars$X1260A.G,
          pt2_vars$X12980G.A, pt2_vars$X4853G.A) * 100,
  mut = c(rep("5140G>A", dim(pt1_vars)[1]), rep("14858G>A", dim(pt1_vars)[1]), 
          rep("1872T>C", dim(pt1_vars)[1]), rep("1260A>G", dim(pt1_vars)[1]), 
          rep("12980G>A", dim(pt2_vars)[1]), rep("4853G>A", dim(pt2_vars)[1]))
)
plot_df$mut <- factor(as.character(plot_df$mut), levels = c("5140G>A", "14858G>A", "4853G>A", "1872T>C", "1260A>G", "12980G>A"))
p1 <- ggplot(plot_df, aes(x = het)) +
  geom_histogram(binwidth = 10, fill = "black") + facet_wrap(~mut ) +
  pretty_plot(fontsize = 7) 
cowplot::ggsave2(p1, file = "../plots/mut_histograms.pdf",width = 3.5, height = 1.8)

# For red number in the supplement
plot_df %>% group_by(mut) %>% summarize(count = sum(het >= 90))






