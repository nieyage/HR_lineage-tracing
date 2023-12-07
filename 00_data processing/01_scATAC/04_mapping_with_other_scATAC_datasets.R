# transfer the previous dataset to signac object 
# run the cellranger for rawdatasets
#PBS -N AR7
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20
#PBS -l mem=30G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data
cellranger-atac count --id=AR7 \
                        --reference=/md01/nieyg/ref/hard-mask/mm10_hard_masked \
                        --fastqs=/md01/nieyg/scATAC-ArchR/uploaddata/Rawdata/AR7 \
                        --sample=AR7-SI-NA-E1 \
                        --localcores=20 \
                        --localmem=64




/md01/nieyg/scATAC-ArchR/uploaddata/ProcessedData

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
setwd("/md01/nieyg/project/lineage_tracing/heart_regeneration/03_HR_scATAC")
# combine peaksets 
# read in peak sets
peaks.AR3 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR3/outs/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.AR7 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR7/outs/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.AR14 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR14/outs/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.P4 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P4/outs/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.P8 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P8/outs/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.P15 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P15/outs/peaks.bed",
  col.names = c("chr", "start", "end"))
# convert to genomic ranges
gr.AR3 <- makeGRangesFromDataFrame(peaks.AR3)
gr.AR7 <- makeGRangesFromDataFrame(peaks.AR7)
gr.AR14 <- makeGRangesFromDataFrame(peaks.AR14)
gr.P4 <- makeGRangesFromDataFrame(peaks.P4)
gr.P8 <- makeGRangesFromDataFrame(peaks.P8)
gr.P15 <- makeGRangesFromDataFrame(peaks.P15)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.AR3, gr.AR7,gr.AR14,gr.P4, gr.P8,gr.P15))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
md.AR3 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR3/outs/singlecell.csv",
  stringsAsFactors = FALSE,sep = ",",header = TRUE, row.names = 1)[-1, ]
md.AR7 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR7/outs/singlecell.csv",
  stringsAsFactors = FALSE,sep = ",",header = TRUE, row.names = 1)[-1, ]
md.AR14 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR14/outs/singlecell.csv",
  stringsAsFactors = FALSE,sep = ",",header = TRUE, row.names = 1)[-1, ]
md.P4 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P4/outs/singlecell.csv",
  stringsAsFactors = FALSE,sep = ",",header = TRUE, row.names = 1)[-1, ]
md.P8 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P8/outs/singlecell.csv",
  stringsAsFactors = FALSE,sep = ",",header = TRUE, row.names = 1)[-1, ]
md.P15 <- read.table(
  file = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P15/outs/singlecell.csv",
  stringsAsFactors = FALSE,sep = ",",header = TRUE, row.names = 1)[-1, ]

# perform an initial filtering of low count cells
md.AR3 <- md.AR3[md.AR3$passed_filters > 500, ]
md.AR7 <- md.AR3[md.AR7$passed_filters > 500, ]
md.AR14 <- md.AR3[md.AR14$passed_filters > 500, ]
md.P4 <- md.P4[md.P4$passed_filters > 500, ]
md.P8 <- md.P4[md.P8$passed_filters > 500, ]
md.P15 <- md.P4[md.P15$passed_filters > 500, ]

# create fragment objects
frags.AR3 <- CreateFragmentObject(
  path = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR3/outs/fragments.tsv.gz",
  cells = rownames(md.AR3))
frags.AR7 <- CreateFragmentObject(
  path = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR7/outs/fragments.tsv.gz",
  cells = rownames(md.AR7))
frags.AR14 <- CreateFragmentObject(
  path = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR14/outs/fragments.tsv.gz",
  cells = rownames(md.AR14))
frags.P4 <- CreateFragmentObject(
  path = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P4/outs/fragments.tsv.gz",
  cells = rownames(md.P4))
frags.P8 <- CreateFragmentObject(
  path = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P8/outs/fragments.tsv.gz",
  cells = rownames(md.P8))
frags.P15 <- CreateFragmentObject(
  path = "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P15/outs/fragments.tsv.gz",
  cells = rownames(md.P15))

AR3.counts <- FeatureMatrix(fragments = frags.AR3,features = combined.peaks,cells = rownames(md.AR3))
AR7.counts <- FeatureMatrix(fragments = frags.AR7,features = combined.peaks,cells = rownames(md.AR7))
AR14.counts <- FeatureMatrix(fragments = frags.AR14,features = combined.peaks,cells = rownames(md.AR14))
AR3_assay <- CreateChromatinAssay(AR3.counts, fragments = frags.AR3)
AR7_assay <- CreateChromatinAssay(AR7.counts, fragments = frags.AR7)
AR14_assay <- CreateChromatinAssay(AR14.counts, fragments = frags.AR14)
P4.counts <- FeatureMatrix(fragments = frags.P4,features = combined.peaks,cells = rownames(md.P4))
P8.counts <- FeatureMatrix(fragments = frags.P8,features = combined.peaks,cells = rownames(md.P8))
P15.counts <- FeatureMatrix(fragments = frags.P15,features = combined.peaks,cells = rownames(md.P15))
P4_assay <- CreateChromatinAssay(P4.counts, fragments = frags.P4)
P8_assay <- CreateChromatinAssay(P8.counts, fragments = frags.P8)
P15_assay <- CreateChromatinAssay(P15.counts, fragments = frags.P15)

AR3 <- CreateSeuratObject(AR3_assay, assay = "ATAC", meta.data=md.AR3)
AR7 <- CreateSeuratObject(AR7_assay, assay = "ATAC", meta.data=md.AR7)
AR14 <- CreateSeuratObject(AR14_assay, assay = "ATAC", meta.data=md.AR14)
P4 <- CreateSeuratObject(P4_assay, assay = "ATAC", meta.data=md.P4)
P8 <- CreateSeuratObject(P8_assay, assay = "ATAC", meta.data=md.P8)
P15 <- CreateSeuratObject(P15_assay, assay = "ATAC", meta.data=md.P15)

saveRDS(AR3,"./AR3_raw.rds")
saveRDS(AR7,"./AR7_raw.rds")
saveRDS(AR14,"./AR14_raw.rds")
saveRDS(P4,"./P4_raw.rds")
saveRDS(P8,"./P8_raw.rds")
saveRDS(P15,"./P15_raw.rds")


# extract gene annotations from EnsDb

library(GenomicFeatures)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(Signac)
AR3<- readRDS("./AR3_raw.rds")
AR7<- readRDS("./AR7_raw.rds")
AR14<- readRDS("./AR14_raw.rds")
P4<- readRDS("./P4_raw.rds")
P8<- readRDS("./P8_raw.rds")
P15<- readRDS("./P15_raw.rds")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
objList <- c(AR3,AR7,AR14,P4,P8,P15)

# calculate the score of NS and TSS
for (i in seq_len(length(objList))) {
    Annotation(objList[[i]]) <- annotations # add the gene information to the object
    objList[[i]] <- NucleosomeSignal(objList[[i]])
    objList[[i]] <- TSSEnrichment(objList[[i]],fast=FALSE)
    }
sample<- c("AR3","AR7","AR14","P4","P8","P15")
library(ggplot2)
# plot TSS and fragrament distribution plot 
for (i in seq_len(length(objList))) {
    dataset.name<-sample[i]
    print(dataset.name)
    filename <- paste0("./01_qc/TSS_distribution_",dataset.name,".pdf")
    pdf(filename)
    objList[[i]]$high.tss<-ifelse(objList[[i]]$TSS.enrichment > 2, 'High', 'Low')
    TSS<-TSSPlot(objList[[i]], group.by = 'high.tss') + labs(title = "NE")
    objList[[i]]$NS <- ifelse(objList[[i]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
    Frag<-FragmentHistogram(object = objList[[i]], group.by = 'NS',region = "chr1-1-195471971")+ labs(title = "NE")
    print(TSS);
    print(Frag);
    dev.off()
    objList[[i]]$pct_reads_in_peaks <- objList[[i]]$peak_region_fragments / objList[[i]]$passed_filters * 100
    objList[[i]]$blacklist_ratio <- objList[[i]]$blacklist_region_fragments / objList[[i]]$peak_region_fragments
    filename <- paste0("./01_qc/QC_NS_TSS_density_",dataset.name,".pdf")
    pdf(filename,width=8,height=4)
    p1<- plot(density(na.omit(objList[[i]]$nucleosome_signal)),xlim=c(0,4))
    p2<- plot(density(na.omit(objList[[i]]$TSS.enrichment)),xlim=c(0,20))
    p3<- plot(density(na.omit(objList[[i]]$nCount_ATAC)),xlim=c(0,200000))
    p4<- plot(density(na.omit(objList[[i]]$nFeature_ATAC)),xlim=c(0,50000))
    p5<- plot(density(na.omit(objList[[i]]$pct_reads_in_peaks)),xlim=c(0,100))
    p6<- plot(density(na.omit(objList[[i]]$blacklist_ratio)),xlim=c(0,20000))
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    dev.off()
    }

for (i in seq_len(length(objList))) {
  # plot QC plot 
  pdf(file = paste("./01_qc/QC_before_Vlnplot_",sample[i],".pdf", sep = ""),width=10,height=6)
  qc<-VlnPlot(object = objList[[i]],
            features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'nucleosome_signal'),
            ncol = 3,
            pt.size = 0.01
          )
  print(qc)
  dev.off();
   }

#############################################
# 3. doubletfinder(ArchR) #
#############################################
library(ArchR)
addArchRThreads(threads = 1) 
addArchRGenome("mm10")
inputFiles <- c("/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR3/outs/fragments.tsv.gz",
  "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR7/outs/fragments.tsv.gz",
  "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/AR14/outs/fragments.tsv.gz",
  "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P4/outs/fragments.tsv.gz",
  "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P8/outs/fragments.tsv.gz",
  "/md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/HR_scATAC_data/P15/outs/fragments.tsv.gz")
names(inputFiles)<- sample
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "02_ArchR",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
projHeme2 <- filterDoublets(projHeme1)





filterDoublets_C4<- projHeme2$cellNames[grep("AR3_C4",projHeme2$cellNames)]
filterDoublets_C4<- gsub("AR3_C4#","",filterDoublets_C4)
filterDoublets_C5<- projHeme2$cellNames[grep("AR3_C5",projHeme2$cellNames)]
filterDoublets_C5<- gsub("AR3_C5#","",filterDoublets_C5)

##############################################
# 4. filter cells (low quality and doublets) #
##############################################

# filter out low quality cells
# To remove doublets,select different cutoff#####
  objList2<-c()
  objList2[[1]]<-subset(x=objList[[1]],
   subset = nCount_ATAC < 100000 &nCount_ATAC > 1000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 25 
    )
  print(objList2[[1]])

  objList2[[2]]<-subset(x=objList[[2]],
   subset = nCount_ATAC < 100000 &nCount_ATAC > 1000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    peak_region_fragments > 1000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 25 
    )
  print(objList2[[2]])

# Subset identified singlets
last_C4<-intersect(colnames(objList2[[1]]),filterDoublets_C4)
AR3_C4_last<- subset(objList2[[1]],cells=last_C4)
last_C5<-intersect(colnames(objList2[[2]]),filterDoublets_C5)
AR3_C5_last<- subset(objList2[[2]],cells=last_C5)
saveRDS(AR3_C4_last,"./03_all_celltype/AR3_C4_scATAC.rds")
saveRDS(AR3_C5_last,"./03_all_celltype/AR3_C5_scATAC.rds")

#####################################
# 5. intergrated object and analyze #
#####################################
combined <- merge(
  x = AR3_C4_last,
  y = AR3_C5_last,
  add.cell.ids = c("AR3_C4", "AR3_C5")
)
combined[["ATAC"]]
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
saveRDS(combined,"./03_all_celltype/combined_AR3_scATAC.rds")







