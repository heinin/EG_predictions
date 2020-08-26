# ======================================
# Comparing features of two sets of E-G predictions
# ======================================

# Work in progress

# Required columns:
# chr
# start [element start pos]
# end [element end pos]
# CellType
# TargetGene
# startTSS [gene TSS pos]

# If a cell type mapping file is provided, column names must match options --pred1name
# and --pred2name.

# ======================================
# Load libraries
# ======================================

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(source("/home/users/hnatri/Ernst_Roadmap/JuicerUtilities.R"))

# ======================================
# Parsing command line arguments
# ======================================

option.list <- list(
  make_option("--pred1", type="character", default="/home/groups/engreitz/Projects/Ernst_Roadmap/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz", help="First set of predictions"),
  make_option("--pred1name", type="character", default="ABC", help="Name of the first set of predictions"),
  make_option("--pred2", type="character", default="/home/groups/engreitz/Projects/Ernst_Roadmap/Liu2017-EnhancerGeneLinksMerged_collapse200bp.tsv.gz", help="Second set of predictions"),
  make_option("--pred2name", type="character", default="Roadmap", help="Name of the second set of predictions"),
  make_option("--promoteractivity", type="character", default="/home/groups/engreitz/Projects/Ernst_Roadmap/GeneTSSActivityQuantile.tsv", help="Promoter activity quantiles for defining expression"),
  make_option("--ABCkeeppromoters", type="logical", default=FALSE, help="Keep promoters for ABC?"),
  make_option("--celltypemap", type="character", default="/home/groups/engreitz/Projects/Ernst_Roadmap/RoadmapABCCellMapping.txt", help="Mapping of matching cell types"),
  make_option("--sharedcelltypes", type="logical", default=TRUE, help="Only compare shared cell types?"),
  make_option("--out.dir", type="character", default="/scratch/users/hnatri/Ernst_Roadmap", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option.list))

# ======================================
# Environment variables
# ======================================

out.dir <- opt$out.dir

# ======================================
# Data
# ======================================

# Enhancer predictions
pred1 <-  fread(opt$pred1)
pred1name <- opt$pred1name
pred2 <- fread(opt$pred2)
pred2name <- opt$pred2name
promoter.activity <- fread(opt$promoteractivity)

# Keeping ABC promoters?
if (pred1name=="ABC" & opt$ABCkeeppromoters==FALSE){
  pred1 <- subset(pred1, class != "promoter")
}

if (pred2name=="ABC" & opt$ABCkeeppromoters==FALSE){
  pred2 <- subset(pred1, class != "promoter")
}

# Only compare shared cell types?
if (opt$sharedcelltypes){
  celltype.mapping <- fread(opt$celltypemap)
  # Only keeping the cell types included in both sets of predictions
  celltype.mapping[,..pred1name]
  pred1 <- pred1[CellType %in% celltype.mapping[,get(pred1name)]]
  pred2 <- pred2[CellType %in% celltype.mapping[,get(pred2name)]]
  colnames(pred2) <- gsub("CellType", pred2name, colnames(pred2))
  
  # Matching cell type names across prediction sets by joining set 2 
  # predictions and mapped celltype names.
  # Setting keys for tables:
  setkeyv(pred2,pred2name)
  setkeyv(celltype.mapping,pred2name)
  pred2.matchingcelltypes <- pred2[celltype.mapping, nomatch=0]
  head(pred2.matchingcelltypes)
  pred2 <- pred2.matchingcelltypes
  colnames(pred2) <- gsub(pred1name, "CellType", colnames(pred2))
  colnames(pred1) <- gsub(pred1name, "CellType", colnames(pred1))
}

# Sanity check, number of shared celltypes
length(unique(pred1$CellType))
length(unique(pred2$CellType))

# Adding some columns to both prediction tables
pred.list <- list(pred1, pred2)
preds <- lapply(pred.list, function(dt) {
  dt$range <- with(dt, paste0(start, "-", end))
  dt$name <- with(dt, paste0(chr, ":", range))
  dt$gene.cell <- with(dt, paste0(TargetGene, ":", CellType))
  dt$enh.cell <- with(dt, paste0(name, ":", CellType))
  dt$enh.gene <- with(dt, paste0(name, ":", TargetGene))
  dt$enhancerbps <- with(dt, end-start)
  dt
})
pred1 <- preds[[1]]
pred2 <- preds[[2]]

# Finding genes that have promoter activity >= 0.4 quantile in each cell type
promoter.activity <- promoter.activity[, .SD, .SDcols = c("name", celltype.mapping[,get(pred1name)])]
promoter.activity.melt <- promoter.activity %>% 
  gather(celltype, expression, -name)
promoter.activity.melt$gene.cell <- with(promoter.activity.melt, paste0(name, ":", celltype))
promoter.activity.melt <- promoter.activity.melt[promoter.activity.melt$expression>=0.4,]

# Keeping expressed genes in each cell type
pred1 <- pred1[gene.cell %in% promoter.activity.melt$gene.cell]
dim(pred1)
pred2 <- pred2[gene.cell %in% promoter.activity.melt$gene.cell]
dim(pred2)

# Testing
set.seed(123)
pred1 <- pred1[sample(1:nrow(pred1), 20000, replace=FALSE) ,]
pred2 <- pred2[sample(1:nrow(pred2), 20000, replace=FALSE) ,]

# ======================================
# Summary statistics
# ======================================

# For both sets of predictions, summarizing features per cell type

pred.list <- list(pred1, pred2)
names(pred.list) <- c(pred1name, pred2name)
CellTypeSummaries <- lapply(pred.list, function(dt) {
  CellTypeSummary <- dt %>% group_by(CellType) %>% summarise(
    mean.e.degree = mean(as.data.table(table(name))$N),
    mean.g.degree.exclude0 = mean(as.data.table(table(TargetGene))$N),
    num.genes.with.gt.1.enh = length(unique(TargetGene)),
    num.connections = n(),
    num.enhancers = length(unique(name)),
    median.distance = as.numeric(median(abs(distance))),
    enh.sequence = sum(end-start),
    median.enhancerbps = as.numeric(median(enhancerbps)))
  as.data.table(CellTypeSummary)
})

names(CellTypeSummaries) <- c(pred1name, pred2name)
CellTypeSummaries.df <- merge(CellTypeSummaries[[1]], CellTypeSummaries[[2]], by="CellType", all=T, suffixes = c(paste0(".", pred1name), paste0(".", pred2name)))

# Saving to a file
file.name <- paste(pred1name, "CellTypeSummary.tsv", sep="_")
write.table(CellTypeSummaries[[pred1name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

file.name <- paste(pred2name, "CellTypeSummary.tsv", sep="_")
write.table(CellTypeSummaries[[pred2name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

file.name <- paste(pred1name, pred2name, "CellTypeSummary.tsv", sep="_")
write.table(CellTypeSummaries.df, file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

# For each gene in a given cell type, finding the total number of enhancers and unique basepairs
GeneCellTypeSummaries <- lapply(pred.list, function(dt) {
  GeneCellTypeSummary <- dt %>% group_by(gene.cell) %>% summarise(
    TargetGene=unique(TargetGene),
    CellType=unique(CellType),
    #gene.cell=gene.cell,
    TotalEnhancers=n(),
    TotalUniqueBases=sum(width(reduce(GRangesFromBed(data.frame(chr=chr, start=start, end=end))))))
  as.data.table(GeneCellTypeSummary)
})

names(GeneCellTypeSummaries) <- c(pred1name, pred2name)
GeneCellTypeSummaries.df <- merge(GeneCellTypeSummaries[[1]], GeneCellTypeSummaries[[2]], by="gene.cell", all=T, suffixes = c(paste0(".", pred1name), paste0(".", pred2name)))

# Saving to a file
file.name <- paste(pred1name, "GeneCellTypeSummary.tsv", sep="_")
write.table(GeneCellTypeSummaries[[pred1name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

file.name <- paste(pred2name, "GeneCellTypeSummary.tsv", sep="_")
write.table(GeneCellTypeSummaries[[pred2name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

#file.name <- paste(pred1name, pred2name, "GeneCellTypeSummary.tsv", sep="_")
#write.table(GeneCellTypeSummaries.df, file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

# For each gene across all cell types, finding the number of enhancers, unique bps,
# number of cell types with enhancer predictions
GeneSummaries <- lapply(pred.list, function(dt) {
  GeneSummary <- dt %>% group_by(TargetGene) %>% summarise(
    TotalEnhancers=n(),
    TotalUniqueBases=sum(width(reduce(GRangesFromBed(data.frame(chr=chr, start=start, end=end))))),
    CellTypesWithPrediction=length(unique(CellType)))
  as.data.table(GeneSummary)
})

names(GeneSummaries) <- c(pred1name, pred2name)
GeneSummaries.df <- merge(GeneSummaries[[1]], GeneSummaries[[2]], by="TargetGene", all=T, suffixes = c(paste0(".", pred1name), paste0(".", pred2name)))

# Writing to a file
file.name <- paste(pred1name, "GeneSummary.tsv", sep="_")
write.table(GeneSummaries[[pred1name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

file.name <- paste(pred2name, "GeneSummary.tsv", sep="_")
write.table(GeneSummaries[[pred2name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

#file.name <- paste(pred1name, pred2name, "GeneSummary.tsv", sep="_")
#write.table(GeneSummaries.df, file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

# For each set of predicitions, finding expressed genes without enhancers in each cell type
GeneSummariesZeroes <- list()

for (pred in c(pred1name, pred2name)){
  dt <- GeneCellTypeSummaries[[pred]]
  genes.wout.enhancers <- matrix(NA, ncol = 2)
  colnames(genes.wout.enhancers) <- c("TargetGene", "CellType")
  
  genes.wout.enhancers <- lapply(unique(dt$CellType), function(celltype)
    setdiff(promoter.activity.melt[celltype==celltype, "name"], dt[CellType==celltype]$TargetGene)
  )
  
  names(genes.wout.enhancers) <- unique(dt$CellType)
  # Collapsing to a table
  genes.wout.enhancers.df1 <- lapply(genes.wout.enhancers, as.data.frame, stringsAsFactors = FALSE)
  genes.wout.enhancers.df2 <- bind_rows(genes.wout.enhancers.df1, .id = "Name")
  colnames(genes.wout.enhancers.df2) <- c("CellType", "TargetGene")
  genes.wout.enhancers.df2$gene.cell <- with(genes.wout.enhancers.df2, paste0(TargetGene, ":", CellType))
  # Adding zeros for genes that don't have any enhancers
  #genes.wout.enhancers.df2$enh.sequence <- 0
  genes.wout.enhancers.df2$TotalEnhancers <- 0
  genes.wout.enhancers.df2$TotalUniqueBases <- 0
  
  # Adding to the list
  GeneSummariesZeroes[[pred]] <- genes.wout.enhancers.df2
  
}

names(GeneSummariesZeroes) <- c(pred1name, pred2name)

# Adding genes with zero enhancers to the summary tables
GeneCellTypeSummariesWith0 <- list()
for (pred in c(pred1name, pred2name)){
  GeneCellTypeSummaryWith0 <- rbind(GeneCellTypeSummaries[[pred]], GeneSummariesZeroes[[pred]])
  GeneCellTypeSummariesWith0[[pred]] <- GeneCellTypeSummaryWith0
}

names(GeneCellTypeSummariesWith0) <- c(pred1name, pred2name)

# Writing to a file
file.name <- paste(pred1name, "GeneCelltypeSummaryWith0.tsv", sep="_")
write.table(GeneCellTypeSummariesWith0[[pred1name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

file.name <- paste(pred2name, "GeneCelltypeSummaryWith0.tsv", sep="_")
write.table(GeneCellTypeSummariesWith0[[pred2name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

# For each gene, finding the number of cell types where the gene is expressed 
# but there are no predicted enhancers
GeneSummariesWith0 <- list()
for (pred in c(pred1name, pred2name)){
  NCellTypesW0Enhancers <- as.data.frame(table(GeneSummariesZeroes[[pred]]$TargetGene))
  colnames(NCellTypesW0Enhancers) <- c("TargetGene", "NExpCelltypesWith0Enhancers")
  
  # Adding to the gene summary table
  GeneSummaryWith0 <- merge(GeneSummaries[[pred]], NCellTypesW0Enhancers, by="TargetGene", all=T)
  GeneSummaryWith0[is.na(GeneSummaryWith0)] <- 0
  #GeneSummaryWith0$NExpCelltypesWith0Enhancers[is.na(GeneSummaryWith0$NExpCelltypesWith0Enhancers)] <- 0
  GeneSummaryWith0$CellTypesWithPredictionIncl0 <- GeneSummaryWith0$CellTypesWithPrediction+GeneSummaryWith0$NExpCelltypesWith0Enhancers
  
  # Adding to the list
  GeneSummariesWith0[[pred]] <- GeneSummaryWith0
}

names(GeneSummariesWith0) <- c(pred1name, pred2name)

# Writing to a file
file.name <- paste(pred1name, "GeneSummaryWith0.tsv", sep="_")
write.table(GeneSummariesWith0[[pred1name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

file.name <- paste(pred2name, "GeneSummaryWith0.tsv", sep="_")
write.table(GeneSummariesWith0[[pred2name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

# Adding expressed genes with no enhancers to the celltype summary
CellTypeSummariesWith0 <- lapply(c(pred1name, pred2name), function(pred) {
  gene0 <- GeneSummariesZeroes[[pred]] %>% group_by(CellType) %>%
    summarise(NumGenesWith0Enh = sum(n()))
  
  CellTypeSummaryWith0 <- merge(CellTypeSummaries[[pred]], as.data.table(gene0))
  CellTypeSummaryWith0$mean.g.degree.include0 <- with(CellTypeSummaryWith0, mean.g.degree.exclude0 * (num.genes.with.gt.1.enh/(num.genes.with.gt.1.enh + NumGenesWith0Enh)))
  
  as.data.table(CellTypeSummaryWith0)
})

names(CellTypeSummariesWith0) <- c(pred1name, pred2name)

# Writing to a file
file.name <- paste(pred1name, "CelltypeSummaryWith0.tsv", sep="_")
write.table(CellTypeSummariesWith0[[pred1name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

file.name <- paste(pred2name, "CelltypeSummaryWith0.tsv", sep="_")
write.table(CellTypeSummariesWith0[[pred2name]], file.path(out.dir, file.name), sep = "\t", quote = F, col.names = T, row.names = F)

# ======================================
# Plotting features
# ======================================

# All enhancers, distance to TSS in a given cell types.

TSSplotData <- lapply(names(pred.list), function(pred) {
  dt <- pred.list[[pred]]
  absdistances <- as.data.frame(abs(dt[,distance]))
  absdistances$group <- pred
  colnames(absdistances) <- c("absdistance", "group")
  
  absdistances
})

TSSplotData.df <- as.data.frame(do.call(rbind, TSSplotData))

# Cumulative density 
distTSSplot <- ggplot(TSSplotData.df, aes(x=absdistance, color=group)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Distance to TSS") + 
  ylab("Cumulative density")

# Saving as a PDF
file.name <- paste(pred1name, pred2name, "distancetoTSS.pdf", sep="_")
ggsave(file.path(out.dir, file.name), 
       distTSSplot,
       width = 6, height = 5)

# Number of enhancers in a given cell type

numEnhancersPlotData <- lapply(names(CellTypeSummariesWith0), function(pred) {
  dt <- CellTypeSummariesWith0[[pred]]
  numEnhancers <- as.data.frame(dt[,num.enhancers])
  numEnhancers$group <- pred
  colnames(numEnhancers) <- c("numEnhancers", "group")
  
  numEnhancers
})

numEnhancersPlotData.df <- as.data.frame(do.call(rbind, numEnhancersPlotData))

# Cumulative density 
numEnhancersPlot <- ggplot(numEnhancersPlotData.df, aes(x=numEnhancers, color=group)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Number of enhancers in a given cell type") + 
  ylab("Cumulative density")

# Saving as a PDF
file.name <- paste(pred1name, pred2name, "enhancersPerCelltype.pdf", sep="_")
ggsave(file.path(out.dir, file.name), 
       numEnhancersPlot,
       width = 6, height = 5)

# Number of enhancer-gene connections in a given cell type

numEGconnectionsPlotData <- lapply(names(CellTypeSummariesWith0), function(pred) {
  dt <- CellTypeSummariesWith0[[pred]]
  numConnections <- as.data.frame(dt[,num.connections])
  numConnections$group <- pred
  colnames(numConnections) <- c("numConnections", "group")
  
  numConnections
})

numEGconnectionsPlotData.df <- as.data.frame(do.call(rbind, numEGconnectionsPlotData))

# Cumulative density 
numEGconnectionsPlot <- ggplot(numEGconnectionsPlotData.df, aes(x=numConnections, color=group)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Number of connections ina given cell type") + 
  ylab("Cumulative density")

# Saving as a PDF
file.name <- paste(pred1name, pred2name, "connectionsPerCelltype.pdf", sep="_")
ggsave(file.path(out.dir, file.name), 
       numEnhancersPlot,
       width = 6, height = 5)

# Mean number of genes regulated by each enhancer in a given cell type

eDegreePlotData <- lapply(names(CellTypeSummariesWith0), function(pred) {
  dt <- CellTypeSummariesWith0[[pred]]
  eDegree <- as.data.frame(dt[,mean.e.degree])
  eDegree$group <- pred
  colnames(eDegree) <- c("mean.e.degree", "group")
  
  eDegree
})

eDegreePlotData.df <- as.data.frame(do.call(rbind, eDegreePlotData))

# Cumulative density 
eDegreePlot <- ggplot(eDegreePlotData.df, aes(x=mean.e.degree, color=group)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Mean number of genes per enhancer in  a given cell type") + 
  ylab("Cumulative density")

# Saving as a PDF
file.name <- paste(pred1name, pred2name, "GperE.PerCelltype.pdf", sep="_")
ggsave(file.path(out.dir, file.name), 
       eDegreePlot,
       width = 6, height = 5)

# Mean number of enhancers connected to each gene in a given cell type, including genes
# with zero enhancers

gDegreePlotData <- lapply(names(CellTypeSummariesWith0), function(pred) {
  dt <- CellTypeSummariesWith0[[pred]]
  gDegree <- as.data.frame(dt[,mean.g.degree.include0])
  gDegree$group <- pred
  colnames(gDegree) <- c("mean.g.degree.include0", "group")
  
  gDegree
})

gDegreePlotData.df <- as.data.frame(do.call(rbind, gDegreePlotData))

# Cumulative density 
gDegreePlot <- ggplot(gDegreePlotData.df, aes(x=mean.g.degree.include0, color=group)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Mean number of enhancers per gene in a given cell type") + 
  ylab("Cumulative density")

# Saving as a PDF
file.name <- paste(pred1name, pred2name, "EperG.PerCelltype.Incl0.pdf", sep="_")
ggsave(file.path(out.dir, file.name), 
       gDegreePlot,
       width = 6, height = 5)

# What is the total number of unique bps included in enhancers (for any gene) 
# across all shared cell types? Including genes with zero enhancers.

geneEnhancerSummaries <- lapply(c(pred1name, pred2name), function(pred) {
  geneEnhancerSummary <- GeneCellTypeSummariesWith0[[pred]] %>% group_by(TargetGene) %>%
    summarise(#Nenhancers=sum(TotalEnhancers),
      UniqueBases=sum(TotalUniqueBases))
  
  as.data.table(geneEnhancerSummary)
})

names(geneEnhancerSummaries) <- c(pred1name, pred2name)

geneEnhancerSummary <- merge(geneEnhancerSummaries[[1]], geneEnhancerSummaries[[2]], by="TargetGene", all=T, suffixes = c(paste0(".", pred1name), paste0(".", pred2name)))
geneEnhancerSummary.melt <- geneEnhancerSummary %>% gather(pred, TotalUniqueBases, -TargetGene)
geneEnhancerSummary.melt$pred <- gsub("UniqueBases.", "", geneEnhancerSummary.melt$pred)

# Cumulative density 
Ebps <- ggplot(geneEnhancerSummary.melt, aes(x=TotalUniqueBases, color=pred)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Unique enhancer bps per gene across shared cell types") + 
  ylab("Cumulative density")

# Saving as a PDF
file.name <- paste(pred1name, pred2name, "uniqueEnhancerbpsPerGeneAllCelltypes.pdf", sep="_")
ggsave(file.path(out.dir, file.name), 
       Ebps,
       width = 6, height = 5)

# What the total number of unique bps included in enhancers (for any gene) 
# in a given cell type? Including genes with zero enhancers.

geneCelltypeEnhancerSummaries <- lapply(c(pred1name, pred2name), function(pred) {
  geneCelltypeEnhancerSummary <- GeneCellTypeSummariesWith0[[pred]] %>% group_by(gene.cell) %>%
    summarise(#Nenhancers=sum(TotalEnhancers),
      UniqueBases=sum(TotalUniqueBases),
      gene.cell=gene.cell)
  
  as.data.table(geneCelltypeEnhancerSummary)
})

names(geneCelltypeEnhancerSummaries) <- c(pred1name, pred2name)

geneCelltypeEnhancerSummary <- merge(geneCelltypeEnhancerSummaries[[1]], geneCelltypeEnhancerSummaries[[2]], by="gene.cell", all=T, suffixes = c(paste0(".", pred1name), paste0(".", pred2name)))
geneCelltypeEnhancerSummary.melt <- geneCelltypeEnhancerSummary %>% gather(pred, UniqueBases, -gene.cell)
geneCelltypeEnhancerSummary.melt$pred <- gsub("UniqueBases.", "", geneCelltypeEnhancerSummary.melt$pred)

# Cumulative density 
Ebps_each_celltype <- ggplot(geneCelltypeEnhancerSummary.melt, aes(x=UniqueBases, color=pred)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Unique enhancer bps per gene in a given cell type") + 
  ylab("Cumulative density")

# Saving as a PDF
file.name <- paste(pred1name, pred2name, "uniqueEnhancerbpsPerGeneGivenCelltype.pdf", sep="_")
ggsave(file.path(out.dir, file.name), 
       Ebps_each_celltype,
       width = 6, height = 5)

# N of enhancers per gene per cell type where the gene is expressed,
# including expressed genes with zero enhancers

geneCelltypeEnhancerExpressedSummaries <- lapply(c(pred1name, pred2name), function(pred) {
  geneCelltypeEnhancerExpressedSummary <- GeneCellTypeSummariesWith0[[pred]] %>% group_by(gene.cell) %>%
    summarise(
      NenhancersPerExpressed=sum(TotalEnhancers)/length(unique(CellType)))
  
  as.data.table(geneCelltypeEnhancerExpressedSummary)
})

names(geneCelltypeEnhancerExpressedSummaries) <- c(pred1name, pred2name)

geneCelltypeEnhancerExpressedSummary <- merge(geneCelltypeEnhancerExpressedSummaries[[1]], geneCelltypeEnhancerExpressedSummaries[[2]], by="gene.cell", all=T, suffixes = c(paste0(".", pred1name), paste0(".", pred2name)))
geneCelltypeEnhancerExpressedSummary.melt <- geneCelltypeEnhancerExpressedSummary %>% gather(pred, NenhancersPerExpressed, -gene.cell)
geneCelltypeEnhancerExpressedSummary.melt$pred <- gsub("NenhancersPerExpressed.", "", geneCelltypeEnhancerExpressedSummary.melt$pred)

# Cumulative density 
enhancerGeneCelltype <- ggplot(geneCelltypeEnhancerExpressedSummary.melt, aes(x=NenhancersPerExpressed, color=pred)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Enhancers per gene / # cell types where the gene is expressed") + 
  ylab("Cumulative density")

# Saving as a PDF
file.name <- paste(pred1name, pred2name, "EnhancersPerGenePerCelltypeExpWith0.pdf", sep="_")
ggsave(file.path(out.dir, file.name), 
       Ebps_each_celltype,
       width = 6, height = 5)

# ======================================
# Finding overlapping and unique elements
# ======================================

# N of enhancers
length(unique(pred1$name))
length(unique(pred2$name))

# To GRanges
pred1.gr <- GRanges(seqnames=pred1$chr,
                    ranges=IRanges(start=pred1$start, end=pred1$end))

pred2.gr <- GRanges(seqnames=pred2$chr,
                    ranges=IRanges(start=pred2$start, end=pred2$end))

# Overlapping elements
pred1.overlapping.pred2.hits <- as.data.table(findOverlaps(pred1.gr, pred2.gr, minoverlap=1))
pred1.overlapping.pred2 <- pred1[unique(pred1.overlapping.pred2.hits$queryHits)]
pred2.overlapping.pred1 <- pred2[unique(pred1.overlapping.pred2.hits$subjectHits)]
dim(pred1.overlapping.pred2)
dim(pred2.overlapping.pred1)
length(unique(pred1.overlapping.pred2$name))
length(unique(pred2.overlapping.pred1$name))

# Unique to the first set of predictions
pred1.unique <- pred1[-(unique(pred1.overlapping.pred2.hits$queryHits))]
dim(pred1.unique)
length(unique(pred1.unique$name))

# Unique to the second set of predictions
pred2.unique <- pred2[-(unique(pred1.overlapping.pred2.hits$subjectHits))]
dim(pred2.unique)
length(unique(pred2.unique$name))

pred1Enhancers <- unique(c(pred1.unique$name, pred1.overlapping.pred2$name))
pred2Enhancers <- unique(c(pred2.unique$name, pred1.overlapping.pred2$name))

# Proportion of shared elements
length(intersect(pred1Enhancers, pred2Enhancers))/length(unique(c(pred1Enhancers, pred2Enhancers)))

# Plotting
filename <- paste(pred1name, pred2name, "elements_Venn.tiff", sep="_")
venn.diagram(
  x = list(pred1Enhancers, pred2Enhancers),
  category.names = c(pred1name , pred2name),
  filename = file.path(out.dir, filename),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  #lwd = 2,
  #lty = 'blank',
  #fill = 'blank',
  
  # Numbers
  cex = .6,
  #fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
  #rotation = 1
)

# ======================================
# # Overlapping and unique enhancer-gene combos across all cell types
# ======================================

# Numbers of unique enhancer-gene combos
pred1$enh.gene <- paste(pred1$name, pred1$TargetGene, sep="_")
pred2$enh.gene <- paste(pred2$name, pred2$TargetGene, sep="_")
length(unique(pred1$enh.gene))
length(unique(pred2$enh.gene))

# All unique genes
genes <- unique(c(pred1$TargetGene, pred2$TargetGene))
length(genes)

# TODO: future_lapply parallelization doesn't work on the cluster?
# Run in parallel
plan(multisession)
options(future.globals.maxSize = +Inf)
# Start the clock
ptm <- proc.time()
# For every gene, find overlapping and unique elements
SharedUniqueElementsGene <- future_lapply(genes, function(TargetGene){
  pred1.gene <- pred1[TargetGene==TargetGene]
  pred2.gene <- pred2[TargetGene==TargetGene]
  
  pred1.gene.gr <- GRanges(seqnames=pred1.gene$chr,
                           ranges=IRanges(start=pred1.gene$start, end=pred1.gene$end))
  
  pred2.gene.gr <- GRanges(seqnames=pred2.gene$chr,
                           ranges=IRanges(start=pred2.gene$start, end=pred2.gene$end))
  
  # Overlapping elements
  pred1.overlapping.pred2.hits <- as.data.table(findOverlaps(pred1.gene.gr, pred2.gene.gr, minoverlap=1))
  
  # If there are no overlapping elements...
  if (nrow(pred1.overlapping.pred2.hits)==0){
    pred1.unique <- pred1.gene
    pred1.unique$group <- paste0(pred1name, ".unique")
    pred2.unique <- pred2.gene
    pred2.unique$group <- paste0(pred2name, ".unique")
    
    merged_groups <- dplyr::bind_rows(pred1.unique, pred2.unique)
    
  } else {
    pred1.overlapping.pred2 <- pred1.gene[unique(pred1.overlapping.pred2.hits$queryHits)]
    pred2.overlapping.pred1 <- pred2.gene[unique(pred1.overlapping.pred2.hits$subjectHits)]
    
    pred1.overlapping.pred2$group <- paste0(pred1name, ".overlapping.", pred2name)
    pred2.overlapping.pred1$group <- paste0(pred2name, ".overlapping.", pred1name)
    
    # Unique to pred1
    pred1.unique <- pred1.gene[-(unique(pred1.overlapping.pred2.hits$queryHits))]
    dim(pred1.unique)
    if (nrow(pred1.unique)>0){
      pred1.unique$group <- paste0(pred1name, ".unique")
    }
    
    # Unique to pred2
    pred2.unique <- pred2.gene[-(unique(pred1.overlapping.pred2.hits$subjectHits))]
    dim(pred2.unique)
    if (nrow(pred2.unique)>0){
      pred2.unique$group <- paste0(pred2name, ".unique")
    }
    
    merged_groups <- dplyr::bind_rows(pred1.overlapping.pred2, pred1.unique, pred2.unique)
    
  }
  
  merged_groups
  
})

# Stop the clock
proc.time() - ptm

names(SharedUniqueElementsGene) <- genes
# Collapsing to a table
SharedUniqueElementsGene_df <- as.data.frame(do.call(rbind, SharedUniqueElementsGene))
head(SharedUniqueElementsGene_df)
dim(SharedUniqueElementsGene_df)

# Writing to a file
filename <- paste(pred1name, pred2name, "_shared_unique_TargetGene.tsv", sep="_")
write.table(SharedUniqueElementsGene_df, file.path(out.dir, filename), sep = "\t", row.names = F)

# Cumulative density plot of distances to the TSS for each groups
distTSSplotSharedUnique <- ggplot(SharedUniqueElementsGene_df, aes(x=abs(distance), color=group)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Distance to TSS") + 
  ylab("Cumulative density")

# Saving as a PDF
filename <- paste(pred1name, pred2name, "shared_unique_TargetGene_distTSSplot.pdf", sep="_")
ggsave(file.path(out.dir, filename), 
       distTSSplot,
       width = 6, height = 5)

# TODO: plot additional features
# N of cell types where the gene is expressed?

# Venn diagram of overlapping enhancer-gene-celltype combos
SharedUniqueElementsGene_df$longname <- paste(SharedUniqueElementsGene_df$chr, SharedUniqueElementsGene_df$start, SharedUniqueElementsGene_df$end, SharedUniqueElementsGene_df$TargetGene, sep="_")

pred1.links <- SharedUniqueElementsGene_df[which(SharedUniqueElementsGene_df$group %in% c("pred1.unique", "pred1.overlapping.pred2")),]$longname
pred2.links <- SharedUniqueElementsGene_df[which(SharedUniqueElementsGene_df$group %in% c("pred2.unique", "pred1.overlapping.pred2")),]$longname

# Proportion of shared enhancer-gene-celltype combos
length(intersect(pred1.links, pred2.links))/length(unique(c(pred1.links, pred2.links)))

filename <- paste(pred1name, pred2name, "EG_Venn.tiff", sep="_")
venn.diagram(
  x = list(pred1.links, pred2.links),
  category.names = c(pred1name , pred2name),
  filename = file.path(out.dir, filename),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  #lwd = 2,
  #lty = 'blank',
  #fill = 'blank',
  
  # Numbers
  cex = .6,
  #fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(27, -27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  inverted=TRUE,
  reverse=TRUE
  #rotation = 1
)

# ======================================
# Finding overlapping and unique element-gene-celltype combinations
# ======================================

# Numbers of unique enhancer-gene-celltype combos
pred1$name.gene.cell <- paste(pred1$name, pred1$gene.cell, sep="_")
pred2$name.gene.cell <- paste(pred2$name, pred2$gene.cell, sep="_")
length(unique(pred1$name.gene.cell))
length(unique(pred2$name.gene.cell))

# Numbers of gene-celltype combinations
length(unique(pred1$gene.cell))
length(unique(pred2$gene.cell))

# All unique gene-celltype combos
genes_celltypes <- unique(c(pred1$gene.cell, pred2$gene.cell))
length(genes_celltypes)

# Run in parallel
plan(multisession)
options(future.globals.maxSize = +Inf)
# Start the clock
ptm <- proc.time()
# For every gene-celltype combination, find overlapping and unique elements
SharedUniqueElementsGeneCell <- future_lapply(genes_celltypes, function(GeneCelltype){
  pred1.gene.cell <- pred1[gene.cell==GeneCelltype]
  pred2.gene.cell <- pred2[gene.cell==GeneCelltype]
  
  pred1.gene.cell.gr <- GRanges(seqnames=pred1.gene.cell$chr,
                                ranges=IRanges(start=pred1.gene.cell$start, end=pred1.gene.cell$end))
  
  pred2.gene.cell.gr <- GRanges(seqnames=pred2.gene.cell$chr,
                                ranges=IRanges(start=pred2.gene.cell$start, end=pred2.gene.cell$end))
  
  # Overlapping elements
  pred1.overlapping.pred2.hits <- as.data.table(findOverlaps(pred1.gene.cell.gr, pred2.gene.cell.gr, minoverlap=1))
  
  # If there are no overlapping elements...
  if (nrow(pred1.overlapping.pred2.hits)==0){
    pred1.unique <- pred1.gene.cell
    pred1.unique$group <- paste0(pred1name, ".unique")
    pred2.unique <- pred2.gene.cell
    pred2.unique$group <- paste0(pred2name, ".unique")
    
    merged_groups <- dplyr::bind_rows(pred1.unique, pred2.unique)
    
  } else {
    pred1.overlapping.pred2 <- pred1.gene.cell[unique(pred1.overlapping.pred2.hits$queryHits)]
    pred2.overlapping.pred1 <- pred2.gene.cell[unique(pred1.overlapping.pred2.hits$subjectHits)]
    
    pred1.overlapping.pred2$group <- paste0(pred1name, ".overlapping.", pred2name)
    pred2.overlapping.pred1$group <- paste0(pred2name, ".overlapping.", pred1name)
    
    # Unique to pred1
    pred1.unique <- pred1.gene.cell[-(unique(pred1.overlapping.pred2.hits$queryHits))]
    dim(pred1.unique)
    if (nrow(pred1.unique)>0){
      pred1.unique$group <- paste0(pred1name, ".unique")
    }
    
    # Unique to pred2
    pred2.unique <- pred2.gene.cell[-(unique(pred1.overlapping.pred2.hits$subjectHits))]
    dim(pred2.unique)
    if (nrow(pred2.unique)>0){
      pred2.unique$group <- paste0(pred2name, ".unique")
    }
    
    merged_groups <- dplyr::bind_rows(pred1.overlapping.pred2, pred1.unique, pred2.unique)
    
  }
  
  merged_groups
  
})

# Stop the clock
proc.time() - ptm

names(SharedUniqueElementsGeneCell) <- genes_celltypes
# Collapsing to a table
SharedUniqueElementsGeneCell_df <- as.data.frame(do.call(rbind, SharedUniqueElementsGeneCell))
head(SharedUniqueElementsGeneCell_df)
dim(SharedUniqueElementsGeneCell_df)

# Writing to a file
filename <- paste(pred1name, pred2name, "_shared_unique_TargetGeneCellType.tsv", sep="_")
write.table(SharedUniqueElementsGeneCell_df, file.path(out.dir, filename), sep = "\t", row.names = F)

# Cumulative density plot of distances to the TSS for each groups
distTSSplotSharedUnique <- ggplot(SharedUniqueElementsGeneCell_df, aes(x=abs(distance), color=group)) + 
  stat_ecdf() + 
  theme_bw() +
  xlab("Distance to TSS") + 
  ylab("Cumulative density")

# Saving as a PDF
filename <- paste(pred1name, pred2name, "shared_unique_TargetGeneCellType_distTSSplot.pdf", sep="_")
ggsave(file.path(out.dir, filename), 
       distTSSplot,
       width = 6, height = 5)

# TODO: plot additional features
# N of cell types where the gene is expressed?

# Venn diagram of overlapping enhancer-gene-celltype combos
SharedUniqueElementsGeneCell_df$longname <- paste(SharedUniqueElementsGeneCell_df$chr, SharedUniqueElementsGeneCell_df$start, SharedUniqueElementsGeneCell_df$end, SharedUniqueElementsGeneCell_df$TargetGene, SharedUniqueElementsGeneCell_df$CellType, sep="_")

pred1.links <- SharedUniqueElementsGeneCell_df[which(SharedUniqueElementsGeneCell_df$group %in% c("pred1.unique", "pred1.overlapping.pred2")),]$longname
pred2.links <- SharedUniqueElementsGeneCell_df[which(SharedUniqueElementsGeneCell_df$group %in% c("pred2.unique", "pred1.overlapping.pred2")),]$longname

# Proportion of shared enhancer-gene-celltype combos
length(intersect(pred1.links, pred2.links))/length(unique(c(pred1.links, pred2.links)))

filename <- paste(pred1name, pred2name, "EGcelltype_Venn.tiff", sep="_")
venn.diagram(
  x = list(pred1.links, pred2.links),
  category.names = c(pred1name , pred2name),
  filename = file.path(out.dir, filename),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  #lwd = 2,
  #lty = 'blank',
  #fill = 'blank',
  
  # Numbers
  cex = .6,
  #fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(27, -27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  inverted=TRUE,
  reverse=TRUE
  #rotation = 1
)
