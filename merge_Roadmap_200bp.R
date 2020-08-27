# Enhancer predictions
roadmap.pred <- fread("/home/groups/engreitz/Projects/Ernst_Roadmap/Liu2017-EnhancerGeneLinksMerged.tsv.gz")

# Adding some columns
roadmap.pred$cell.target <- with(roadmap.pred, paste0(CellType, ";", GeneSymbol))
roadmap.pred$cell.target.TSS <- with(roadmap.pred, paste0(cell.target, ";", startTSS))
roadmap.pred <- roadmap.pred[,!c("strand")]

# Merging neighboring ranges, datatable to df to granges:
roadmap.pred.df <- as.data.frame(roadmap.pred)
roadmap.pred.gr <- makeGRangesFromDataFrame(roadmap.pred.df,
                                            seqnames.field="chrElement",
                                            start.field="startElement",
                                            end.field="endElement",
                                            keep.extra.columns=TRUE,
                                            ignore.strand=TRUE)

merged.gr <- unlist(reduce(split(roadmap.pred.gr, roadmap.pred.gr$cell.target.TSS)))
merged.gr$cell.target.TSS <- names(merged.gr)

# Replacing rownames
names(merged.gr) <- seq(1:length(names(merged.gr)))

# To a data.table
merged.dt <- data.table(as.data.frame(merged.gr))

# Splitting celltype-gene-tss, adding columns
merged.dt[, c("CellType", "TargetGene", "startTSS") := tstrsplit(cell.target.TSS, ";", fixed=TRUE)]
merged.dt$startTSS <- as.numeric(merged.dt$startTSS)
merged.dt$range <- with(merged.dt, paste0(start, "-", end))
merged.dt$name <- with(merged.dt, paste0(seqnames, ":", range))
merged.dt$distance <- with(merged.dt, startTSS-start)
merged.dt$absdistance <- with(merged.dt, abs(distance))
merged.dt$gene.cell <- with(merged.dt, paste0(TargetGene, ":", CellType))
merged.dt$enh.cell <- with(merged.dt, paste0(name, ":", CellType))
merged.dt$enhancerbps <- with(merged.dt, end-start)

colnames(merged.dt) <- gsub("seqnames", "chr", colnames(merged.dt))

write.table(merged.dt, "/home/groups/engreitz/Projects/Ernst_Roadmap/Liu2017-EnhancerGeneLinksMerged_collapse200bp.tsv", sep="\t", row.names=F, col.names=T, quote=F)
