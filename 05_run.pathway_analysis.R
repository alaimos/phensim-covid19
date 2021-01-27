library("org.Hs.eg.db")
library(readr)
library(dplyr)
library(HGNChelper)
library(openxlsx)
library(SPIA)
library(GOfuncR)
library(ReactomePA)
library(qvalue)
library(ggplot2)

degs.base  <- "data/LFCS/virus_blanco_melo/%s.csv"
cell.lines <- c("A549-ACE2_MOI_0.2", "A549-ACE2_MOI_2", "A549_MOI_0.2", "A549_MOI_2", "Calu-3", "NHBE")
p.value.th <- 0.05

### Change JAVA path as needed (Java 13 or higher)
run.mithril <- function (de.df, id.all) {
    tmp.input.file  <- tempfile()
    tmp.output.file <- tempfile()
    df.mith <- rbind(de.df, data.frame(id=id.all[!(id.all %in% de.df$id)], logFC=0))
    write.table(df.mith, tmp.input.file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    ### Change JAVA path as needed (Java 13 or higher)
    MITHRIL.COMMAND     <- "C:\\Java\\bin\\java -jar \"%s/MITHrIL2.jar\" mithril -verbose -m -i %s -enrichment-evidence-type STRONG -seed 3214 -o %s"
    run.command <- sprintf(MITHRIL.COMMAND, getwd(), tmp.input.file, tmp.output.file)
    ecode <- system(run.command, intern = FALSE, wait = TRUE)
    if (ecode != 0) {
        cat("An error occurred in MITHrIL. Exit code: ", ecode, "\n")
        return (NULL)
    }
    if (!file.exists(tmp.output.file)) {
        cat("An error occurred in MITHrIL. No output produced\n")
        return (NULL)
    }
    if (file.exists(tmp.input.file)) unlink(tmp.input.file)
    mith.output <- read_delim(tmp.output.file, "\t", escape_double = FALSE, na = "NA", trim_ws = TRUE)
    mith.output <- na.omit(mith.output)
    class(mith.output) <- "data.frame"
    colnames(mith.output) <- c("Pathway.Id", "Pathway.Name","Raw.Accumulator","Impact.Factor","Probability.Pi","Probability.Network","Corrected.Accumulator","pValue","Adjusted.pValue")
    mith.output$Pathway.Name <- gsub("\\s+\\-\\s+Enriched", "", mith.output$Pathway.Name, perl = TRUE)
    if (file.exists(tmp.output.file)) unlink(tmp.output.file)
    return (mith.output)
}

if (!dir.exists("results/PathwayAnalysis")) dir.create("results/PathwayAnalysis", recursive = TRUE)

results <- setNames(vector(mode = "list", length = length(cell.lines)), cell.lines)

# ### Uncomment this line if you wish to build plots using pre-computed results
# results <- readRDS("results/pathway.analysis.result.rds")

for (i in cell.lines) {
    print(i)
    lfcs <- read.csv(sprintf(degs.base, i), row.names=1, sep=";")
    corr.symb <- checkGeneSymbols(rownames(lfcs), unmapped.as.na = FALSE, species = "human")
    map <- setNames(strsplit(corr.symb$Suggested.Symbol, split = "\\s+///\\s+", perl = TRUE), corr.symb$x)
    df.map <- na.omit(data.frame(name=rep(names(map), sapply(map, length)), name_new=unname(unlist(map))))
    map <- mapIds(org.Hs.eg.db, df.map$name_new, 'ENTREZID', 'SYMBOL', multiVals = "list")
    df.map <- df.map %>% inner_join(na.omit(data.frame(name=rep(names(map), sapply(map, length)), id=unname(unlist(map)))), by = c("name_new"="name"))
    df.map <- df.map[,c("name", "id")]
    lfcs.map <- data.frame(name=rownames(lfcs), lfcs)
    lfcs.map <- lfcs.map %>% inner_join(df.map, by="name")
    id.all   <- unique(lfcs.map$id)
    nm.all   <- unique(lfcs.map$name)
    lfcs.map <- lfcs.map[lfcs.map$adj.P.Val < p.value.th & abs(lfcs.map$logFC) > 0.6,]
    de.nm    <- unique(lfcs.map$name)
    de.df    <- lfcs.map %>% dplyr::select(id, logFC) %>% group_by(id) %>% summarise_all(mean)
    de.vec   <- setNames(de.df$logFC, de.df$id)
    input.hy <- data.frame(gene_ids=de.nm, is_candidate=1)
    ### Comment from these lines if you wish build plots using pre-computed results
    res.gsea <- suppressWarnings(gsePathway(geneList = sort(de.vec, decreasing = TRUE), pvalueCutoff = 1))
    res.reac <- as.data.frame(res.gsea)
    res.hy   <- go_enrich(input.hy, n_randset=100)
    res.mith <- run.mithril(de.df = de.df, id.all = id.all)
    res.spia <- spia(de=de.vec, all = id.all)
    results[[i]] <- list(
        reactome=res.reac,
        go=res.hy,
        mithril=res.mith,
        spia=res.spia
    )
    save(i, results, file = "tmp.output.RData")
    ### The section to comment ends here!!
    x <- enrichPathway(gene=names(de.vec),pvalueCutoff=0.05, readable=T)
    tmp <- results[[i]]$mithril
    colnames(tmp)[1] <- "ID"
    colnames(tmp)[2] <- "Description"
    tmp$GeneRatio <- NA
    tmp$BgRatio <- NA
    colnames(tmp)[8] <- "pvalue"
    colnames(tmp)[9] <- "p.adjust"
    colnames(tmp)[7] <- "Accumulator"
    colnames(tmp)[4] <- "Impact"
    tmp$qvalue <- qvalue(tmp$pvalue)$qvalues
    tmp$geneID <- NA
    tmp <- tmp[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Accumulator", "Impact")]
    tmp <- tmp[tmp$Description != "Metabolic pathways",]
    tmp$Impact <- round(tmp$Impact)
    tmp1 <- tmp[tmp$p.adjust < 0.05,]
    tmp1 <- tmp1[order(tmp1$Accumulator),]
    nneg <- min(12,nrow(tmp1[tmp1$Accumulator < 0,]))
    npos <- min(12,nrow(tmp1[tmp1$Accumulator > 0,]))
    tmp2 <- tmp
    tmp <- rbind(tmp1[1:nneg,], tmp1[(nrow(tmp1)-npos):nrow(tmp1),])
    x@result <- tmp
    p <- dotplot(x, x="Accumulator", size="Impact", showCategory=25, font.size=14)
    ggsave(filename = file.path("results/PathwayAnalysis/figures", paste0(i,"_mithril_dotplot.png")), plot = p, width = 30, height = 20, units = "cm", dpi = 600)
}

### Comment from these lines if you wish build plots using pre-computed results
if (file.exists("tmp.output.RData")) {
    unlink("tmp.output.RData")
}
saveRDS(results, file = "results/pathway.analysis.result.rds")
### The section to comment ends here!!

rm(list = ls())

results <- readRDS("results/pathway.analysis.result.rds")

for (i in names(results)) {
    r  <- results[[i]]
    fn <- file.path("results/PathwayAnalysis/tables", paste0(i,".xlsx"))
    wb <- createWorkbook()
    addWorksheet(wb, "MITHrIL")
    addWorksheet(wb, "SPIA")
    addWorksheet(wb, "Reactome")
    addWorksheet(wb, "GO")
    r$mithril$Raw.Accumulator <- NULL
    writeData(wb, "MITHrIL", r$mithril)
    r$spia$KEGGLINK <- NULL
    writeData(wb, "SPIA", r$spia)
    r$reactome$rank <- NULL
    r$reactome$leading_edge <- NULL
    r$reactome$core_enrichment <- NULL
    writeData(wb, "Reactome", r$reactome)
    r$go$results$raw_p_underrep <- NULL
    r$go$results$FWER_underrep <- NULL
    writeData(wb, "GO", r$go$results)
    saveWorkbook(wb, file = fn, overwrite = TRUE)
}




