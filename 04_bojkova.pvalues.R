library(xlsx)
library(VennDiagram)

combine <- function (p1 = NULL, p2 = NULL) {
    tm = na.omit(c(p1, p2))
    if (!all(tm >= 0 & tm <= 1)) {
        stop("values of p1 and p2 have to be >=0 and <=1 or NAs")
    }
    ## Deals with machine error
    p1[p1 <= .Machine$double.eps] <- .Machine$double.eps
    p2[p2 <= .Machine$double.eps] <- .Machine$double.eps
    p1[p1 >= 0.9999999] <- 0.9999999
    p2[p2 >= 0.9999999] <- 0.9999999
    comb = pnorm((qnorm(p1) + qnorm(p2))/sqrt(2))
    comb[is.na(p1)] <- p2[is.na(p1)]
    comb[is.na(p2)] <- p1[is.na(p2)]
    return(comb)
}


infected.cell.lines    <- c("bojkova_proteome_2hrs", "bojkova_proteome_6hrs", "bojkova_proteome_10hrs", "bojkova_proteome_24hrs")
phensim.base           <- "data/simulations/virus_bojkova/%s.tsv"
proteome.base          <- "data/LFCS/virus_bojkova/%s.tsv"
bojkova_sel_pathways   <- readRDS("data/selected_pathways.rds")

APPEND <- FALSE
for (i in infected.cell.lines) {
    print(i)
    list.of.all.vars <- ls()
    INPUT.FILE.PHENSIM      <- sprintf(phensim.base, i)
    INPUT.FILE.PROTEOME.TSV <- sprintf(proteome.base, i)
    phensim.output          <- read.delim(INPUT.FILE.PHENSIM)
    proteome                <- read.table(INPUT.FILE.PROTEOME.TSV, sep="\t", header = FALSE)
    phensim.output$Input    <- phensim.output$Node.Id %in% proteome[[1]]
    
    bojkova_pathways <- unique(phensim.output[,1:2])
    tbl.results <- as.data.frame(matrix(NA, nrow=nrow(bojkova_pathways), ncol=16,dimnames = list(bojkova_pathways[[1]],c(
        "PathwayId",
        "PathwayName",
        "Selected",
        "UPINPUT",
        "DOWNINPUT",
        "UPPRED",
        "DOWNPRED",
        "NONPRED",
        "UPPERT",
        "DOWNPERT",
        "COMPLETEPERT",
        "PV",
        "ADJ",
        "PVG",
        "PVComb",
        "FDRComb"
    ))))
    tbl.results$PathwayId   <- bojkova_pathways[[1]]
    tbl.results$PathwayName <- bojkova_pathways[[2]]
    tbl.results$Selected    <- bojkova_pathways[[1]] %in% bojkova_sel_pathways[[1]]
    
    all.nodes <- length(unique(phensim.output$Node.Id))
    all.pert.nodes <- length(unique(phensim.output$Node.Id[phensim.output$Average.Node.Perturbation != 0]))
    sel.genes <- list("path:hsa00071"=NULL,"path:hsa00520"=NULL,"path:hsa00010"=NULL,"path:hsa00020"=NULL,"path:hsa00230"=NULL,"path:hsa01200"=NULL,"path:hsa00240"=NULL)
    for (p in bojkova_pathways[[1]]) {
        tmp <- phensim.output[phensim.output[[1]] == p,]
        if (nrow(tmp) > 0) {
            all.nodes.p      <- length(unique(tmp$Node.Id))
            all.pert.nodes.p <- length(unique(tmp$Node.Id[tmp$Average.Node.Perturbation != 0]))
            tbl.results[p,4] <- length(which(tmp$Node.Id != "CV19" & tmp$Input & tmp$Average.Node.Perturbation > 0))
            tbl.results[p,5] <- length(which(tmp$Node.Id != "CV19" & tmp$Input & tmp$Average.Node.Perturbation < 0))
            tbl.results[p,6] <- length(which(tmp$Node.Id != "CV19" & !tmp$Input & tmp$Average.Node.Perturbation > 0))
            tbl.results[p,7] <- length(which(tmp$Node.Id != "CV19" & !tmp$Input & tmp$Average.Node.Perturbation < 0))
            tbl.results[p,8] <- length(which(tmp$Node.Id != "CV19" & !tmp$Input & tmp$Average.Node.Perturbation == 0))
            tbl.results[p,9] <- sum(tmp$Average.Node.Perturbation[tmp$Node.Id != "CV19" & tmp$Average.Node.Perturbation > 0])
            tbl.results[p,10] <- sum(tmp$Average.Node.Perturbation[tmp$Node.Id != "CV19" & tmp$Average.Node.Perturbation < 0])
            tbl.results[p,11] <- unique(tmp$Average.Pathway.Perturbation)
            tbl.results[p,12] <- unique(tmp$Pathway.p.value)
            tbl.results[p,13] <- unique(tmp$Pathway.Adjusted.p.value)
            tbl.results[p,14] <- phyper(q=all.pert.nodes.p - 1, m=all.nodes.p, n = all.nodes - all.nodes.p, k = all.pert.nodes, lower.tail = FALSE) 
            if (p %in% names(sel.genes)) {
                sel.genes[[p]] <- tmp$Node.Id[tmp$Node.Id != "CV19" & tmp$Average.Node.Perturbation != 0]
            }
        }
    }
    tbl.results[,12][tbl.results[,12] > 1 ] <- 1
    tbl.results[,15] <- combine(tbl.results[,12], tbl.results[,14])
    tbl.results[,16] <- p.adjust(tbl.results[,15], method = "fdr")
    venn.diagram(sel.genes[1:5], paste0("results/venn.",i,"_1.tiff"), fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), alpha = 0.50, cat.cex = 0.4, cat.pos = c(0,10,120,120,40))
    venn.diagram(sel.genes[3:7], paste0("results/venn.",i,"_2.tiff"), fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), alpha = 0.50, cat.cex = 0.4, cat.pos = c(0,10,120,120,40))
    write.xlsx(tbl.results, file = "results/output_bojkova_pvalues.xlsx", sheetName = i, append = APPEND)
    APPEND <- TRUE
    rm(list = setdiff(ls(), list.of.all.vars))    
}