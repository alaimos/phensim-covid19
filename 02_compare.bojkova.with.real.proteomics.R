library("org.Hs.eg.db")
library(dplyr)
library(HGNChelper)
library(openxlsx)
library(ggplot2)
library(ggpubr)

infected.cell.lines    <- c("bojkova_proteome_2hrs", "bojkova_proteome_6hrs", "bojkova_proteome_10hrs", "bojkova_proteome_24hrs")
phensim.base           <- "data/simulations/virus_bojkova/%s.tsv"
proteome.base          <- "data/LFCS/virus_bojkova/%s.csv"
non.exp.base           <- "data/simulations/virus_bojkova/%s_non_exp.txt"
FINAL.OUTPUT.FILE      <- "results/bojkova_comparison_results.xlsx"
COLUMN.PROTEOME        <- 1

mtx.results <- as.data.frame(matrix(NA, nrow = length(infected.cell.lines), ncol = 7, 
                      dimnames = list(infected.cell.lines, 
                                      c("Accuracy", "Predicted.PPV", "Predicted.SENS", "Predicted.SPEC", "Predicted.Percentage", "Common.Proteins", "All.Proteins"))))
for (i in infected.cell.lines) {
    list.of.all.vars <- ls()
    INPUT.FILE.PHENSIM      <- sprintf(phensim.base, i)
    INPUT.FILE.PROTEOME.CSV <- sprintf(proteome.base, i)
    NON.EXP                 <- sprintf(non.exp.base, i)

    phensim.output <- na.omit(read.delim(INPUT.FILE.PHENSIM, check.names = FALSE))
    proteome <- read.csv(INPUT.FILE.PROTEOME.CSV, row.names=1, sep=";", header = FALSE)
    all.proteins <- length(unique(rownames(proteome)))
    all.genes <- unique(phensim.output[,3])
    all.zero.genes <- unique(phensim.output[phensim.output[[7]] == 0,3])
    tmp.df <- unique(phensim.output[,c(3,7)])
    perts.complete    <- factor(sign(setNames(tmp.df[[2]], tmp.df[[1]])), levels = c(1,-1,0))
    orig.data.all     <- setNames(proteome[,COLUMN.PROTEOME], rownames(proteome))
    common.all        <- intersect(names(perts.complete), names(orig.data.all))
    orig.data.all     <- orig.data.all[common.all]
    perts.complete    <- perts.complete[common.all]
    tmp.orig.data     <- sign(orig.data.all)
    orig.data.all     <- factor(tmp.orig.data, levels = c(1,-1,0))
    rm(tmp.orig.data, tmp.df)
    table.all  <- table(real=orig.data.all[common.all], predicted=perts.complete[common.all])

    accuracy  <- function (table) {
        return (sum(diag(table))/sum(table))
    }
    
    sens.spec <- function (table, pos, neg=((1:nrow(table))[-pos])) {
        n    <- nrow(table)
        rest <- neg #(1:n)[-pos]
        TP   <- table[pos,pos]
        FP   <- sum(table[rest, pos])
        FN   <- sum(table[pos, rest])
        TN   <- sum(diag(table)[-pos])
        return (c(
            PPV=(TP/(TP+FP)),
            SENS=(TP/(TP+FN)),
            SPEC=(TN/(TN+FP))
        ))
    }
    mtx.results[i,1]   <- accuracy(table.all[1:2,1:2])
    mtx.results[i,2:4] <- sens.spec(table.all, 1, 2)
    mtx.results[i,5]   <- length(which(perts.complete!=0)) / length(perts.complete)
    mtx.results[i,6]   <- length(perts.complete)
    mtx.results[i,7]   <- all.proteins
    print(table.all)
    rm(list = setdiff(ls(), list.of.all.vars))    

}
write.xlsx(data.frame(mtx.results), file = FINAL.OUTPUT.FILE, col.names = TRUE, row.names = TRUE)