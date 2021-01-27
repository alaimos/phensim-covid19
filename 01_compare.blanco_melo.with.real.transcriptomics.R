library("org.Hs.eg.db")
library(dplyr)
library(HGNChelper)
library(openxlsx)

infected.cell.lines    <- c("A549-ACE2_MOI_0.2", "A549-ACE2_MOI_2", "A549_MOI_0.2", "A549_MOI_2", "Calu-3", "NHBE")
phensim.base           <- "data/simulations/virus_blanco_melo/%s.tsv"
proteome.base          <- "data/LFCS/virus_blanco_melo/%s.csv"
non.exp.base           <- "data/simulations/virus_blanco_melo/%s_non_exp.txt"
FINAL.OUTPUT.FILE      <- "results/blanco_melo_comparison_results.xlsx"
COLUMN.PROTEOME        <- 1

mtx.results <- matrix(NA, nrow = length(infected.cell.lines), ncol = 9, 
                      dimnames = list(infected.cell.lines, 
                                      c("Overall.Accuracy", "Predicted.Accuracy", "Predicted.PPV", "Predicted.SENS", "Predicted.SPEC", "Non.Predicted.PPV", "Non.Predicted.SENS", "Non.Predicted.SPEC", "Non.Predicted.Accuracy")))
for (i in infected.cell.lines) {
    list.of.all.vars <- ls()
    INPUT.FILE.PHENSIM      <- sprintf(phensim.base, i)
    INPUT.FILE.PROTEOME.CSV <- sprintf(proteome.base, i)
    NON.EXP                 <- sprintf(non.exp.base, i)

    phensim.output <- read.delim(INPUT.FILE.PHENSIM, check.names = FALSE)
    phensim.output <- phensim.output[!is.na(phensim.output[[3]]),]
    suppressMessages(suppressWarnings({
        proteome <- read.csv(INPUT.FILE.PROTEOME.CSV, row.names=1, sep=";")
        corr.symb <- checkGeneSymbols(rownames(proteome), unmapped.as.na = FALSE, species = "human")
        map <- setNames(strsplit(corr.symb$Suggested.Symbol, split = "\\s+///\\s+", perl = TRUE), corr.symb$x)
        df.map <- na.omit(data.frame(name=rep(names(map), sapply(map, length)), name_new=unname(unlist(map))))
        map <- mapIds(org.Hs.eg.db, df.map$name_new, 'ENTREZID', 'SYMBOL', multiVals = "list")
        df.map <- df.map %>% inner_join(na.omit(data.frame(name=rep(names(map), sapply(map, length)), id=unname(unlist(map)))), by = c("name_new"="name"))
        df.map <- df.map[,c("name", "id")]
        prot.map <- data.frame(name=rownames(proteome), proteome)
        prot.map <- suppressWarnings(prot.map %>% inner_join(df.map, by="name") %>% group_by(id) %>% summarise_all(mean))
        tmp <- prot.map[,colnames(proteome)]
        rownames(tmp) <- prot.map$id
        class(tmp) <- "data.frame"
        proteome <- tmp
        rm(map, df.map, prot.map, tmp)
    }))
    
    all.genes <- unique(phensim.output[,3])
    all.ep    <- unique(phensim.output[phensim.output[[5]] == "Yes",3])
    all.zero.genes <- unique(phensim.output[phensim.output[[7]] == 0,3])
    all.zero.ep    <- unique(phensim.output[phensim.output[[7]] == 0 & phensim.output[[5]] == "Yes",3])
    tmp.df <- unique(phensim.output[,c(3,7)])
    perts.complete    <- factor(sign(setNames(tmp.df[[2]], tmp.df[[1]])), levels = c(1,-1,0))
    tmp.df <- unique(phensim.output[phensim.output[[5]] == "Yes", c(3,7)])
    perts.complete.ep <- factor(sign(setNames(tmp.df[[2]], tmp.df[[1]])), levels = c(1,-1,0))
    orig.data.all     <- setNames(proteome[,COLUMN.PROTEOME], rownames(proteome))
    common.all        <- intersect(names(perts.complete), names(orig.data.all))
    orig.data.all     <- orig.data.all[common.all]
    perts.complete    <- perts.complete[common.all]

    tmp.orig.data     <- sign(orig.data.all)
    tmp.orig.data[abs(orig.data.all) <= 0.6] <- 0
    orig.data.all     <- factor(tmp.orig.data, levels = c(1,-1,0))
    rm(tmp.orig.data, tmp.df)
    common.ep  <- intersect(names(perts.complete.ep), names(orig.data.all))
    table.all  <- table(real=orig.data.all[common.all], predicted=perts.complete[common.all])
    table.ep   <- table(real=orig.data.all[common.ep], predicted=perts.complete.ep[common.ep])
    print(i)
    print(table.all)
    print(1 - (table.all[3,3] / length(which(perts.complete[common.all] == 0))))
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
    mtx.results[i,1]   <- accuracy(table.all)
    mtx.results[i,2]   <- accuracy(table.all[1:2,1:2])
    mtx.results[i,3:5] <- sens.spec(table.all, 1, 2)
    mtx.results[i,6:8] <- sens.spec(table.all, 3)
    mtx.results[i,9]   <- table.all[3,3] / length(which(orig.data.all[common.all] == 0))
    rm(list = setdiff(ls(), list.of.all.vars))    

}

write.xlsx(mtx.results, file = FINAL.OUTPUT.FILE, col.names = TRUE, row.names = TRUE)