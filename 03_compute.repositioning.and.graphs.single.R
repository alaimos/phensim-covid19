library(ggplot2)
library(dplyr)

### List of important pathways to plot in the small figure
imp.pathways   <- c("path:hsa04668","path:hsa04630","path:hsa04060","path:hsa04062","path:hsa04620",
                    "path:hsa04064","path:hsa04621","path:hsa04024","path:hsa04217","path:hsa04625",
                    "path:hsa04380","path:hsa04610","path:hsa04071","path:hsa04622","path:hsa04151",
                    "path:hsa04658","path:hsa04066","path:hsa01100","path:hsa00010")

infected.base  <-  "data/simulations/virus_blanco_melo/%s.tsv"
drugs          <- c("Methylpred", "Metformin", "HCQ-CQ", "Acalabrutinib", "2DG", "Everolimus", "Dexamethasone")
drugs.names    <- c("Methylprednisolone", "Metformin", "HCQ", "Acalabrutinib", "2DG", "Everolimus", "Dexamethasone")
drugs.base     <- "data/simulations/drugs/%s/%s.tsv"
# ### Uncomment the following two lines to compute repositioning and graphs on a single cell line
# cell.line      <- "A549-ACE2_MOI_2"
# drug.cell.line <- "A549-ACE2"


read.phensim  <- function (INPUT.FILE.PHENSIM, selected.pathways) {
    phensim.output <- na.omit(read.delim(INPUT.FILE.PHENSIM, check.names = FALSE, stringsAsFactors = FALSE))
    if (!is.null(selected.pathways)) {
        phensim.output <- phensim.output[phensim.output[[1]] %in% selected.pathways, ]
    }
    pathway.output <- unique(phensim.output[,c(1,11)])
    pathway.output <- setNames(pathway.output[[2]], pathway.output[[1]])
    genes.output <- unique(phensim.output[,c(3,7)])
    genes.output <- setNames(genes.output[[2]], genes.output[[1]])
    ep.output <- unique(phensim.output[phensim.output[[5]] == "Yes",c(3,7)])
    ep.output <- setNames(ep.output[[2]], ep.output[[1]])
    ep.by.pa <- tapply(1:nrow(phensim.output), phensim.output[[1]], function(x) ((phensim.output[x,])[phensim.output[[5]][x] == "Yes",c(3,7)]))
    ep.by.pa <- lapply(ep.by.pa, function(x)(setNames(x[[2]], x[[1]])))
    gn.by.pa <- tapply(1:nrow(phensim.output), phensim.output[[1]], function(x) ((phensim.output[x,])[,c(3,7)]))
    gn.by.pa <- lapply(gn.by.pa, function(x)(setNames(x[[2]], x[[1]])))
    pathways.names <- unique(phensim.output[,c(1,2)])
    pathways.names <- setNames(pathways.names[[2]], pathways.names[[1]])
    return (list(
        pa=pathway.output,
        ge=genes.output,
        ep=ep.output,
        ep.by.pa=ep.by.pa,
        gn.by.pa=gn.by.pa,
        pn=pathways.names
    ))
}

raw.corr <- function () {
    x <- inf.data$ep
    y <- drug.data$ep
    c  <- intersect(names(x), names(y))
    x  <- x[c]
    y  <- y[c]
    xm <- mean(x)
    ym <- mean(y)
    x  <- x - xm
    y  <- y - ym
    sp <- sqrt(sum(x ^ 2)) * sqrt(sum(y ^ 2))
    a <- sum(x * y)
    return ((a / sp))
}

label.bold <- function(src, bolder) {
    b.vec <- rep("plain", length(src))
    b.vec[src %in% bolder] <- "bold"
    return (b.vec)
}

partial.corr <- function (pp) {
    x <- unlist(inf.data$ep.by.pa)
    y <- unlist(drug.data$ep.by.pa)
    c  <- intersect(names(x), names(y))
    x  <- x[c]
    y  <- y[c]
    xm <- mean(x)
    ym <- mean(y)
    x  <- x - xm
    y  <- y - ym
    sp <- sqrt(sum(x ^ 2)) * sqrt(sum(y ^ 2))
    t <- paste(pp,names(inf.data$ep.by.pa[[pp]]),sep = ".")
    a <- sum(x[intersect(t, c)] * y[intersect(t,c)])
    return ((a / sp))
}

n.tot        <- 20

df.data.all <- data.frame(name=drugs.names, Correlation=0, row.names = drugs)

for (drug in drugs) {
    cat(drug,"\n")
    dir.create(paste0("results/repositioning/details/",drug,"/"), recursive = TRUE)
    inf.data  <- read.phensim(sprintf(infected.base, cell.line), NULL)
    drug.data <- read.phensim(sprintf(drugs.base, drug, drug.cell.line), NULL)

    pcorr <- na.omit(sapply(intersect(names(inf.data$ep.by.pa), names(drug.data$ep.by.pa)), function(pp) {
        return (partial.corr(pp))
    }))

    ocorr <- na.omit(sapply(intersect(names(inf.data$ep.by.pa), names(drug.data$ep.by.pa)), function(pp) {
        a <- na.omit(inf.data$ep.by.pa[[pp]])
        b <- na.omit(drug.data$ep.by.pa[[pp]])
        c <- intersect(names(a),names(b))
        return (cor(a[c], b[c], method = "pearson"))
    }))

    df.data    <- data.frame(pathway=names(pcorr), name=inf.data$pn[names(pcorr)], Correlation=unname(pcorr))
    df.data.ep <- df.data[order(df.data$Correlation),]
    df.data.all[drug, "Correlation"] <- sum(df.data$Correlation)

    tmp.neg <- df.data.ep[df.data.ep$Correlation < 0, ]
    tmp.neg <- tmp.neg[order(tmp.neg$Correlation, decreasing = FALSE),]
    tmp.pos <- df.data.ep[df.data.ep$Correlation > 0, ]
    tmp.pos <- tmp.pos[order(tmp.pos$Correlation, decreasing = TRUE),]
    if (nrow(tmp.neg) > 0) tmp.neg <- tmp.neg[1:min(n.tot,nrow(tmp.neg)),] else tmp.neg <- NULL
    if (nrow(tmp.pos) > 0) tmp.pos <- tmp.pos[1:min(n.tot,nrow(tmp.pos)),] else tmp.pos <- NULL
    df.data.ep.filt <- rbind(tmp.neg, tmp.pos)
    dff <- setdiff(imp.pathways, df.data.ep.filt$pathway)
    if (length(dff) > 0) {
        df.data.ep.filt <- rbind(df.data.ep.filt, df.data.ep[df.data.ep$pathway %in% dff,])
    }
    df.data.ep.filt <- df.data.ep.filt[order(df.data.ep.filt$Correlation),]
    important.pathways <- df.data.ep[df.data.ep$pathway %in% imp.pathways, "name"]

    df.data.ep <- df.data.ep[,c("name", "Correlation")] %>% group_by(name) %>% summarise_all("sum")
    qqq <- quantile(abs(df.data.ep$Correlation), 0.85)
    others<- sum(df.data.ep[abs(df.data.ep$Correlation) <= qqq,]$Correlation)
    df.data.ep <- df.data.ep[abs(df.data.ep$Correlation) > qqq,]
    df.data.ep <- rbind(df.data.ep, data.frame(name="Other pathways", Correlation=others))
    
    ppl <- ggplot(df.data.ep, aes(x = reorder(name, Correlation), y = Correlation)) +
        geom_col(aes(fill=Correlation)) +
        scale_fill_gradient2(low = "blue", mid="gray", high = "orange", midpoint = 0) +
        labs(x="", y="Partial Correlation") + geom_hline(yintercept=0) +
        theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                              axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    ggsave(filename = paste0("results/repositioning/details/",drug,"/",cell.line,"_without_titles.png"), plot = ppl, width = 35, height = 20, units = "cm", dpi=600)
    
    ppl <- ggplot(df.data.ep, aes(x = reorder(name, Correlation), y = Correlation)) +
        geom_col(aes(fill=Correlation)) +
        scale_fill_gradient2(low = "blue", mid="gray", high = "orange", midpoint = 0) +
        labs(x="", y="Partial Correlation") + geom_hline(yintercept=0) +
        theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                              axis.text.x = element_text(size = 10, angle = 300, hjust = 0, vjust = 0.3),
                              axis.title.x=element_blank())
    ggsave(filename = paste0("results/repositioning/details/",drug,"/",cell.line,"_with_titles.png"), plot = ppl, width = 100, height = 20, units = "cm", dpi=600, limitsize = FALSE)
    
    ppl <- ggplot(df.data.ep.filt, aes(x = reorder(name, Correlation), y = Correlation)) +
        geom_col(aes(fill=Correlation)) +
        scale_fill_gradient2(low = "blue", mid="gray", high = "orange", midpoint = 0) +
        labs(x="", y="Partial Correlation") + geom_hline(yintercept=0) +
        theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                              axis.text.x = element_text(size = 16, angle = 290, hjust = 0, vjust = 0.3,face=label.bold(df.data.ep.filt$name, important.pathways)),
                              axis.title.x=element_blank())
    ## These lines are used to make paper figures that have the same y axis so we can line them up
    if (drug == "Methylpred" && cell.line == "A549-ACE2_MOI_0.2") {
        MY.YLIM <- range(df.data.ep.filt$Correlation)
    }
    if (drug == "Methylpred" && cell.line == "A549-ACE2_MOI_2") {
        ppl <- ppl + ylim(c(MY.YLIM))
    }
    ggsave(filename = paste0("results/repositioning/details/",drug,"/",cell.line,"_top_and_important.png"), plot = ppl, width = 30, height = 30, units = "cm", dpi=600)
}

if (!dir.exists("results/repositioning/summary/")) dir.create(paste0("results/repositioning/summary/"), recursive = TRUE)
pp1 <- ggplot(df.data.all, aes(reorder(name, Correlation), Correlation)) +
    geom_col(aes(fill=Correlation)) +
    scale_fill_gradient2(low = "green", mid="gray", high = "red", midpoint = 0) +
    labs(x="", y="Pearson Correlation") + geom_hline(yintercept=0) +
    theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                          axis.title.x=element_blank()) #,axis.text.x=element_blank(),axis.ticks.x=element_blank())
ggsave(paste0("results/repositioning/summary/",cell.line,".png"), pp1, width = 7, height = 4, dpi=600)



