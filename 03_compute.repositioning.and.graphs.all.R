

cell.lines      <- c("A549-ACE2_MOI_0.2", "A549-ACE2_MOI_2", "A549_MOI_0.2", "A549_MOI_2", "Calu-3", "NHBE")
drug.cell.lines <- c("A549-ACE2",         "A549-ACE2",       "A549",         "A549",       "Calu-3", "NHBE")
# These line is used to make aligned paper figures (same y axis range)
# MY.YLIM         <- NULL

for (iiiii in 1:length(cell.lines)) {
    my.vars <- ls()
    cell.line      <- cell.lines[iiiii]
    drug.cell.line <- drug.cell.lines[iiiii]
    print(cell.line)
    source("./03_compute.repositioning.and.graphs.single.R")
    rm(list = setdiff(ls(), my.vars))
}


