library(sarlacc)
library(igraph)

UMI1 <- readRDS("UMI1.rds")
UMI2 <- readRDS("UMI2.rds")
reads <- readRDS("chopreads.rds")

groups <- readRDS("groups_final.rds")
granges <- readRDS("granges_final.rds")
flip <- granges[[2]]
granges <- granges[[1]]

#granges <- alignPrep(sam = "//nfs/research2/marioni/florian/ONT_mouse/newfull.sam", restricted = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"))
#overlaps <- findOverlaps(granges[[1]], minoverlap = 100L)
#groups <- clusterReads(overlaps)

#saveRDS(granges, "granges_final.rds")
#saveRDS(groups, "groups_final.rds")

#UMI1 <- UMI1[names(UMI1) %in% names(granges)]
#UMI2 <- UMI2[names(UMI2) %in% names(granges)]
#UMI1 <- UMI1[order(names(UMI1))]
#UMI2 <- UMI2[order(names(UMI2))]

#members <- umiGroup(UMI1, UMI2, groups, threshold=5, flip=flip)
members <- readRDS("members_final.rds")

reads <- reads[names(reads) %in% names(granges)]
reads <- reads[order(names(reads))]

#saveRDS(members, "members_final.rds")
#saveRDS(reads, "reads_final.rds")
#saveRDS(flip, "flip.rds")
#saveRDS(UMI1, "UMI1_fin.rds")
#saveRDS(UMI2, "UMI2_fin.rds")


msa <- multiReadAlign(reads = reads, groups = members, flip = flip)
conseq <- consensusReadSeq(alignments = msa)

saveRDS(msa, "msa.rds")
saveRDS(conseq, "conseq.rds")














