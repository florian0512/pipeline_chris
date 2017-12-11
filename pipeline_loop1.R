

library(sarlacc)
library(ShortRead)

es_reads <- readFastq("//nfs/research2/marioni/florian/ONT_mouse/mouse_es_cell.fastq", withIds=TRUE)

read <- sread(es_reads)

UMI_1 <- NULL
UMI_2 <- NULL
chop_reads <- NULL
chop_qual <- NULL
a=1
i=100000
while(i < length(read)){
    align_data <- adaptorAlign(adaptor1 = "AGATGTGTATAAGAGACAGNNNNNNNNNNNNAGTACATCGAGA", adaptor2 = "AGATGTGTATAAGAGACAGNNNNNNNNNNNNAGTACATCGAGA", reads = read[a:i], 
                               quality = quality(getClass(quality(es_reads)))[a:i], gapOpening=1, gapExtension=5, match=2, mismatch=-5, tolerance = 100)
    a <- a+100000
    i <- i+100000
    
    chop_data <- chopReads(align_data, essential1 = TRUE, essential2 = TRUE, score1 = -25, score2 = -25)
    
    UMI1 <- umiExtract(align.stats = chop_data$adaptor1)
    UMI2 <- umiExtract(align.stats = chop_data$adaptor2)
    
    UMI_1 <- append(UMI_1, UMI1)
    UMI_2 <- append(UMI_2, UMI2)
    chop_reads <- append(chop_reads, chop_data$reads)
    chop_qual <- append(chop_qual, chop_data$quality)

}
    align_data <- adaptorAlign(adaptor1 = "AGATGTGTATAAGAGACAGNNNNNNNNNNNNAGTACATCGAGA", adaptor2 = "AGATGTGTATAAGAGACAGNNNNNNNNNNNNAGTACATCGAGA", reads = read[a:length(read)],
                           quality = quality(getClass(quality(es_reads)))[a:length(es_reads)], gapOpening=1, gapExtension=5, match=2, mismatch=-5, tolerance = 100)

    chop_data <- chopReads(align_data, essential1 = TRUE, essential2 = TRUE, score1 = -25, score2 = -25)

    UMI1 <- umiExtract(align.stats = chop_data$adaptor1)
    UMI2 <- umiExtract(align.stats = chop_data$adaptor2)

    UMI_1 <- append(UMI_1, UMI1)
    UMI_2 <- append(UMI_2, UMI2)
    chop_reads <- append(chop_reads, chop_data$reads)
    chop_qual <- append(chop_qual, chop_data$quality)

saveRDS(chop_reads, "chopreads.rds")
saveRDS(UMI_1, "UMI1.rds")
saveRDS(UMI_2, "UMI2.rds")
writeXStringSet(chop_reads, "//nfs/research2/marioni/florian/ONT_mouse/escell.fastq", format = "fastq", qualities = chop_qual)


# Mapping is now made by using minimap2 with settings to remove secondary alignments from the output and enhance mapping for nanopore.
# bsub -e min1.err -o min1.out -M 32000 "/nfs/research2/marioni/florian/minimap2/minimap2 -N=0 -ax map-ont //nfs/research2/marioni/florian/RefSeq/Mus_musculus.GRCm38.dna.sm.toplevel.fa escell.fastq > escell.sam" 

