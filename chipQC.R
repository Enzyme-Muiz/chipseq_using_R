library(ChIPQC)
samples<- read.csv("samplesheet.csv")
samples
dim(samples)
library(BiocParallel)
register(SerialParam(), default=TRUE)

register(SerialParam())
experiment = ChIPQC(samples, annotation = "hg19", chromosomes= c("chr18", "chr19", "chr2"), consensus=FALSE)
experiment = ChIPQCsample("./Input_1subread/test2_chrsorted.bam.bam")


experiment = ChIPQCsample("./sorted_rsubread.bam")
experiment1 = ChIPQCsample("./sorted_rsubread.bam", peaks="./macsresult1/NA_summits.bed")


ChIPQCreport(experiment)








samples<- read.csv("samplesheet1.csv")
samples
dim(samples)
library(BiocParallel)
register(SerialParam(), default=TRUE)

register(SerialParam())
experiment = ChIPQC(samples, annotation = "hg19")




library(ChIPQC)
library(GenomeInfoDb)
data("blacklist_hg19")
seqlevels(blacklist.hg19)
seqlevelsStyle(blacklist.hg19) <- "ncbi"
seqlevels(blacklist.hg19)



ExamplesExp <- ChIPQC(samples, annotation = "hg19", blacklist = blacklist.hg19)






remotes::install_github("imbforge/encodeChIPqc")






################################

library(data.table)
bam_file<- GenomicAlignments::readGAlignments("mia_38_3259.bam")
aln <- data.table(
  strand=as.factor(BiocGenerics::as.vector(strand(bam_file))),
  seqnames=as.factor(BiocGenerics::as.vector(seqnames(bam_file))),
  pos=ifelse(strand(bam_file) == "+", start(bam_file), end(bam_file))
)

readsPerPosition <- aln[,list(count=.N), by=list(strand, seqnames, pos)]$count

PBC <- sum(readsPerPosition == 1) / length(readsPerPosition)
PBC
