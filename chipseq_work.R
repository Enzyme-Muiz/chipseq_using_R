####genomation
library(genomation)

require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
prom <- promoters(genes, upstream=2000, downstream=2000)
head(prom)
length(prom)
seqnames(prom)
ranges(prom)






##window data
q= prom[seqnames(prom)=="chr12"]
length(q)
#####target data
#########
bed_to_granges <- function(file, header = FALSE){
  df <- read.table(file,
                   header = header,
                   stringsAsFactors = FALSE)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  library("GenomicRanges")
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}
###########
cage_like= bed_to_granges("3378_MIA1_chrsorted_summits.bed")


####scorematrix
sm = ScoreMatrix(target = cage_like, windows = q)
sm.scaled = scaleScoreMatrix(sm)

heatMatrix(sm.scaled, xcoords = c(-1000, 1000), col = c("lightgray", "blue"))

plotMeta(sm, xcoords = c(-2000, 2000))


peaksatProm = genes[subjectHits(prom)]
gene.ids1 <- unique(names(prom))











library(genomation)
data("cage")
data(promoters)
sm = ScoreMatrix(target = cage, windows = promoters)
sm.scaled = scaleScoreMatrix(sm)
heatMatrix(sm.scaled, xcoords = c(-1000, 1000), col = c("lightgray", "blue"))
plotMeta(sm, xcoords = c(-1000, 1000))




##########sp1 only genes






























library(DiffBind)
setwd("D:/azzahra_11/other_samples/prepared")
samples<- read.csv("diffbind_on_SP1_cll_2.csv")
samples
#samples<- samples[ c(1,4,6,7,10,12),]
#SP1<- dba(sampleSheet="diffbind_on_SP1_cll.csv")
SP1<- dba(sampleSheet= samples)

#SP1<- dba.blacklist(SP1, blacklist="TRUE")

SP1
dev.off()
plot(SP1)


SP1_OLAP_RATE<- dba.overlap(SP1, mode= DBA_OLAP_RATE)
SP1_OLAP_RATE
plot(SP1_OLAP_RATE, type= "b")

consensus<- dba.peakset(SP1,minOverlap=5, bRetrieve = TRUE)
width(consensus)
plot(SP1)
SP1
SP1<- dba.count(SP1)
dev.off()
SP1
plot(SP1)
dev.off()

dba.plotPCA(SP1)
dba.plotPCA(SP1, DBA_CONDITION)

#SP1 <- dba.normalize(SP1)
#SP1.model <- dba.contrast(SP1, categories=DBA_CONDITION, minMembers = 2)
SP1.model <- dba.contrast(SP1, categories=c(DBA_CONDITION, DBA_TISSUE))


#SP1


#SP1.model<- dba.contrast(SP1)
#SP1.model<- dba.contrast(SP1.model)
SP1.analysis <- dba.analyze(SP1.model, method= DBA_ALL_METHODS)
SP1.analysis
dev.off()
plot(SP1.analysis, contrast=1)
dev.off()
SP1.db<- dba.report(SP1.analysis)
SP1.db
dba.plotVolcano(SP1.analysis)

subset(SP1.db, Fold > 0, select= Fold)
#subset(gr, strand == "+" & score > 5, select = score)

require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
prom <- promoters(genes, upstream=2000, downstream=2000)


peaksatProm = findOverlaps(SP1.db, prom)

peaksatProm_genes = genes[subjectHits(peaksatProm)]
gene.ids <- unique(names(peaksatProm_genes))
df<- data.frame(gene.ids)
write.csv(df, file="SP1_regulated_genes.csv", quote=FALSE,
          row.names=FALSE)
SP1_genes<- read.csv("sp1_regulated_genes.csv")
library(org.Hs.eg.db)
library(annotate)

gene_symbols<- c()
#for (x in gene.ids)
for (x in SP1_genes[[1]])
{jj<- as.character(x)
sym<- annotate::getSYMBOL(jj, data='org.Hs.eg')
gene_symbols<- append(gene_symbols, sym)
}


genes_downregulated_in_mia<- unname(annotate::getSYMBOL(as.character(gene.ids), data='org.Hs.eg'))
print(genes_downregulated_in_mia)


library(RITANdata)
library(RITAN)
names(geneset_list)
resources<- c("ReactomePathways")
e<- term_enrichment(genes_downregulated_in_mia, resources= resources, all_symbols =  cached_coding_genes)
e
