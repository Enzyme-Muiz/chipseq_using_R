JJ<- "C:/Users/rajim/Downloads/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
buildindex(basename="GRCh38index",reference=JJ)

align(index="C:/Users/rajim/Documents/index17/GRCh37index",readfile1="D:/Muizdeen/Azzahra/D3629UT_1.fq")

library(Rsubread)
library("Rsamtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")


sortBam("D:/Muizdeen/trimmed_aligned_sorted_unique_chr.bam", "D:/Muizdeen/trimmed_aligned_sorted_unique_chr1")

indexBam("D:/Muizdeen/trimmed_aligned_sorted_unique_chr1.bam")


sortBam("D:/Muizdeen/unsortedbscore(ER) <- as.numeric(ER$name)
 
wa/bwa.bam", "D:/Muizdeen/unsortedbwa/sorted")

indexBam("D:/Muizdeen/unsortedbwa/sorted.bam")










require(ChIPpeakAnno) 
require(GenomicRanges)
require(rtracklayer)

extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")





# import our peak data
peaks1 <- import.bed("D:/azzahra_11/other_samples/prepared/3381_MIA1_chrsorted_peaks.narrowPeak", extraCols=extraCols_narrowPeak)
peaks2 <- import.bed("D:/azzahra_11/other_samples/prepared/3381_UT1_chrsorted_peaks.narrowPeak", extraCols=extraCols_narrowPeak)

intersect_peaks<- intersect(peaks1, peaks2)
length(intersect_peaks)


require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene



union.features <- assignChromosomeRegion(peaks1, TxDb=txdb, nucleotideLevel=FALSE)

pie(union.features$percentage)

bp <- barplot(union.features$percentage, ylab="%")
text(bp, union.features$percentage, signif(union.features$percentage, 4),
     pos=1)

genes <- genes(txdb)
prom <- promoters(genes, upstream=2000, downstream=2000)
#####peaks1 which are unique to to peaks1 as against peaks2
###and also with high qval
range(peaks1$qVal)
peaks1[peaks1$qVal>5]
peaks1







peaksatProm <- subsetByOverlaps(peaks1, prom)
length(peaksatProm) / length(peaks1) * 100
length(peaksatProm)
peaksatProm = findOverlaps(peaks1, prom)
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


unname(annotate::getSYMBOL(as.character(gene.ids), data='org.Hs.eg'))



peaks.top <- peaks[order(peaks$score, decreasing = TRUE)][1:100]
require(BSgenome.Hsapiens.UCSC.hg19)
peaks.seq <- getSeq(Hsapiens, peaks.top)
peaks.seq
writeXStringSet(peaks.seq, "./untreated.fasta")




peaksMIA <- import.bed("peaksMIA.bed", extraCols=extraCols_narrowPeak)


#FOXA1.top <- FOXA1[order(FOXA1$score, decreasing = TRUE)][1:100]
####score_plot

mlPvals <- score(peaks1)
hist(mlPvals, xlab="-log_10(p-value)", col="gray")

#####mia_vs_ut_3259
bp <- barplot(
    c(length(peaks1), length(peaks2)), 
    names=c("mia", "ut")
    )
# add actual values as text lables to the plot
text(bp, c(length(peaks1), length(peaks2)), labels=c(length(peaks1),
                                                length(peaks2)), pos=1)






######venn diagram

peaks1.uniq <- setdiff(peaks1, peaks2)
peaks2.uniq <- setdiff(peaks2, peaks1)

venn <- Venn(SetNames=c("MIA", "UT"),
                           Weight=c(
                                 '10'=length(peaks1.uniq),
                                 '11'=length(intersect_peaks),
                                '01'=length(peaks2.uniq)
                             )
              )
 plot(venn)
 
 
 
 ######################
 ###unique peaks1
 file1="C:/Users/rajim/Desktop/chipseq_result/whole_mia_3381_38.csv"
 peaksatProm <- subsetByOverlaps(peaks1, prom)
 length(peaksatProm) / length(peaks1) * 100
 length(peaksatProm)
 peaksatProm = findOverlaps(peaks1, prom)
 peaksatProm_genes = genes[subjectHits(peaksatProm)]
 gene.ids <- unique(names(peaksatProm_genes))
 df<- data.frame(gene.ids)
 write.csv(df, file=file1, quote=FALSE,
           row.names=FALSE)
 SP1_genes<- read.csv(file1)
 library(org.Hs.eg.db)
 library(annotate)
 
 gene_symbols1<- c()
 #for (x in gene.ids)
 for (x in SP1_genes[[1]])
 {jj<- as.character(x)
 sym<- unname(annotate::getSYMBOL(jj, data='org.Hs.eg'))
 gene_symbols1<- append(gene_symbols1, sym)
 }
 
 
 #unname(annotate::getSYMBOL(as.character(gene.ids), data='org.Hs.eg'))
 df<- data.frame(genes= gene_symbols1)
 write.csv(df, file=file1, quote=FALSE,
           row.names=FALSE)
gene_symbols1
 
 
 ####peaks2
 file2= "C:/Users/rajim/Desktop/chipseq_result/whole_UT_3381_38.csv"
 peaksatProm <- subsetByOverlaps(peaks2, prom)
 length(peaksatProm) / length(peaks2) * 100
 length(peaksatProm)
 peaksatProm = findOverlaps(peaks2, prom)
 peaksatProm_genes = genes[subjectHits(peaksatProm)]
 gene.ids <- unique(names(peaksatProm_genes))
 df<- data.frame(gene.ids)
 write.csv(df, file= file2, quote=FALSE,
           row.names=FALSE)
 SP1_genes<- read.csv(file2)
 library(org.Hs.eg.db)
 library(annotate)
 
 gene_symbols2<- c()
 #for (x in gene.ids)
 for (x in SP1_genes[[1]])
 {jj<- as.character(x)
 sym<- unname(annotate::getSYMBOL(jj, data='org.Hs.eg'))
 gene_symbols2<- append(gene_symbols2, sym)
 }
 gene_symbols2
 
 #unname(annotate::getSYMBOL(as.character(gene.ids), data='org.Hs.eg'))
 df<- data.frame(genes= gene_symbols2)
 write.csv(df, file=file2, quote=FALSE,
           row.names=FALSE)
 
 
 
 
########venn_for_mRNA_only_whole(without unique)
gene_symbols1.uniq <- setdiff(gene_symbols1, gene_symbols2)
gene_symbols2.uniq <- setdiff(gene_symbols2, gene_symbols1)
 intersect_peaks<- intersect(gene_symbols1, gene_symbols2)
 venn <- Venn(SetNames=c("MIA", "UT"),
              Weight=c(
                '10'=length( gene_symbols1.uniq),
                '11'=length(intersect_peaks),
                '01'=length( gene_symbols2.uniq)
              )
 )
 plot(venn)
 #####after_not-unique-gene-locations_mrnas_then_unique_to_ut
 file3= "C:/Users/rajim/Desktop/chipseq_result/whole_UT_3381_unique.csv"
 gene_symbols3<- setdiff(gene_symbols2.uniq, gene_symbols1.uniq)
 df<- data.frame(genes= gene_symbols3)
 write.csv(df, file=file3, quote=FALSE,
           row.names=FALSE)
 
 gene_symbols3
 
 
 
 ##############
 file3= "C:/Users/rajim/Desktop/chipseq_result/UT_3259_unique.csv"
 gene_symbols3<- setdiff(gene_symbols2, gene_symbols1)
 df<- data.frame(genes= gene_symbols3)
 write.csv(df, file=file3, quote=FALSE,
           row.names=FALSE)
 
 gene_symbols3
 
 
 ##########
 #####intersecting_genes
 file1= "C:/Users/rajim/Desktop/chipseq_result/whole_UT_3259_unique.csv"
 file2= "C:/Users/rajim/Desktop/chipseq_result/whole_UT_3424_unique.csv"
 file3= "C:/Users/rajim/Desktop/chipseq_result/whole_UT_3511_unique.csv"
 file4= "C:/Users/rajim/Desktop/chipseq_result/whole_UT_3512_unique.csv"
 file5= "C:/Users/rajim/Desktop/chipseq_result/whole_UT_3378_unique.csv"
 file6= "C:/Users/rajim/Desktop/chipseq_result/whole_UT_3381_unique.csv"
 
 
 
 
 SP1_genes_3259<- read.csv(file1)
 SP1_genes_3424<- read.csv(file2)
 SP1_genes_3511<- read.csv(file3)
 SP1_genes_3512<- read.csv(file4)
 SP1_genes_3378<- read.csv(file5)
 SP1_genes_3381<- read.csv(file6)
 length(SP1_genes_3259$genes)
 length(SP1_genes_3424$genes)
 length(SP1_genes_3511$genes)
 length(SP1_genes_3512$genes)
 length(SP1_genes_3378$genes)
 length(SP1_genes_3381$genes)
 
 
 
 
#my_genes<- intersect(intersect(SP1_genes_3259$genes, SP1_genes_3424$genes), SP1_genes_3512$genes)
 #length(my_genes)
 #my_genes<- SP1_genes_3378$genes
 #my_genes<- SP1_genes_3424$genes
 
 my_genes<- intersect(SP1_genes_3259$genes, SP1_genes_3424$genes) 
 #my_genes<- union(SP1_genes_3259$genes, SP1_genes_3424$genes) 
 #my_genes<- gene_symbols3
 file9= "C:/Users/rajim/Desktop/chipseq_result/intersect_3424_3259.csv"
 file10= "C:/Users/rajim/Desktop/chipseq_result/union_3424_3259.csv"
 my_genes_intersect<- intersect(SP1_genes_3259$genes, SP1_genes_3424$genes) 
 my_genes_union<- union(SP1_genes_3259$genes, SP1_genes_3424$genes) 
 
 df<- data.frame(genes= my_genes_intersect)
 write.csv(df, file=file9, quote=FALSE,
           row.names=FALSE)
 
 df<- data.frame(genes= my_genes_union)
 write.csv(df, file=file10, quote=FALSE,
           row.names=FALSE)
 length(my_genes_union)
 
 
 #####RITAN and RITANdata
 library(RITANdata)
 library(RITAN)
 length(my_genes)
 names(geneset_list)
 #resources <- c("ReactomePathways", "MSigDB_Hallmarks")
 resources <- c("ReactomePathways")
 #resources<- "Blood_Translaiton_Modules"
 #resources<- "KEGG_filtered_canonical_pathways"
 e <- term_enrichment(my_genes_union, resources = resources, all_symbols = cached_coding_genes)
 length(my_genes_union)
 jj<- e[c(1:10),]
 jj
 jj$name<- gsub(pattern=c("ReactomePathways."), replacement = "", x=jj$name
                )
 library(ggplot2)
 
 ggplot(jj, aes(x=-log10(q), y=name)) +geom_bar(stat= "identity", fill="red")+
        theme(panel.background = element_rect(fill="lightblue", colour = "lightblue",
        size =0.5, linetype = "solid"), plot.margin = unit(c(1, 0.0001, 1, 0.0001), "cm"))+ ylab("pathways")
 
show_active_genesets_hist()
plot(e)
 
 
 
 gggg
 gggg
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 peaks.top <- peaks[order(peaks$score, decreasing = TRUE)][1:100]
 require(BSgenome.Hsapiens.UCSC.hg19)
 
 
 
 
 
 
 
 
 
 peaks1
 
 gg
 gg
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 peaks1 <- import.bed("./macs711/NA_peaks.narrowPeak", extraCols=extraCols_narrowPeak)
 intersect_peaks<- intersect(peaks, peaks2)
 require(TxDb.Hsapiens.UCSC.hg19.knownGene)
 txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
 union.features <- assignChromosomeRegion(peaks, TxDb=txdb, nucleotideLevel=FALSE)
 
 pie(union.features$percentage)
 
 bp <- barplot(union.features$percentage, ylab="%")
 text(bp, union.features$percentage, signif(union.features$percentage, 4),
      pos=1)
 
 genes <- genes(txdb)
 prom <- promoters(genes, upstream=2000, downstream=2000)
 
 peaksatProm <- subsetByOverlaps(peaks, prom)
 length(peaksatProm) / length(peaks) * 100
 length(peaksatProm)
 peaksatProm1 = findOverlaps(peaks1, prom)
 peaksatProm_genes1 = genes[subjectHits(peaksatProm1)]
 gene.ids1 <- unique(names(peaksatProm_genes1))
 df<- data.frame(gene.ids)
 write.csv(df, file="SP1_regulated_genes.csv", quote=FALSE,
           row.names=FALSE)
 SP1_genes<- read.csv("SP1_regulated_genes.csv")
 library(org.Hs.eg.db)
 library(annotate)
 
 gene_symbols<- c()
 #for (x in gene.ids)
 for (x in SP1_genes[[1]])
 {jj<- as.character(x)
 sym<- annotate::getSYMBOL(jj, data='org.Hs.eg')
 gene_symbols<- append(gene_symbols, sym)
 }
 
 
 unname(annotate::getSYMBOL(as.character(gene.ids), data='org.Hs.eg'))
 
 
 
 peaks.top <- peaks[order(peaks$score, decreasing = TRUE)][1:100]
 require(BSgenome.Hsapiens.UCSC.hg19)
 peaks.seq <- getSeq(Hsapiens, peaks.top)
 peaks.seq
 writeXStringSet(peaks.seq, "./untreated.fasta")
 
 
 
 
 peaksMIA <- import.bed("peaksMIA.bed", extraCols=extraCols_narrowPeak)
 
 
 #FOXA1.top <- FOXA1[order(FOXA1$score, decreasing = TRUE)][1:100]
 
 
 
 
 venn <- Venn(SetNames=c("UT", "MIA"),
              Weight=c(
                '10'=length(peaks.uniq),
                '11'=length(intersect),
                '01'=length(peaksMIA.uniq)
              )
 )
 plot(venn)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 #BiocManager::install("EnrichedHeatmap")
 #BiocManager::install("genomation")
 #remotes::install_github("davetang/bedr")
 library(bedr)
 library("genomation")
 library("EnrichedHeatmap")
 require(TxDb.Hsapiens.UCSC.hg38.knownGene)
 txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
 genes <- genes(txdb)
 prom <- promoters(genes, upstream=2000, downstream=2000)
 
 peaksatProm <- subsetByOverlaps(peaks, prom)
 