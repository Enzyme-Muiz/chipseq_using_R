library(Rsubread)
library("Rsamtools")
list_one<- list.files("C:/Users/rajim/Desktop/azzahra_11/3259_readline1", full.names = TRUE)
list_two<- list.files("C:/Users/rajim/Desktop/azzahra_11/3259_readline2", full.names = TRUE)



list_one<- list.files("D:/azzahra_11/NEW_2", full.names = TRUE)
#list_two<- list.files("C:/Users/rajim/Desktop/azzahra_11/3259_readline2", full.names = TRUE)

#input_samples
one<- "C:/Users/rajim/Desktop/azzahra_11/Input_1.fq"
two<- "C:/Users/rajim/Desktop/azzahra_11/Input_2.fq"
align(index="C:/Users/rajim/Documents/index18/GRCh38index",readfile1=one, readfile2=two, output_format="BAM",type="dna", nsubreads
      = 10, nthreads=16, output_file="input_CLL.bam")



###########


align(index="C:/Users/rajim/Documents/index17/GRCh37index",readfile1="D:/Muizdeen/D3629UT_1.fq.gz", readfile2="D:/Muizdeen/D3629UT_2.fq.gz", output_format="BAM",nthreads=16, output_file="rsubread.bam")
align(index="C:/Users/rajim/Documents/index17/GRCh37index",readfile1="D:/Muizdeen/encode_chipseq.fastq.gz",output_format="BAM",minFragLength=25, nthreads=16, output_file="rsubread.bam")
align(index="C:/Users/rajim/Documents/index17/GRCh37index",readfile1="D:/Muizdeen/D3712MIA_1.fq.gz", readfile2="D:/Muizdeen/D3712MIA_2.fq.gz", output_format="BAM",nthreads=16, output_file="rsubread11.bam")
align(index="C:/Users/rajim/Documents/index17/GRCh37index",readfile1="D:/Muizdeen/D3711UT_1.fq.gz", readfile2="D:/Muizdeen/D3711UT_2.fq.gz", output_format="BAM",nthreads=16, output_file="rsubread711.bam")



align(index="C:/Users/rajim/Documents/index18/GRCh38index",readfile1=list_one[1], readfile2=list_one[2], output_format="BAM",type="dna", nsubreads
= 10, nthreads=16, output_file="UT_38_3259.bam")



sortBam("./rsubread711.bam", "./sorted_rsubread711")


indexBam("./sorted_rsubread711.bam")








sortBam("./UT_38_3259_chr.bam", "./UT_38_3259_chr_sorted")


indexBam("./UT_38_3259_chr_sorted.bam")

####sort and index input_CLL
sortBam("./input_CLL_chr.bam", "./input_CLL_chrsorted")


indexBam("./input_CLL_chrsorted.bam")


