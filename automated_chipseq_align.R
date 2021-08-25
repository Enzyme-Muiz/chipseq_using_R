for (dirs in list.dirs('.', recursive = FALSE)[c(8, 9, 10)])
{fastqs<- list.files(dirs, full.names = TRUE)
  Rsubread::align(index="C:/Users/rajim/Documents/index18/GRCh38index",readfile1=fastqs[1], readfile2=fastqs[2], output_format="BAM",type="dna", nsubreads
        = 10, nthreads=16, output_file=paste0(stringr::str_remove(basename(fastqs[1]), ".fastq"),".bam"))
  
}

for (files2 in list.files('.', full.names = TRUE))
{
  if (grepl('_chr.bam$', files2))
  {
    Rsamtools::sortBam(files2, paste0("./",tools::file_path_sans_ext(files2), "sorted"))
    Rsamtools::indexBam(paste0("./",tools::file_path_sans_ext(files2), "sorted.bam"))
    }
}
