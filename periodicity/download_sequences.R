# load rentrez package for fetching FASTA files from NCBI database
library(rentrez)
library(stringr)

#Implement arguments
args = commandArgs(trailingOnly=TRUE)
#Search and fetch sequences from NCBI database as FASTA file
search <- entrez_search(db="nucleotide", term = args[1], use_history = TRUE)

for(seq_start in seq(1, search$count, 100)){
  recs <- entrez_fetch(db="nucleotide", web_history = search$web_history,
                       rettype="fasta", retmax= 100, retstart=seq_start)
  cat(recs, file= paste0(args[2]), append=TRUE)
  cat(seq_start+99, "sequences downloaded\r")
}

#Create readme text file with
date <- strftime(Sys.Date(),"%y-%m-%d")
read_me_txt <- c(paste0("search query: ", args[1]), paste0("date: ", date))
writeLines(read_me_txt, paste0(args[2], ".readme"))
