suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
  library(gsignal)
  library(readr)
  library(biostrext)
})

# functions
# Calculate spectral density of respective period for each periodicity motif using periodicDNA and gstrings built-in functions
normalizeHistogram <- function(hist, roll_smoothed.h = 7) {
  # this function calculates the background signal with windowsize roll_smoothed.h and subtracts the bg from the meaned histogram
  h <- hist
  h <- h / sum(h, na.rm = TRUE)
  smoothed.h <- c(zoo::rollmean(h, k = roll_smoothed.h), rep(0, roll_smoothed.h-1))
  norm.h <- h - smoothed.h
  norm.h[is.na(norm.h)] <- 0
  # overwrite the first 2 distances, since they are artefacts from the background subtractions
  # you can still see the same artefacts at the end of the histogram but these are not used for periodicity calculations
  norm.h[1:2] <- 0
  # add the distance 0 and set it at the maximum to simulate the dinucleotide existing at position 0 for periodicity measurements
  new.norm.h <- rep(0, length(norm.h) + 1)
  new.norm.h[2:length(new.norm.h)] <- norm.h
  max_dist <- max(norm.h[1:50]) # only look at the beginning of the histogram (these values are more relevant IMO)
  new.norm.h[1] <- max_dist*1.3
  return(new.norm.h)
}

calculate_histogram <- function(dinucleotide, new_seqs){
  dists <- lapply(seq(1, length(new_seqs)), function(x){
    dists <- Biostrings::vmatchPattern(dinucleotide, new_seqs[x], max.mismatch = 0, fixed = FALSE)[[1]] %>%
      IRanges::start() %>%
      stats::dist() %>%
      c()
    return(dists[dists <= 300])
  }) %>% unlist()
  max_dist <- max(dists)
  hist <- hist(dists, breaks = seq(1, max_dist+1, 1), plot = FALSE)$counts
  hist <- zoo::rollmean(hist, k = 3 , na.pad = TRUE, align = 'center') # averaging window is size 3
  norm_hist <- normalizeHistogram(hist)
  norm_hist_tibble <- tibble(distance = seq(0, length(norm_hist)-1), norm_counts = norm_hist, dinucleotide = dinucleotide)
  return(norm_hist_tibble)
  gc()
}

calculate_distances_Xdist <- function(dinuc1, dinuc2, sequence){
  gc()
  pos_dinuc1 <- Biostrings::matchPattern(dinuc1, sequence, max.mismatch = 0, fixed = FALSE) %>%
    IRanges::start() %>%
    c()
  
  pos_dinuc2 <- Biostrings::matchPattern(dinuc2, sequence, max.mismatch = 0, fixed = FALSE) %>%
    IRanges::start() %>%
    c()
  
  dists <- outer(pos_dinuc1, pos_dinuc2, "-") %>% abs() %>% c() %>% tibble() %>% dplyr::filter(. < 300) %>% pull(.)
  
  return(dists)
}

calculate_histogram_from_dists <- function(dists){
  hist <- hist(dists, breaks = max(dists), plot = FALSE)$counts
  hist <- zoo::rollmean(hist, k = 3 , na.pad = TRUE, align = 'center') # averaging window is size 3
  norm_hist <- normalizeHistogram(hist)
  norm_hist_tibble <- tibble(distance = seq(0, length(norm_hist)-1), norm_counts = norm_hist)
  return(norm_hist_tibble)
}

# fit function
fit_sine <- function(df, l_start = 20, A_start = 1e-3, H_start = 20){
  nlsfit <- nls(
    norm_counts ~ exp(-1*(log(2)*distance)/abs(H)) * abs(A) * sin((((2*pi*distance)/abs(l))) + (pi/2)),
    data.frame(df),
    start=list(l=l_start, A=A_start, H=H_start),
    control = nls.control(maxiter = 10000, minFactor = 1/2048))
  gof <- nlsfit$m$deviance() / sd(df %>% pull(norm_counts))
  
  return(list(GOF=gof, fit=nlsfit))
}


#implement arguments

# 1: in fasta file
# 2: out file basename
# 3: family
# 4: genus
args = commandArgs(trailingOnly=TRUE)

hist_csv <- paste0(args[2], args[4], ".Xdist.HIST", ".csv")

fam <- args[3]
gen <- args[4]

seq <- readDNAStringSet(args[1])

# prepare sequences
split_len <- 10000

#If long sequences are present (>10000 bp) in DNAStringSet split long sequences into shorter sequence chunks (=< 10000) and append them to a new list of sequences, else append all (short) sequences (<10000 bp) of DNAStringSet to new list of sequences
short_seqs <- seq[which(width(seq) <= split_len)]
long_seqs <- seq[which(width(seq) > split_len)]
my_list <- list()
if(length(long_seqs) > 0){
  long_seqs <- do.call(c, long_seqs)
  for (i in 1:length(long_seqs)){
    split_long_seqs <- split_in_windows(long_seqs[[i]], split_len, overlap = 0)
    my_list[[i]] <- split_long_seqs
  }
  split_long_seqs <- do.call(c, my_list)
  
  # filter out the split sequences that are now shorter
  split_long_seqs <- split_long_seqs[which(width(split_long_seqs) >= split_len)]
  
  short_and_split_long_seqs <- DNAStringSetList(short_seqs, split_long_seqs)
  short_and_split_long_seqs <- unlist(short_and_split_long_seqs)
  rm(short_seqs)
  rm(long_seqs)
  invisible(gc())
} else {
  short_and_split_long_seqs <- do.call(c, short_seqs)
  rm(short_seqs)
  rm(long_seqs)
  invisible(gc())
}

#From new list of sequences, create new DNAStringSet of same width as Dependoparvoviridae DNAStringSet (i.e dataset of smallest width) for comparability
subsampled_seqs <- DNAStringSet()

cutoff <-  100000
counter <- 1
chosen_ind <- list()
# randomly choose sequences out of the pool until the length reaches the cutoff
while(sum(width(subsampled_seqs)) < cutoff){
  ind <- sample.int(length(short_and_split_long_seqs), 1, replace = F)
  # if this index has not been chosen before add a new sequence to the list
  # if this index has already been selected once, skip it
  if (!ind %in% chosen_ind){
    subsampled_seqs[[counter]] <- short_and_split_long_seqs[[ind]]
    chosen_ind <- append(chosen_ind, ind)
    counter <- counter + 1
  } else if (length(chosen_ind) == length(short_and_split_long_seqs)){
    break
  } else {
    next
  }
}

# Calculation
# also add the reverse complementary of every chosen sequence
new_seqs <- DNAStringSet()
new_seqs <- list()
for (x in seq(length(subsampled_seqs))){
  new_seqs <- append(new_seqs, subsampled_seqs[x])
  #new_seqs <- append(new_seqs, Biostrings::reverseComplement(subsampled_seqs[x]))
}

# what dinucleotides
dinucleotides_anchor <- list("YY", "RR", "SS", "WW")
dinucleotides_2 <- list("YY", "RR", "SS", "WW")


# calculate the histogram
hist <- tibble()
for (dinuc1 in dinucleotides_anchor){
  for (dinuc2 in dinucleotides_2){
    
    f1 <- lapply(seq, function(x){return(calculate_distances_Xdist(dinuc1, dinuc2, x))}) %>%
      unlist(use.names = F) %>% 
      calculate_histogram_from_dists() %>% 
      mutate(dinuc1 = dinuc1, dinuc2 = dinuc2)
    
    hist <- bind_rows(hist, f1)
  }
}

# output
write.csv(hist, file = hist_csv)