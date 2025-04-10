suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
  library(gsignal)
  library(readr)
  library(biostrext)
})

# functions
# Calculate spectral density of respective period for each periodicity motif using periodicDNA and gstrings built-in functions
normalize_histogram <- function(hist, roll_smoothed.h = 7) {
  # this function calculates the background signal with windowsize roll_smoothed.h and subtracts the bg from the meaned histogram
  h <- hist / sum(hist, na.rm = TRUE)
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
  norm_hist <- normalize_histogram(hist)
  norm_hist_tibble <- tibble(distance = seq(0, length(norm_hist)-1), norm_counts = norm_hist, dinucleotide = dinucleotide)
  return(norm_hist_tibble)
  gc()
}

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

hist_csv <- paste0(args[2], ".HIST", ".csv")
curve_csv <- paste0(args[2], ".CURVE", ".csv")

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
  new_seqs <- append(new_seqs, Biostrings::reverseComplement(subsampled_seqs[x]))
}

# what dinucleotides
dinucleotides <- list(
  "AA", "TT", "TA", "AT", "WW",
  "GG", "CC", "CG", "GC", "SS",
  "TC", "CT", "YY",
  "AG", "GA", "RR",
  "YR", "RY", "WS", "SW", "YS", "SY", "YW", "WY", "WR", "RW")


# calculate the histogram
hist_master <- lapply(dinucleotides, calculate_histogram, new_seqs=new_seqs) %>%
  bind_rows() %>% 
  dplyr::filter(distance < 295)

# maximum tries for the curve fitting
max_tries <- 300

# fit the curve
df <- tibble()
failed <- list()
for (dnoi in dinucleotides){
  max_dist <- 80
  # initialisation parameters
  A_start <- 1e-3
  l_start <- 20
  H_start <- 20
  # try to fit the curve to the data of the dinucleotide 100 times
  # every time re-initialize the parameters
  tries <- 1
  fit_list <- list()
  GOF_list <- rep(NA, max_tries)
  while(tries < (max_tries+1)){
    # init the default fit that is put out if no fit can be found with the chosen parameters
    f <- list("fit" = NA, "GOF" = Inf)
    # try to fit the sine curve
    try(f <- fit_sine(
      hist_master %>% dplyr::filter(dinucleotide == dnoi, distance < max_dist),
      A_start = A_start,
      l_start = l_start,
      H_start = H_start)
      )
    # try to extract the Goodness of Fit
    try(GOF <- f$GOF)
    
    # append the fit and GOF to the list of fits and GOFs
    # add a failsave for very small periods
    if (!any(is.na(f$fit))){
      if (f$fit$m$getPars()["l"] > 3){
        fit_list[[tries]] <- f$fit
        GOF_list[tries] <- f$GOF
      }
      else{
        fit_list[[tries]] <- NA
        GOF_list[tries] <- Inf
      }}
      else{
        fit_list[[tries]] <- NA
        GOF_list[tries] <- Inf
      }
    
    # init random variables
    A_start = sample(seq(5e-4, 1e-2, 1e-4), 1)
    l_start = runif(n=1, min = 5, max = 30)
    H_start = runif(n=1, min = 20, max = 50)
    # counter
    tries <- tries + 1
  }
  
  # out of the found fits, get the fit index with the lowest GOF
  best_fit <- which(unlist(GOF_list) == min(unlist(GOF_list)))
  # sometimes there are multiple with the same GOF. In that case just choose the first one
  if (length(best_fit) > 1){
    best_fit <- best_fit[1]
  }
  
  GOF <- GOF_list[best_fit]
  
  # failsave if none of the 100 tries succeeded in fitting a curve
  if (GOF == Inf){
    failed <- append(failed, dnoi)
    list_df <- c(
      "Family"=fam,
      "Genus"=gen,
      "dinucleotide"=dnoi,
      "nGOF"=NA,
      "period"=NA,
      "amplitude"=NA,
      "HL"=NA
    )
    df <- df %>% bind_rows(list_df)
    next
  }
  
  # get the best model
  f <- fit_list[[best_fit]]
  # predict the data with the best fitting curve
  hist_master_pred <- hist_master %>% 
    dplyr::filter(dinucleotide == dnoi, distance < max_dist) %>% 
    mutate(
      pred_alg = stats::predict(f, newdata = hist_master %>% dplyr::filter(dinucleotide == dnoi, distance < max_dist)))
  # get the normalized Goodness Of Fit, Period, half-life, and amplitude
  nGOF <- GOF / sd(hist_master_pred %>% pull(norm_counts))
  period <- unname(f$m$getPars()["l"])
  amplitude <- unname(f$m$getPars()["A"])
  HL <- unname(f$m$getPars()["H"])
  # plot and save
  p <- hist_master_pred %>% 
    ggplot(aes(x=distance)) +
    geom_line(aes(y=norm_counts), color="black") +
    geom_line(aes(y=pred_alg), color="red") +
    xlim(0, max_dist) +
    ggtitle(paste(dnoi, round(nGOF, 2)), subtitle=paste(round(period, 2), format(amplitude, scientific = T)))
  ggsave(paste0(args[2], ".", dnoi, ".HISTCURVE", ".png"), plot=p)
  # output df
  list_df <- c(
    "Family"=fam,
    "Genus"=gen,
    "dinucleotide"=dnoi,
    "GOF"=GOF,
    "nGOF"=nGOF,
    "period"=abs(period),
    "amplitude"=abs(amplitude),
    "HL"=abs(HL)
  )
  df <- df %>% bind_rows(list_df)
}

# format the table
df <- df %>% 
  mutate(
    nGOF = as.double(nGOF),
    period = as.double(period),
    amplitude = as.double(amplitude),
    HL = as.double(HL)
  )

# output
write.csv(hist_master, file = hist_csv)
write.csv(df, file = curve_csv)