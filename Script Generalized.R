#Load libraries
library(Biostrings)
library(dplyr)
library(ggplot2)

#Load data necessary data
promoter_data = read.csv(file = "SL4 3k promoters.csv")

#User inputs
motif1 = "CKTV" #Can include degenerate nucleotides
motif2 = "BAMG"
input_promoters = "input_Cont.txt" #Path to input file
promoter_size = 1000 # Can be between 0 and 3000
total_motif_length = 13
gap_width = 5

# Calculate the length of the first motif
motif1_length <- nchar(motif1)

# Call the functions
hit_results1_UR <- select_truncate_and_find_motifs(promoter_data = promoter_data,
                                                   input_promoters = input_promoters,
                                                   promoter_size = promoter_size,
                                                   motif = motif1,
                                                   total_motif_length = total_motif_length)

hit_results2_UR <- search_second_motif(hit_results1_UR, gap_width, motif1_length, motif2)

hit_results1_Cont <- select_truncate_and_find_motifs(promoter_data = promoter_data,
                                                   input_promoters = input_promoters,
                                                   promoter_size = promoter_size,
                                                   motif = motif1,
                                                   total_motif_length = total_motif_length)

hit_results2_Cont <- search_second_motif(hit_results1_Cont, gap_width, motif1_length, motif2)


length(unique(hit_results2_UR$LocusID))
length(unique(hit_results2_Cont$LocusID))
cat(unique(hit_results2_UR$LocusID), sep = ", ")
length(unique(hit_results1_UR$LocusID))

# Define the function to select, truncate, and find the first motif in promoters
select_truncate_and_find_motifs <- function(promoter_data, input_promoters, promoter_size, motif, total_motif_length) {
  # Load the promoters of interest
  promoters_of_interest <- read.table(file = input_promoters, sep = "\n", stringsAsFactors = FALSE)
  
  # Filter and truncate the promoters
  truncated_promoters <- promoter_data %>%
    filter(LocusID %in% promoters_of_interest$V1) %>%
    mutate(Sequence = substr(Sequence, nchar(Sequence) - promoter_size + 1, nchar(Sequence))) %>%
    select(LocusID, Sequence)
  
  # Convert to Biostrings object
  sequences <- DNAStringSet(truncated_promoters$Sequence)
  names(sequences) <- truncated_promoters$LocusID
  
  # Find motifs using Biostrings, allowing for degenerate nucleotides
  hits <- vmatchPattern(motif, sequences, fixed = FALSE)
  
  # Extract information about the hits
  hit_info <- data.frame(LocusID = character(), Location = integer(), Pattern = character(), stringsAsFactors = FALSE)
  
  for (i in seq_along(hits)) {
    if (length(hits[[i]]) > 0) {
      for (j in seq_along(hits[[i]])) {
        locus <- names(hits)[i]
        start <- start(hits[[i]])[j]
        end <- end(hits[[i]])[j]
        pattern_start <- max(1, start)
        pattern_end <- min(width(sequences[i]), start + (total_motif_length - length(motif))) # Wider pattern
        pattern_seq <- as.character(subseq(sequences[i], pattern_start, pattern_end))
        
        hit_info <- rbind(hit_info, data.frame(LocusID = locus, Location = start, Pattern = pattern_seq))
      }
    }
  }
  return(hit_info)
}

# Define the function to search for the second motif within the hit regions
search_second_motif <- function(hit_results, gap_width, motif1_length, motif2) {
  # Create the truncatedPattern by removing the first X nucleotides from Pattern
  hit_results <- hit_results %>%
    mutate(truncatedPattern = substr(Pattern, motif1_length + gap_width + 1, nchar(Pattern)))
  
  # Convert the truncated patterns to Biostrings object
  truncated_patterns <- DNAStringSet(hit_results$truncatedPattern)
  names(truncated_patterns) <- hit_results$LocusID
  
  # Find the second motif within the truncated patterns
  hits <- vmatchPattern(motif2, truncated_patterns, fixed = FALSE)
  
  # Filter the hit_results to keep only rows with motif2 in truncatedPattern
  keep_rows <- sapply(seq_along(hits), function(i) length(hits[[i]]) > 0)
  filtered_hits <- hit_results[keep_rows, ]
  
  # Remove the truncatedPattern variable
  filtered_hits <- filtered_hits %>%
    select(-truncatedPattern)
  
  return(filtered_hits)
}
