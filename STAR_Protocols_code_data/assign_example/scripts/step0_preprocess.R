# Begin script - this script is for pre-processing the mapped barcode reads
# Summarize reads from all lanes for each replicate

# Define a utility function to remove duplicated strains and keep the max sum
remove_duplicates <- function(df, name_col, sum_col) {
  # Ensure that the dataframe has the required columns
  if (!all(c(name_col, sum_col) %in% names(df))) {
    stop(paste("The dataframe must contain columns", name_col, "and", sum_col))
  }
  # Identify the duplicated strain names
  duplicate_values <- unique(df[[name_col]][duplicated(df[[name_col]])])
  print("Duplicated strains (Keeping the one with highest read sum):")
  print(duplicate_values)
  result <- df[!df[[name_col]] %in% duplicate_values, ]
  # Loop through each unique ID and keep the row with the highest SUM
  for (id in duplicate_values) {
    # Subset the data frame for the current ID
    subset_df <- df[df[[name_col]] == id, ]
    # Find the row with the maximum SUM value
    max_row <- subset_df[which.max(subset_df[[sum_col]]), ]
    # Append the max_row to the result data frame
    result <- rbind(result, max_row)
  }
  return(result)
}

# Set working directory
work_dir <- "C.albicans_in_vitro_pipeline/STAR_Protocols_code_data/assign_example/"
setwd(work_dir)

########################################
# PARAMETER SETTING & DATA PREPROCESSING
########################################

# Set important parameters - please change based on needs
low_sum_reads_filter <- 150 # Remove strains with very low or inconsistent read counts across samples
bool_remove_duplicates <- TRUE # decide whether to keep unique occurrence of strain
map_dir <- "./input/map_replicate_reads.txt"
read_dir <- "./output/MappingOutput-2024-11-17/merged_results.txt"
save_sum_raw_dir <- "./output/sum_data_numeric.txt"
save_sum_final_dir <- "./output/sum_data.txt"

# Read the meta map that guides the replicate column names to refer to
df_map <- read.csv(file=map_dir, header=TRUE, sep='\t')

# Read the merged_results.txt file generated by sum_alignments.sh
# Setting check.names = FALSE will ensure R does not alter the column names
# Otherwise, column names that start with a number or contain spaces or other 
# special characters will be prefixed with “X” and thus cannot be mapped well.
df_reads <- read.csv(file=read_dir, header=TRUE, sep='\t', check.names = FALSE)

# Initialize another dataframe to collect the summarized reads
# The basic info columns can be changed based on the new input data
df_sum <- df_reads[, c("orf19", "Feature", "Common", "Description", 
                       "GRACE.Plate", "GRACE.Position")]
colnames(df_sum) <- c("orf19", "feature", "common", "description", "plate", "position")

# The UP and DN replicates must be input in correct order
# and align with each other well (i.e., A_condition_UP <-> A_condition_DOWN)
reps_up <- unique(df_map$UP)
reps_dn <- unique(df_map$DN)

# Check if the lengths are the same
# If the condition is FALSE, R will raise an error and stop execution
stopifnot(length(reps_up) == length(reps_dn))

# Iterate the process through all unique replicates
for (i in 1:length(reps_up))  {
  # Retrieve the UP and DN LANE information for each replicate
  rep_up <- reps_up[[i]]
  rep_dn <- reps_dn[[i]]
  rep_up_lanes <- df_map[df_map$UP == rep_up, "UP_LANE"]
  rep_dn_lanes <- df_map[df_map$DN == rep_dn, "DN_LANE"]
  # Assign the sum of all lanes to be the replicate column in df_sum
  df_sum[[rep_up]] <- rowSums(df_reads[, rep_up_lanes])
  df_sum[[rep_dn]] <- rowSums(df_reads[, rep_dn_lanes])
}

# Optional: run the line below to save a raw version before filtering strains
#write.table(df_sum, file = save_sum_raw_dir, row.names = FALSE, sep='\t')

# Get the sum of reads across all replicates, and filter by the cutoff
df_sum$sum <- rowSums(df_sum[, c(reps_up, reps_dn)], na.rm=TRUE) 

# Data pruning to remove strains with very low/inconsistent read counts
df_sum_f <- df_sum[df_sum$sum > low_sum_reads_filter, ]

# If enabled, remove duplicates and only keep the unique strain with max sum
if (bool_remove_duplicates) {
  df_sum_f <- remove_duplicates(df_sum_f, name_col='orf19', sum_col='sum')
}
df_sum_f <- df_sum_f[order(df_sum_f$feature), ] # Sort by feature name

# Export the final txt file
df_sum_f <- df_sum_f[, !names(df_sum_f) %in% "sum"]
write.table(df_sum_f, file = save_sum_final_dir, row.names = FALSE, sep='\t')