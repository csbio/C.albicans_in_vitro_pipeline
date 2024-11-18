# Begin script - this script is for scoring with moderated t-test

# Set working directory
work_dir <- "C.albicans_in_vitro_pipeline/STAR_Protocols_code_data/score_example/"
setwd(work_dir)

# Source utility functions
# - Might get a negligible warning 
# - because BiocManager modifies the repository settings for Bioconductor and 
# - CRAN packages, but it gets overridden by the user’s default CRAN repository 
# - (set by getOption("repos")).
source("utils.R")

########################################
# PARAMETER SETTING & DATA PREPROCESSING
########################################

# Set important parameters
bool_compare_median_or_not <- TRUE # TRUE if we want to compare the LFCs to the median of mean LFC; FALSE if we want to compare to 0
bool_normalize_variance <- TRUE # TRUE if we want to apply variance normalization to the LFC values before doing the moderated t-test
nodox_qc_filter <- 50 # QC threshold - for each mutant, if its NO_DOX read is lower than this value, it will be replaced as NA to be filtered form downstream analysis

# Read the metatable that guides the column names to refer to
meta_dir = "./input/metatable.txt"
df_meta <- read.csv(file=meta_dir, header=TRUE, sep='\t')

# Get the column names from the meta table
orf19_col = unique(df_meta$orf19_column)[1]
plate_col = unique(df_meta$plate_column)[1]
feature_col = unique(df_meta$feature_column)[1]
common_col = unique(df_meta$common_column)[1]

# Start the pipeline by pre-processing the input data to get LFC
read_dir = "./input/sum_data.txt"
df_all <- read.csv(file=read_dir, header=TRUE, sep='\t')

################################################################################
# 1. Moderated t-test for single condition evaluation
################################################################################

# Conduct one sample moderated t-test on every condition, integrating UP and DN tags

# Set output directory
output_dir = "./output/output_scores/"
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

# Create dataframe to summarize all results for global comparison
# Add orf19/feature/common name to the output dataframe
# Users can choose to remove or add some of the columns based on own meta table
df_summary <- data.frame(orf19 = df_all[, orf19_col], plate = df_all[, plate_col], 
                         feature = df_all[, feature_col], common = df_all[, common_col], 
                         row.names = rownames(df_all))

# Iterate the scoring process through conditions
groups <- unique(df_meta$condition)
for (group in groups) {
  print(paste("Condition:", group))
  group_dir = paste0(output_dir, group, '/')
  if (!dir.exists(group_dir)) { dir.create(group_dir, recursive = TRUE) }
  
  df_meta_group <- df_meta[df_meta$condition==group, ]
  df_meta_group_dox <- df_meta_group[df_meta_group$DOX==TRUE, ]
  df_meta_group_nodox <- df_meta_group[df_meta_group$DOX==FALSE, ]
  conditions_dox_list <- c(df_meta_group_dox$UP, df_meta_group_dox$DN)
  conditions_no_dox_list <- c(df_meta_group_nodox$UP, df_meta_group_nodox$DN)
  
  # For each group (condition), process LFC values and score them
  df_lfc <- data_preprocess_get_lfc(df_all, conditions_dox_list, conditions_no_dox_list, no_dox_cutoff=nodox_qc_filter)
  
  # Store the LFC data before variance normalization
  write.table(df_lfc, paste0(group_dir, group, "_LFC.txt"),
              sep='\t', row.names = FALSE, quote = FALSE)
  
  cond_output_dir <- paste0(group_dir, group, "_")
  
  # Case 1: Compare mean LFC of every mutant to 0
  # Only need to copy commands as below but change the parameter bool_compare_median_or_not to FALSE at the beginning of this script
  # Then define new output directory with correct names
  
  # Case 2: Compare mean LFC of every mutant to the median of mean LFC values for all mutants
  df_modt_output <- get_mod_t_test_updn(df_lfc, lfc_up_list=paste0(df_meta_group_nodox$UP, "_LFC"),
                                        lfc_dn_list=paste0(df_meta_group_nodox$DN, "_LFC"),
                                        condition_name=group, num_valid_rep=3,
                                        compare_to_median_of_meanLFC=bool_compare_median_or_not,
                                        norm_var=bool_normalize_variance, output_direc=cond_output_dir)
  
  # Concatenate all results for global comparison
  # Removed strains that were not in the test will also be included but displayed as NAs in mean and fdr
  df_summary[rownames(df_modt_output), paste0(group, "_mean")] <- df_modt_output$mean
  df_summary[rownames(df_modt_output), paste0(group, "_fdr")] <- df_modt_output$fdr
}

# Output the summarized moderated t-test results
# Can change the output directory below and make sure we mark correctly (compare to 0 or median of means)
summary_direc <- paste(output_dir, "summary_mod_t_test_compare_median_normvar.csv", sep="")
print(paste("Saved final LFC scores to", summary_direc, sep=" "))
write.csv(df_summary, summary_direc)



################################################################################
# 2. Moderated t-test for condition interaction effect compared to reference standard
################################################################################

# Conduct one sample moderated t-test on the difference between every condition and reference, integrating UP and DN tags

# Define output directory:
dLFC_output_dir <- "./output/output_each_vs_REF/"
if (!dir.exists(dLFC_output_dir)) { dir.create(dLFC_output_dir, recursive = TRUE) }

# Create dataframe to summarize all results for global comparison
# Add orf19/feature/common name to the output dataframe
# Users can choose to remove or add some of the columns based on own meta table
df_summary_dLFC <- data.frame(orf19 = df_all[, orf19_col], plate = df_all[, plate_col], 
                              feature = df_all[, feature_col], common = df_all[, common_col], 
                              row.names = rownames(df_all))

# Initialize lists to collect dLFC dataframes for conditions (only needed for QC)
up_df_dlfc_list <- c()
dn_df_dlfc_list <- c()

# Utilize all the biological replicates (standard YNB) to form the reference tag lists and provide more precise reference
ref_group_name = 'YNB'
df_meta_ref_group <- df_meta[df_meta$condition==ref_group_name, ]
df_meta_ref_dox <- df_meta_ref_group[df_meta_ref_group$DOX==TRUE, ]
df_meta_ref_nodox <- df_meta_ref_group[df_meta_ref_group$DOX==FALSE, ]
ref_dox_list <- c(df_meta_ref_dox$UP, df_meta_ref_dox$DN)
ref_no_dox_list <- c(df_meta_ref_nodox$UP, df_meta_ref_nodox$DN)
df_lfc_ref <- data_preprocess_get_lfc(df_all, ref_dox_list, ref_no_dox_list, no_dox_cutoff=nodox_qc_filter)

REF_UP_LFC <- paste0(df_meta_ref_nodox$UP, "_LFC")
REF_DN_LFC <- paste0(df_meta_ref_nodox$DN, "_LFC")

# Iterate the scoring process through conditions
groups <- unique(df_meta$condition)
dLFC_groups <- groups[groups != ref_group_name]
for (group in dLFC_groups) {
  print(paste("Condition:", group))
  group_dir = paste0(dLFC_output_dir, group, '/')
  if (!dir.exists(group_dir)) { dir.create(group_dir, recursive = TRUE) }
  
  df_meta_group <- df_meta[df_meta$condition==group, ]
  df_meta_group_dox <- df_meta_group[df_meta_group$DOX==TRUE, ]
  df_meta_group_nodox <- df_meta_group[df_meta_group$DOX==FALSE, ]
  conditions_dox_list <- c(df_meta_group_dox$UP, df_meta_group_dox$DN)
  conditions_no_dox_list <- c(df_meta_group_nodox$UP, df_meta_group_nodox$DN)
  
  # For each group (condition), process LFC values and score them
  df_lfc <- data_preprocess_get_lfc(df_all, conditions_dox_list, conditions_no_dox_list, no_dox_cutoff=nodox_qc_filter)
  
  # Get the moderated t-test output for dLFC
  cond_dLFC_output_dir <- paste0(group_dir, group, "_")
  df_dLFC_modt_output <- get_dLFC_mod_t_test(df_lfc, df_lfc_ref, lfc_up_list=paste0(df_meta_group_nodox$UP, "_LFC"),
                                             lfc_dn_list=paste0(df_meta_group_nodox$DN, "_LFC"),
                                             ref_lfc_up_list=REF_UP_LFC, ref_lfc_dn_list=REF_DN_LFC,
                                             condition_name=paste0(group, "_vs_REF"), num_valid_rep=3,
                                             norm_var=bool_normalize_variance,
                                             output_direc=cond_dLFC_output_dir)
  
  # Concatenate all results for global comparison
  # Removed strains that were not in the test will also be included but displayed as NAs in mean and fdr
  df_summary_dLFC[rownames(df_dLFC_modt_output), paste0(group, "_mean")] <- df_dLFC_modt_output$mean
  df_summary_dLFC[rownames(df_dLFC_modt_output), paste0(group, "_fdr")] <- df_dLFC_modt_output$fdr
}

# Output the summarized moderated t-test results
dLFC_summary_direc <- paste(dLFC_output_dir, "summary_mod_t_test_dLFC_normvar.csv", sep="")
print(paste("Saved final dLFC scores to", dLFC_summary_direc, sep=" "))
write.csv(df_summary_dLFC, dLFC_summary_direc)