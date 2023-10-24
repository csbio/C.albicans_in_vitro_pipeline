######
# UTILITY FUNCTIONS
######


#' Log-normalizes reads.
#' 
#' Log2-normalizes reads for a column with a given pseudocount and scaling factor
#' 
#' @param df List of read counts.
#' @param cf1 Scaling factor (default 1e6).
#' @param cf2 Pseudocount (default 1).
#' @return Log- and depth-normalized read counts.
#' @export
normalize_reads <- function(df, cf1 = 1e6, cf2 = 1) {
  return(log2((df / sum(df, na.rm = TRUE)) * cf1 + cf2))
}


#' Normalize the LFC values across replicates.
#' 
#' Assume consistent variance for mice replicates (same width of distribution)
#' αr: median(a1, a2, a3,...an) - median of the standard deviation from all mice's sd
#' αm: standard deviation of a single mouse column
#' LFC^ =  LFC * αr / αm
#' Effect: normalize achieve consistent standard deviation across replicates
#' 
#' @param df dataframe of LFC values.
#' @param name_LFC_list List with input replicate column names (usually 1 tag (UP or DN) only).    e.g. c('A_FBS_UP_LFC', 'B_FBS_UP_LFC',  'C_FBS_UP_LFC')
#' @return variance-normalized LFC dataframe for a tag group
#' @export
normalize_variance <- function(df, name_LFC_list) {
  df_tag <- df[name_LFC_list]
  am_list <- c()
  for (tag in name_LFC_list) {
    am <- sd(df_tag[,tag], na.rm=TRUE)
    am_list <- c(am_list, am)
  }
  ar <- median(am_list, na.rm=TRUE)
  for (tag in name_LFC_list) {
    am <- sd(df_tag[,tag], na.rm=TRUE)
    df_tag[,tag] <- df_tag[,tag] * ar / am
  }
  return(df_tag)
}


#' Pre-process the input dataframe for every condition into LFC data,
#' including assigning NAs to entries with NO DOX reads < 50, and Log2-normalizes 
#' reads with a given pseudocount and scaling factor
#' 
#' @param df_all Dataframe of the raw read counts.
#' @param conditions_dox_list List of dox column names.
#' @param conditions_no_dox_list List of no-dox column names.
#' @param no_dox_cutoff No-dox read count cutoff under which the value will be replaced with NA in the table.
#' @return Pre-processed and normalized LFC data.
#' @export
data_preprocess_get_lfc <- function(df_all, conditions_dox_list, conditions_no_dox_list, no_dox_cutoff = 50) {
  df_all_nodox <- df_all[conditions_no_dox_list]
  df_all_dox <- df_all[conditions_dox_list]
  
  # Assign NAs to where No DOX reads < 50
  df_all_nodox[df_all_nodox<no_dox_cutoff] <- NA
  
  # Normalize reads
  for (col in conditions_no_dox_list) {
    df_all_nodox[,col] <- normalize_reads(df_all_nodox[,col])
  }
  for (col in conditions_dox_list) {
    df_all_dox[,col] <- normalize_reads(df_all_dox[,col])
  }
  
  # Calculate Log2 fold change
  df_lfc <- df_all_nodox - df_all_dox
  colnames(df_lfc) <- paste(conditions_no_dox_list, "_LFC", sep="")
  df_lfc$orf19 <- df_all$orf19
  df_lfc$plate <- df_all$plate
  return(df_lfc)
}


#' Function to subtract the median of means from LFC values of UP tag or DN tag data only
#' 
#' Goal: Compare the mean LFC of a mutant to the median of the mean LFC values of all mutants via moderated t-test
#' Technically that median of means is an estimate with some variance associated, 
#' but given how many values we are estimating from, that variance is negligible relative to the variance in our LFCs. 
#' So, in practice, a one sample t-test should be fine.
#'
#' @param df_cond: Dataframe that contains three columns of LFC values for a single condition in a single tag (e.g. df_cond_up) 
#' @return Dataframe in the same format with median of means subtracted
#' @export
subtract_median_from_LFC <- function(df_cond) {
  median_of_LFCmeans <- median(rowMeans(df_cond), na.rm=TRUE)
  df_cond <- df_cond - median_of_LFCmeans
  return(df_cond)
}


#' Function to generate the QC PCC heatmaps for one input dataframe with replicates' LFC values
#'
#' @param df: Dataframe that contains all replicates' LFC values for a certain conditions    e.g. df_cond_up
#' @param filename: Name of the file based on the condition and tag    e.g. FBS_UP
#' @return None (will generate a PCC heatmap for an input dataframe)
#' @export
plot_qc_correlation <- function(df, filename, show_num=TRUE, font_size=10, 
                                bk_min=0, color_interval=50, input_type="LFC") {
  # Change the label from LFC to dLFC if given input_type as "dLFC"
  if (!(input_type == "LFC")) {
    columns_to_plot <- sub("LFC", input_type, colnames(df))
    colnames(df) <- columns_to_plot
  }
  
  # Plot the heatmap for Pearson Correlation coefficients among the replicates
  cor_matrix <- cor(df, use = "pairwise.complete.obs", method = "pearson")
  heatmap_file <- paste0(filename, "_", input_type, "_PCC_heatmap", ".pdf")
  
  # Define a color palette for the heatmap
  #n_colors <- 100  # Number of colors
  #color_palette <- colorRampPalette(brewer.pal(8, "Blues"))
  
  pdf(heatmap_file)
  
  # For cases with high PCCs, we raise the lower bound of plotting
  bk <- seq(bk_min, 1.0, by = 0.01)
  #bk <- c(seq(0,0.49,by=0.01), seq(0.5,1,by=0.01))
  
  pheatmap(cor_matrix, border_color = "white", breaks = bk, 
           #color=colorRampPalette(rev(brewer.pal(n=7,name= "Blues")))(50),
           color=colorRampPalette(rev(brewer.pal(n=7,name= "YlGnBu")))(color_interval),
           display_numbers = show_num, cluster_rows = FALSE, 
           cluster_cols = FALSE, fontsize = font_size)
  dev.off()
}


#' Function to generate the QC scatter plots for one input dataframe with replicates' LFC values
#'
#' @param df: Dataframe that contains all replicates' LFC values for a certain conditions    e.g. df_cond_up
#' @param filename: Name of the file based on the condition and tag    e.g. FBS_UP
#' @return None (will generate pairwise scatter plots for an input dataframe)
#' @export
plot_qc_scatter <- function(df, filename, input_type="LFC") {
  # Set the file format
  file_format <- "pdf"
  
  # Set the column names to be plotted
  columns_to_plot <- colnames(df) 
  
  # Change the label from LFC to dLFC if given input_type as "dLFC"
  if (!(input_type == "LFC")) {
    columns_to_plot <- sub("LFC", input_type, columns_to_plot) 
    colnames(df) <- columns_to_plot
  }
  
  # Generate pairwise scatter plots and save them to files
  for (i in 1:length(columns_to_plot)) {
    for (j in 1:length(columns_to_plot)) {
      if (i == j) next
      
      # Set the file name
      file_name <- paste0(filename, columns_to_plot[i], "_vs_", columns_to_plot[j], ".", file_format)
      
      # Calculate the correlation coefficient
      x <- df[, columns_to_plot[i]]
      y <- df[, columns_to_plot[j]]
      cor_coef <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
      
      # Store file
      pdf(file_name)
      
      # Plot the scatter plot for the current pair of columns
      plot(x, y, xlab = columns_to_plot[i], ylab = columns_to_plot[j],
           main=paste("Pearson Correlation:", round(cor_coef, 2)))
      
      # Close the PDF or PNG file
      dev.off()
    }
  }
}


#' Function to generate the LFC QC plots for each condition, UP & DN tags separately
#'
#' @param df_lfc: Dataframe that contains all LFC values for all conditions    e.g. df_lfc
#' @param lfc_up_list: List of UP tag columns names of a certain condition (length must be 3 in the order of A, B, C)    e.g. c('A_FBS_UP_LFC', 'B_FBS_UP_LFC',  'C_FBS_UP_LFC')
#' @param lfc_dn_list: List of DN tag columns names of a certain condition (length must be 3 in the order of A, B, C)   e.g. c('A_FBS_DN_LFC', 'B_FBS_DN_LFC',  'C_FBS_DN_LFC')
#' @param condition_name: Name of condition    e.g. "FBS"
#' @param output_direc: Directory for output files    e.g. "/Users/qc/"
#' @returna list containing one dataframe for UP and one for DN under this condition
#' @export
generate_qc_plots <- function(df_lfc, lfc_up_list, lfc_dn_list, condition_name, 
                              output_direc, input_type="LFC", bk_min=0) {
  heatmap_dir <- paste0(output_direc, input_type, "_heatmap/")
  scatter_dir <- paste0(output_direc, input_type, "_scatter/")
  if (!dir.exists(heatmap_dir)) { dir.create(heatmap_dir, recursive = TRUE) }
  if (!dir.exists(scatter_dir)) { dir.create(scatter_dir, recursive = TRUE) }
  
  # Take the subset columns for a single condition
  # Seperate UP & DN tags
  df_cond_up <- df_lfc[lfc_up_list]
  df_cond_dn <- df_lfc[lfc_dn_list]
  
  # Plot heatmap
  plot_qc_correlation(df_cond_up, paste0(heatmap_dir, condition_name, "_UP"), bk_min=bk_min, input_type=input_type)
  plot_qc_correlation(df_cond_dn, paste0(heatmap_dir, condition_name, "_DN"), bk_min=bk_min, input_type=input_type)
  
  # Plot scatter
  plot_qc_scatter(df_cond_up, scatter_dir, input_type=input_type)
  plot_qc_scatter(df_cond_dn, scatter_dir, input_type=input_type)
  
  return(list("lfc_up" = df_cond_up, "lfc_dn" = df_cond_dn))
}


#' Function to generate the dLFC QC plots for each condition, UP & DN tags separately
#'
#' @param df_lfc: Dataframe that contains all LFC values for all conditions    e.g. df_lfc
#' @param lfc_up_list: List of UP tag columns names of a certain condition
#' @param lfc_dn_list: List of DN tag columns names of a certain condition
#' @param ref_lfc_up_list: List of UP tag columns names of the reference condition (no length requirement as we will take the average)    e.g. c('A_CONT_UP_LFC', 'A_CONT_UP_LFC',  'A_CONT_UP_LFC')
#' @param ref_lfc_dn_list: List of DN tag columns names of the reference condition (no length requirement as we will take the average)   e.g. c('A_CONT_DN_LFC', 'A_CONT_DN_LFC',  'A_CONT_DN_LFC')
#' @param condition_name: Name of condition    e.g. "FBS_vs_REF"
#' @param output_direc: Directory for output files    e.g. "/Users/qc/"
#' @return a list containing one dataframe for UP and one for DN under this condition
#' @export
generate_qc_plots_dLFC <- function(df_lfc, df_lfc_ref, lfc_up_list, lfc_dn_list,
                                   ref_lfc_up_list, ref_lfc_dn_list, 
                                   condition_name, output_direc, 
                                   bk_min=0.5, input_type="dLFC") {
  heatmap_dir <- paste0(output_direc, input_type, "_heatmap/")
  scatter_dir <- paste0(output_direc, input_type, "_scatter/")
  if (!dir.exists(heatmap_dir)) { dir.create(heatmap_dir, recursive = TRUE) }
  if (!dir.exists(scatter_dir)) { dir.create(scatter_dir, recursive = TRUE) }
  
  # Take the subset columns for a single condition
  # Seperate UP & DN tags at first for possible median subtraction (will merge later)
  df_cond_up <- df_lfc[lfc_up_list]
  df_cond_dn <- df_lfc[lfc_dn_list]
  df_ref_up <- df_lfc_ref[ref_lfc_up_list]
  df_ref_dn <- df_lfc_ref[ref_lfc_dn_list]
  
  # For each tag, calculate the differential LFC values by subtracting the mean control LFCs from the three condition replicates separately
  # Assumption: each row is corresponding to the same mutant in both datasets
  # Note: it is not an A-A, B-B, C-C relationship
  df_cond_dLFC_up <- df_cond_up - rowMeans(df_ref_up, na.rm=TRUE)
  df_cond_dLFC_dn <- df_cond_dn - rowMeans(df_ref_dn, na.rm=TRUE)
  
  # Plot heatmap
  plot_qc_correlation(df_cond_dLFC_up, paste0(heatmap_dir, condition_name, "_UP"), bk_min=bk_min, color_interval=100, input_type=input_type)
  plot_qc_correlation(df_cond_dLFC_dn, paste0(heatmap_dir, condition_name, "_DN"), bk_min=bk_min, color_interval=100, input_type=input_type)
  
  # Plot scatter
  plot_qc_scatter(df_cond_dLFC_up, scatter_dir, input_type=input_type)
  plot_qc_scatter(df_cond_dLFC_dn, scatter_dir, input_type=input_type)
  
  return(list("dlfc_up" = df_cond_dLFC_up, "dlfc_dn" = df_cond_dLFC_dn))
}


#' Function to conduct one sample moderated t-test on the LFC values for a single condition, integrating UP & DN tags
#'
#' @param df_lfc: Dataframe that contains all LFC values for all conditions    e.g. df_lfc
#' @param lfc_up_list: List of UP tag columns names of a certain condition
#' @param lfc_dn_list: List of DN tag columns names of a certain condition
#' @param condition_name: Name of condition    e.g. "FBS"
#' @param compare_to_median_of_meanLFC: Whether we compare the mean to the median of the mean LFC values of all mutants, or to 0 (e.g. TRUE)
#' @param norm_var: Whether we apply variance normalization to the LFC values before doing the moderated t-test (e.g. TRUE)
#' @param output_direc: Directory for output files    e.g. "/Users/mod_t_test_results/"
#' @return Moderated t-test results for a single condition
#' @export
get_mod_t_test_updn <- function(df_lfc, lfc_up_list, lfc_dn_list, condition_name, 
                                num_valid_rep=3, compare_to_median_of_meanLFC=TRUE, norm_var=TRUE, output_direc="mod_t_test_results/"){
  # Create the output directory if it does not exist
  # Create a folder to store strains that are not qualified to enter the statistical test
  remove_direc <- paste(output_direc, "removed_strains/", sep="")
  if (!dir.exists(remove_direc)) { dir.create(remove_direc, recursive = TRUE) }
  # Create a folder to store results from the moderated t-test for each condition
  output_direc <- paste(output_direc, "results_separate_conditions/", sep="")
  if (!dir.exists(output_direc)) { dir.create(output_direc, recursive = TRUE) }
  
  # Take the subset columns for a single condition
  # Seperate UP & DN tags at first for possible median subtraction (will merge later)
  df_cond_up <- df_lfc[lfc_up_list]
  df_cond_dn <- df_lfc[lfc_dn_list]
  
  # Normalize the variance per tag if specified (TRUE by default)
  if (isTRUE(norm_var)){
    df_cond_up <- normalize_variance(df_cond_up, lfc_up_list)
    df_cond_dn <- normalize_variance(df_cond_dn, lfc_dn_list)
  }
  
  # Will enter the if condition if compare_to_median_of_meanLFC is TRUE, which means we want to compare the mean LFCs to the median of mean LFCs
  if (isTRUE(compare_to_median_of_meanLFC)) {
    df_cond_up <- subtract_median_from_LFC(df_cond_up)
    df_cond_dn <- subtract_median_from_LFC(df_cond_dn)
  } 
  
  # Merge UP & DN tags into a single dataframe
  # df_cond <- data.frame(df_cond_up[,1], df_cond_dn[,1], df_cond_up[,2], df_cond_dn[,2], df_cond_up[,3], df_cond_dn[,3])
  # colnames(df_cond) <- c(lfc_up_list[1], lfc_dn_list[1], lfc_up_list[2], lfc_dn_list[2], lfc_up_list[3], lfc_dn_list[3])
  # rownames(df_cond) <- df_lfc$orf19
  
  # Update merging method
  df_cond <- data.frame(df_cond_up[,1], df_cond_dn[,1])
  colnames(df_cond) <- c(lfc_up_list[1], lfc_dn_list[1])
  for (i in 2:length(lfc_up_list)) {
    df_cond[, lfc_up_list[i]] <- df_cond_up[,i]
    df_cond[, lfc_dn_list[i]] <- df_cond_dn[,i]
  }
  rownames(df_cond) <- df_lfc$orf19
  
  # Check number of NAs per mutant (must has >=3 reps), filter unqualified rows and record them
  sum_NA_list <- rowSums(is.na(df_cond))  # Will return a named num structure
  remove_list <- names(sum_NA_list[sum_NA_list>(length(lfc_up_list)*2-num_valid_rep)])
  print(paste("Removed strains with insufficient data from condition", condition_name, sep=" "))
  print(remove_list)
  df_cond$plate <- df_lfc$plate
  df_remove <- df_cond[(rownames(df_cond) %in% remove_list), ]
  df_remove$mean <- rowMeans(df_remove[1:(dim(df_remove)[2]-1)], na.rm=TRUE)
  remove_direc <- paste(remove_direc,'mod_t_test_',condition_name,'.csv',sep='')
  write.csv(df_remove, remove_direc, row.names = TRUE)
  
  # Take the filtered dataframe to the moderated t-test
  df_cond <- df_cond[!(rownames(df_cond) %in% remove_list), ]
  
  # Define block: A B C for UP and DN === c(1,1,2,2,3,3)
  # This means the input data should have columns as: A_UP, A_DN, B_UP, B_DN, C_UP, C_DN, ...
  brep <- rep(1:length(lfc_up_list), each = 2) 
  #print(brep)
  multipleHypothesisMethod <- 'BH'
  
  # Directly compare to 0, do one sample moderated t-test
  # The last column is "plate" which should not be included --- df_cond[1:(dim(df_cond)[2]-1)]
  input_matrix <- data.matrix(df_cond[1:(dim(df_cond)[2]-1)])
  
  dupcor <- duplicateCorrelation(input_matrix, design=NULL, block=brep)
  #print(dupcor$consensus.correlation)
  
  # moderated t-test
  fit <- lmFit(input_matrix, design=NULL, block=brep, correlation=dupcor$consensus)
  efit <- eBayes(fit)
  #topTable(efit, coef=1)
  df_cond["mean"] <- fit$Amean
  df_cond["pval"] <- efit$p.value
  df_cond["fdr"] <- p.adjust(efit$p.value, method = multipleHypothesisMethod)

  output_direc <- paste(output_direc, 'mod_t_test_', condition_name, '.csv', sep='')
  write.csv(df_cond, output_direc, row.names = TRUE)
  return(df_cond)
}


#' Function to conduct moderated t-test on the differential LFC values for a condition vs control, integrating UP & DN tags
#' Currently we calculate the differential LFC values as (condition_LFCs - mean(control_LFCs)), which will result in 3 residual values for the test.
#'
#' @param df_lfc: Dataframe that contains all LFC values for all conditions    e.g. df_lfc
#' @param lfc_up_list: List of UP tag columns names of a certain condition
#' @param lfc_dn_list: List of DN tag columns names of a certain condition
#' @param ref_lfc_up_list: List of UP tag columns names of the reference condition (no length requirement as we will take the average)    e.g. c('A_CONT_UP_LFC', 'A_CONT_UP_LFC',  'A_CONT_UP_LFC')
#' @param ref_lfc_dn_list: List of DN tag columns names of the reference condition (no length requirement as we will take the average)   e.g. c('A_CONT_DN_LFC', 'A_CONT_DN_LFC',  'A_CONT_DN_LFC')
#' @param condition_name: Name of condition    e.g. "FBS_vs_REF"
#' @param norm_var: Whether we apply variance normalization to the LFC values before doing the moderated t-test (e.g. TRUE)
#' @param output_direc: Directory for output files    e.g. "/Users/mod_t_test_dLFC_results/"
#' @return Moderated t-test results for a condition vs control
#' @export
get_dLFC_mod_t_test <- function(df_lfc, df_lfc_ref, lfc_up_list, lfc_dn_list, ref_lfc_up_list, ref_lfc_dn_list, condition_name,
                                num_valid_rep=3, norm_var=TRUE, output_direc="mod_t_test_dLFC_results/"){
  # Create the output directory if it does not exist
  # Create a folder to store strains that are not qualified to enter the statistical test
  remove_direc <- paste(output_direc, "removed_strains/", sep="")
  if (!dir.exists(remove_direc)) { dir.create(remove_direc, recursive = TRUE) }
  # Create a folder to store results from the moderated t-test for each condition
  output_direc <- paste(output_direc, "results_separate_conditions/", sep="")
  if (!dir.exists(output_direc)) { dir.create(output_direc, recursive = TRUE) }
  
  # Take the subset columns for a single condition
  # Seperate UP & DN tags at first for possible median subtraction (will merge later)
  df_cond_up <- df_lfc[lfc_up_list]
  df_cond_dn <- df_lfc[lfc_dn_list]
  df_ref_up <- df_lfc_ref[ref_lfc_up_list]
  df_ref_dn <- df_lfc_ref[ref_lfc_dn_list]
  
  # Normalize the variance per tag if specified (TRUE by default)
  if (isTRUE(norm_var)){
    df_cond_up <- normalize_variance(df_cond_up, lfc_up_list)
    df_cond_dn <- normalize_variance(df_cond_dn, lfc_dn_list)
    df_ref_up <- normalize_variance(df_ref_up, ref_lfc_up_list)
    df_ref_dn <- normalize_variance(df_ref_dn, ref_lfc_dn_list)
  }
  
  # For each tag, calculate the differential LFC values by subtracting the mean control LFCs from the three condition replicates separately
  # Assumption: each row is corresponding to the same mutant in both datasets
  # Note: it is not an A-A, B-B, C-C relationship
  df_cond_dLFC_up <- df_cond_up - rowMeans(df_ref_up, na.rm=TRUE)
  df_cond_dLFC_dn <- df_cond_dn - rowMeans(df_ref_dn, na.rm=TRUE)
  
  # Merge UP & DN tags into a single dataframe
  # df_cond <- data.frame(df_cond_dLFC_up[,1], df_cond_dLFC_dn[,1], df_cond_dLFC_up[,2], df_cond_dLFC_dn[,2], df_cond_dLFC_up[,3], df_cond_dLFC_dn[,3])
  # colnames(df_cond) <- c(lfc_up_list[1], lfc_dn_list[1], lfc_up_list[2], lfc_dn_list[2], lfc_up_list[3], lfc_dn_list[3])
  # rownames(df_cond) <- df_lfc$orf19
  
  # Update merging method
  df_cond <- data.frame(df_cond_dLFC_up[,1], df_cond_dLFC_dn[,1])
  colnames(df_cond) <- c(lfc_up_list[1], lfc_dn_list[1])
  for (i in 2:length(lfc_up_list)) {
    df_cond[, lfc_up_list[i]] <- df_cond_dLFC_up[,i]
    df_cond[, lfc_dn_list[i]] <- df_cond_dLFC_dn[,i]
  }
  rownames(df_cond) <- df_lfc$orf19
  
  # Check number of NAs per mutant (must has >=3 reps), filter unqualified rows and record them
  sum_NA_list <- rowSums(is.na(df_cond))  # Will return a named num structure
  remove_list <- names(sum_NA_list[sum_NA_list>(length(lfc_up_list)*2-num_valid_rep)])
  print(paste("Removed strains from condition", condition_name, sep=" "))
  print(remove_list)
  df_cond$plate <- df_lfc$plate
  df_remove <- df_cond[(rownames(df_cond) %in% remove_list), ]
  df_remove$mean <- rowMeans(df_remove[1:(dim(df_remove)[2]-1)], na.rm=TRUE)
  remove_direc <- paste(remove_direc,'mod_t_test_dLFC_',condition_name,'.csv',sep='')
  write.csv(df_remove, remove_direc, row.names = TRUE)
  
  # Take the filtered dataframe to the moderated t-test
  df_cond <- df_cond[!(rownames(df_cond) %in% remove_list), ]
  
  # Define block: A B C for UP and DN === c(1,1,2,2,3,3)
  # This means the input data should have columns as: A_UP, A_DN, B_UP, B_DN, C_UP, C_DN, ...
  brep <- rep(1:length(lfc_up_list), each = 2) 
  #print(brep)
  multipleHypothesisMethod <- 'BH'
  
  # Directly compare to 0, do one sample moderated t-test
  # The last column is "plate" which should not be included --- df_cond[1:(dim(df_cond)[2]-1)]
  input_matrix <- data.matrix(df_cond[1:(dim(df_cond)[2]-1)])
  
  dupcor <- duplicateCorrelation(input_matrix, design=NULL, block=brep)
  #print(dupcor$consensus.correlation)
  
  # moderated t-test
  fit <- lmFit(input_matrix, design=NULL, block=brep, correlation=dupcor$consensus)
  efit <- eBayes(fit)
  #topTable(efit, coef=1)
  df_cond["mean"] <- fit$Amean
  df_cond["pval"] <- efit$p.value
  df_cond["fdr"] <- p.adjust(efit$p.value, method = multipleHypothesisMethod)
  
  output_direc <- paste(output_direc, 'mod_t_test_dLFC_', condition_name, '.csv', sep='')
  write.csv(df_cond, output_direc, row.names = TRUE)
  return(df_cond)
}


#' Get the hit genes.
#' 
#' Get the hit genes given the effect size cutoff and FDR cutoff
#' 
#' @param df Dataframe of integrated moderated t-test output (each condition has an effect size column and FDR column)
#' @param lfc_col LFC column name
#' @param fdr_col FDR column name
#' @param lfc_cutoff LFC cutoff
#' @param fdr_cutoff FDR cutoff
#' @return None (plot and save the saturation plot)
#' @export
get_hits <- function(df, lfc_col, fdr_col, lfc_cutoff, fdr_cutoff, gene_col='orf19') {
  genes_LFC <- df[gene_col][df[lfc_col]>=lfc_cutoff]
  genes_FDR <- df[gene_col][df[fdr_col]<=fdr_cutoff]
  hits <- intersect(genes_LFC, genes_FDR)
  hits <- hits[!is.na(hits)]
  return(hits)
}


calculate_group_correlation <- function(df, column_groups) {
  # Initialize the correlation matrix
  num_groups <- length(column_groups)
  correlation_results <- matrix(NA, nrow = num_groups, ncol = num_groups)
  rownames(correlation_results) <- colnames(correlation_results) <- names(column_groups)
  
  # Calculate average correlation within each group and between groups
  for (i in 1:num_groups) {
    for (j in 1:num_groups) {
      if (i == j) {
        # Calculate average correlation within the group
        group_columns <- column_groups[[i]]
        group_data <- df[, group_columns]
        avg_cor <- mean(cor(group_data, use = "pairwise.complete.obs", method = "pearson"))
        correlation_results[i, j] <- avg_cor
      } else {
        # Calculate correlation between groups
        group1_columns <- column_groups[[i]]
        group2_columns <- column_groups[[j]]
        group1_data <- df[, group1_columns]
        group2_data <- df[, group2_columns]
        cor_val <- cor(group1_data, group2_data, use = "pairwise.complete.obs", method = "pearson")
        avg_cor <- mean(cor_val)
        correlation_results[i, j] <- avg_cor
      }
    }
  }
  
  # Return the correlation matrix
  return(correlation_results)
}


group_list <- function(input_list, subgroup_size=3) {
  grouped_list <- split(input_list, rep(1:(length(input_list) %/% subgroup_size), each = subgroup_size, length.out = length(input_list)))
  return(grouped_list)
}


# Apply greedy algorithm to sort the input dictionary of different list of mutants
collect_ordered_lists <- function(sorted_hits) {

  # Initialize the ordered lists with the longest list
  ordered_lists <- sorted_hits[1]
  common_values <- sorted_hits[[1]]
  
  # Initialize the remaining items as the sorted dictionary items
  remaining_items <- sorted_hits[-which(names(sorted_hits) == names(sorted_hits[1]))]
  
  # Iterate until there are no remaining items
  while (length(remaining_items) > 0) {
    # Find the item that adds the most unique elements to the ordered lists
    common_values_len_list <- c()
    for (i in 1:length(remaining_items)) {
      current_list <- remaining_items[[i]]
      new_values <- setdiff(current_list, common_values)
      new_common_values <- union(common_values, new_values)
      common_values_len_list <- c(common_values_len_list, length(new_common_values))
    }
    # Get the index of the largest element
    index_of_largest <- which.max(common_values_len_list)
    current_item <- remaining_items[index_of_largest]
    common_values <- union(common_values, current_item)
    ordered_lists <- c(ordered_lists, current_item)
    
    # Remove the processed item from the remaining items
    remaining_items <- remaining_items[-index_of_largest]#which(names(remaining_items) == names(current_item))]
  }
  return(ordered_lists)
}


#' Saturation plot at baseline replicate level (e.g., YNB).
#' 
#' This version is sorting the dictionary based on the new additional unique hits
#' Greedily plot how many new hits are being identified with new replicates added
#' Note that this is before the statistical test and serves as a check
#' We should merge the technical replicates (e.g., A/B/C) into the basic unit
#' 
#' @param df Dataframe of the LFC values for all baseline replicates (only LFC columns)
#' @param replicates All replicates' names
#' @param output_direc: Directory for output files    e.g. "/Users/qc/"
#' @param gene_col: Name of the column for all the genes
#' @param lfc_cutoff: LFC cutoff
#' @param fdr_cutoff: FDR cutoff
#' @param baseline_cond: Name of baseline condition
#' @return None (plot and save the saturation plot)
#' @export
saturation_plot_baseline_replicates <- function(df, replicates, output_direc="", 
                                                           gene_col="orf19", lfc_cutoff=2,
                                                           baseline_cond="YNB") {
  # Collect the hits for each replicate given an LFC cutoff
  keys <- list()
  values <- list()
  for (i in 1:length(replicates)) {
    lfc_col <- replicates[[i]]
    # Average the replicates
    average_lfc <- rowMeans(df[,lfc_col])
    indices <- which(unlist(average_lfc) >= lfc_cutoff)
    replicate_hits <- df[indices, gene_col]
    replicate_hits <- replicate_hits[!is.na(replicate_hits)]
    # print(lfc_col)
    # print(length(replicate_hits))
    keys <- c(keys, list(lfc_col))
    values <- c(values, list(replicate_hits))
  }
  # Create a dictionary-like structure using names()
  all_hits <- setNames(values, keys)
  print(length(all_hits))
  sorted_hits <- all_hits[order(lengths(all_hits), decreasing = TRUE)]
  sorted_hits <- collect_ordered_lists(sorted_hits)
  #print(sorted_hits)
  
  # Start from the longest list
  longest_list <- sorted_hits[[1]]
  common_values <- longest_list
  
  # Initialize vectors for storing data
  iteration <- c(1)
  common_values_length <- c(length(common_values))
  
  # Iterate over the remaining lists in descending order
  for (i in 2:length(sorted_hits)) {
    current_list <- sorted_hits[[i]]
    new_values <- setdiff(current_list, common_values)
    common_values <- union(common_values, new_values)
    
    # Print the common values in each iteration
    print(paste("Common values after", i-1, "iterations:", length(common_values)))
    
    # Store iteration number and length of common_values
    iteration <- c(iteration, i)
    common_values_length <- c(common_values_length, length(common_values))
  }
  # Plot the length of common_values versus iterations
  pdf(output_direc)  # Save the plot into a PDF file
  plot(iteration, common_values_length, type = "b", 
       xlab = paste0(baseline_cond, " replicates"), ylab = "Number of total hits",
       main = "Number of total hits given new replicates (add one at a time)",
       xlim = c(0, max(iteration)), ylim=c(0, max(common_values_length)))
  dev.off() # Close the PDF device
  return(common_values_length)
}
