# Begin script - this script is only for generating QC stats and plots

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
bool_qc <- TRUE # TRUE if we want to generate the QC plots
nodox_qc_filter <- 50 # QC threshold - for each mutant, if its NO_DOX read is lower than this value, it will be replaced as NA to be filtered form downstream analysis

# Read the metatable that guides the column names to refer to
meta_dir = "./input/metatable_example.txt"
df_meta <- read.csv(file=meta_dir, header=TRUE, sep='\t')

# Note: this QC step does not require any plate or strain information

# Start the pipeline by pre-processing the input data to get LFC
df_all <- read.csv("./input/protocol_example_data.tsv", sep='\t')

# Pre-define the column names of the LFC values based on input, UP/DN separately
groups <- unique(df_meta$condition)
column_groups_up <- c()
column_groups_dn <- c()
for (group in groups) {
  df_meta_group <- df_meta[df_meta$condition==group, ]
  df_meta_group_nodox <- df_meta_group[df_meta_group$DOX==FALSE, ]
  condition_no_dox_up <- df_meta_group_nodox$UP
  condition_no_dox_dn <- df_meta_group_nodox$DN
  column_groups_up[[group]] <- paste0(condition_no_dox_up, "_LFC")
  column_groups_dn[[group]] <- paste0(condition_no_dox_dn, "_LFC")
}

################################################################################
# 1. Individual condition evaluation
# For one sample moderated t-test on each condition, integrating UP and DN tags
################################################################################

# Set QC directory
qc_output_dir = "./output/qc/"
if (!dir.exists(qc_output_dir)) { dir.create(qc_output_dir, recursive = TRUE) }

# Initialize lists to collect LFC dataframes for conditions (only needed for QC)
up_df_lfc_list <- c()
dn_df_lfc_list <- c()

# Iterate the scoring process through conditions
for (group in groups) {
  print(paste("Condition:",group))
  df_meta_group <- df_meta[df_meta$condition==group, ]
  df_meta_group_dox <- df_meta_group[df_meta_group$DOX==TRUE, ]
  df_meta_group_nodox <- df_meta_group[df_meta_group$DOX==FALSE, ]
  conditions_dox_list <- c(df_meta_group_dox$UP, df_meta_group_dox$DN)
  conditions_no_dox_list <- c(df_meta_group_nodox$UP, df_meta_group_nodox$DN)
  
  # For each group (condition), process LFC values and score them
  df_lfc <- data_preprocess_get_lfc(df_all, conditions_dox_list, conditions_no_dox_list, no_dox_cutoff=nodox_qc_filter)
  
  # # Plot saturation plot for baseline condition (YNB)
  # # Note: If needed, please adjust the gap parameters based on your need
  # if (group=="YNB") {
  #   print("Saturation plot for the baseline condition (YNB)")
  #   replicates_UP <- group_list(column_groups_up$YNB, subgroup_size=3)
  #   replicates_DN <- group_list(column_groups_dn$YNB, subgroup_size=3)
  #   
  #   sorted_hits_rep_len_up <- saturation_plot_baseline_replicates(df_lfc, replicates_UP,
  #                                       output_direc=paste0(qc_output_dir, "YNB_UP_LFC_replicates_saturation_plot.pdf"),
  #                                       gene_col="orf19", lfc_cutoff=2, baseline_cond=group)
  #   sorted_hits_rep_len_dn <- saturation_plot_baseline_replicates(df_lfc, replicates_DN,
  #                                       output_direc=paste0(qc_output_dir, "YNB_DN_LFC_replicates_saturation_plot.pdf"),
  #                                       gene_col="orf19", lfc_cutoff=2, baseline_cond=group)
  #   iteration <- c(1:length(sorted_hits_rep_len_up))
  #   pdf(paste0(qc_output_dir, "YNB_LFC_replicates_saturation_plot.pdf"))  # Save the plot into a PDF file
  #   gap.plot(iteration, sorted_hits_rep_len_up, type = "b", pch = 1, 
  #            gap = c(5, 210), brw = 0.02,
  #            ytics=c(0, 220, 230, 240, 250), xtics=seq(0,max(iteration),by=1),
  #        xlim = c(0, max(iteration)),
  #        ylim = c(0, max(c(sorted_hits_rep_len_up, sorted_hits_rep_len_dn))),
  #        xlab = paste0(group, " replicates"), ylab = "Number of total hits (LFC >= 2)",
  #        main = "Number of total hits given new replicates (add one at a time)")
  #   #lines(iteration, sorted_hits_rep_len_dn, type = "b", pch = 2)
  #   gap.plot(iteration, sorted_hits_rep_len_dn, type = "b", pch = 2, 
  #            gap = c(5, 210), brw = 0.02, add=TRUE,
  #            ytics=c(0, 220, 230, 240, 250), xtics=seq(0,max(iteration),by=1),
  #            xlim = c(0, max(iteration)),
  #            ylim = c(0, max(c(sorted_hits_rep_len_up, sorted_hits_rep_len_dn))))
  #   axis.break(2, 5, breakcol="snow", style="gap")
  #   axis.break(2, 5.5, breakcol="black", style="slash")
  #   axis.break(4, 5.5, breakcol="black", style="slash")
  #   legend("topleft", legend = c("UP tag", "DN tag"), pch=c(1,2))
  #   dev.off() # Close the PDF device
  # }
  
  # Generate QC plots (PCC heatmap & scatter plots)
  if (isTRUE(bool_qc)) {
    print("Generating PCC heatmap & scatter plots for the condition")
    df_lfc_updn <- generate_qc_plots(df_lfc, lfc_up_list=paste0(df_meta_group_nodox$UP, "_LFC"),
                                     lfc_dn_list=paste0(df_meta_group_nodox$DN, "_LFC"),
                                     condition_name=group, output_direc=qc_output_dir, input_type="LFC", bk_min=0.5)
    up_df_lfc_list <- c(up_df_lfc_list, df_lfc_updn$lfc_up)
    dn_df_lfc_list <- c(dn_df_lfc_list, df_lfc_updn$lfc_dn)
  }
  
}

# Plot global LFC PCC heatmap
if (isTRUE(bool_qc)) {
  print("Generating global LFC PCC heatmap & Within / Between correlation analysis")
  df_cond_up_all <- do.call(cbind, up_df_lfc_list)
  df_cond_dn_all <- do.call(cbind, dn_df_lfc_list)
  
  # Plot global QC heatmap
  plot_qc_correlation(df_cond_up_all, paste0(qc_output_dir, "global_heatmap_UP"),
                      show_num=FALSE, font_size=6, input_type="LFC", bk_min=0.5)
  plot_qc_correlation(df_cond_dn_all, paste0(qc_output_dir, "global_heatmap_DN"),
                      show_num=FALSE, font_size=6, input_type="LFC", bk_min=0.5)
  
  # Within / Between correlation analysis
  # Define the column groups
  column_groups_up <- list()
  column_groups_dn <- list()
  groups <- unique(df_meta$condition)
  for (group in groups) {
    print(paste("Condition:",group))
    df_meta_group <- df_meta[df_meta$condition==group, ]
    df_meta_group_nodox <- df_meta_group[df_meta_group$DOX==FALSE, ]
    lfc_up_list <- paste0(df_meta_group_nodox$UP, "_LFC")
    lfc_dn_list <- paste0(df_meta_group_nodox$DN, "_LFC")
    column_groups_up <- c(column_groups_up, list(lfc_up_list))
    column_groups_dn <- c(column_groups_dn, list(lfc_dn_list))
  }
  names(column_groups_up) <- groups
  names(column_groups_dn) <- groups
  
  # column groups defined at the top
  df_up_wb <- calculate_group_correlation(df_cond_up_all, column_groups_up)
  df_dn_wb <- calculate_group_correlation(df_cond_dn_all, column_groups_dn)
  
  bk <- seq(0.5, 1.0, by = 0.01)
  
  pdf(paste0(qc_output_dir, "wb_up.pdf"))
  pheatmap(df_up_wb, border_color = "white", breaks = bk,
           color=colorRampPalette(rev(brewer.pal(n=7,name= "YlGnBu")))(50),
           display_numbers = TRUE, cluster_rows = FALSE, 
           cluster_cols = FALSE)
  dev.off()
  
  pdf(paste0(qc_output_dir, "wb_dn.pdf"))
  pheatmap(df_dn_wb, border_color = "white", breaks = bk,
           color=colorRampPalette(rev(brewer.pal(n=7,name= "YlGnBu")))(50),
           display_numbers = TRUE, cluster_rows = FALSE, 
           cluster_cols = FALSE)
  dev.off()
}



################################################################################
# 2. Condition interaction effect compared to the reference standard
# For one sample moderated t-test on the difference between each condition and reference
# integrating UP and DN tags
################################################################################

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
  print(paste("Condition:",group))
  df_meta_group <- df_meta[df_meta$condition==group, ]
  df_meta_group_dox <- df_meta_group[df_meta_group$DOX==TRUE, ]
  df_meta_group_nodox <- df_meta_group[df_meta_group$DOX==FALSE, ]
  conditions_dox_list <- c(df_meta_group_dox$UP, df_meta_group_dox$DN)
  conditions_no_dox_list <- c(df_meta_group_nodox$UP, df_meta_group_nodox$DN)
  
  # For each group (condition), process LFC values and score them
  df_lfc <- data_preprocess_get_lfc(df_all, conditions_dox_list, conditions_no_dox_list, no_dox_cutoff=nodox_qc_filter)
  
  # Generate QC plots for dLFC (PCC heatmap & scatter plots)
  if (isTRUE(bool_qc)) {
    df_dlfc_updn <- generate_qc_plots_dLFC(df_lfc, df_lfc_ref, lfc_up_list=paste0(df_meta_group_nodox$UP, "_LFC"),
                                           lfc_dn_list=paste0(df_meta_group_nodox$DN, "_LFC"),
                                           ref_lfc_up_list=REF_UP_LFC, ref_lfc_dn_list=REF_DN_LFC, 
                                           condition_name=paste0(group, "_vs_REF"), output_direc=qc_output_dir, 
                                           bk_min=0.0, input_type="dLFC")
    up_df_dlfc_list <- c(up_df_dlfc_list, df_dlfc_updn$dlfc_up)
    dn_df_dlfc_list <- c(dn_df_dlfc_list, df_dlfc_updn$dlfc_dn)
  }
}

# Plot global dLFC PCC heatmap
if (isTRUE(bool_qc)) {
  df_dLFC_cond_up_all <- do.call(cbind, up_df_dlfc_list)
  df_dLFC_cond_dn_all <- do.call(cbind, dn_df_dlfc_list)
  plot_qc_correlation(df_dLFC_cond_up_all, paste0(qc_output_dir, "global_heatmap_UP"),
                      show_num=FALSE, font_size=6, bk_min=0.0, color_interval=100, input_type="dLFC")
  plot_qc_correlation(df_dLFC_cond_dn_all, paste0(qc_output_dir, "global_heatmap_DN"),
                      show_num=FALSE, font_size=6, bk_min=0.0, color_interval=100, input_type="dLFC")
}