library(pheatmap)
library(RColorBrewer)
library(plotrix)

# Begin script - this script is only for generating QC stats and plots

# Set working directory
work_dir = "C.albicans_in_vitro_pipeline/src/score"
setwd(work_dir)

# Source utility functions
source("utils.R")

######
# PARAMETER SETTING & DATA PREPROCESSING
######

# Set important parameters
bool_qc <- TRUE # TRUE if we want to generate the QC plots
nodox_qc_filter <- 50 # QC threshold - for each mutant, if its NO_DOX read is lower than this value, it will be replaced as NA to be filtered form downstream analyis

# Set metatable directory
meta_dir = "../../data/input/metatable_Dec2022.txt"

# Read the metatable that guides the column names to refer to
df_meta <- read.csv(file=meta_dir, header=TRUE, sep='\t')

# Get the column name of the plates
plate_col = unique(df_meta$plate_column)[1]

# Start the pipeline by pre-processing the input data to get LFC
df_all <- read.csv("../../data/input/in_vitro_raw_read_data.tsv", sep='\t')

column_groups_up <- list(
  FBS = c("A_FBS_UP_LFC",
          "B_FBS_UP_LFC",
          "C_FBS_UP_LFC",
          "D_FBS_UP_LFC",
          "E_FBS_UP_LFC",
          "F_FBS_UP_LFC"),
  Iron = c("A_FeFree_UP_LFC",
           "B_FeFree_UP_LFC",
           "C_FeFree_UP_LFC",
           "D_Iron_UP_LFC",
           "E_Iron_UP_LFC",
           "F_Iron_UP_LFC"),
  Temp37 = c("A_37_UP_LFC",
             "B_37_UP_LFC",
             "C_37_UP_LFC",
             "D_Temp37_UP_LFC",
             "E_Temp37_UP_LFC",
             "F_Temp37_UP_LFC"),
  YPD = c("A_YPD_UP_LFC",
          "B_YPD_UP_LFC",
          "C_YPD_UP_LFC",
          "D_YPD_UP_LFC",
          "E_YPD_UP_LFC",
          "F_YPD_UP_LFC"),
  NaCl = c("A_NaCl_UP_LFC",
           "B_NaCl_UP_LFC",
           "C_NaCl_UP_LFC",
           "NaCl_E_S2_L001_UP_LFC",
           "NaCl_D_S1_L001_UP_LFC",
           "NaCl_F_S3_L001_UP_LFC"),
  SDS = c("A_SDS_UP_LFC",
          "B_SDS_UP_LFC",
          "C_SDS_UP_LFC",
          "D_SDS_UP_LFC",
          "E_SDS_UP_LFC",
          "F_SDS_UP_LFC"),
  Sorbitol = c("A_Sorbitol_UP_LFC",
               "B_Sorbitol_UP_LFC",
               "C_Sorbitol_UP_LFC",
               "D_Sorbitol_UP_LFC",
               "E_Sorbitol_UP_LFC",
               "F_Sorbitol_UP_LFC"),
  YNB = c("A_CONT_UP_LFC",
          "B_CONT_UP_LFC",
          "C_CONT_UP_LFC",
          "A_YNB_BIO_UP_LFC",
          "B_YNB_BIO_UP_LFC",
          "C_YNB_BIO_UP_LFC",
          "D_YNB1_UP_LFC",
          "E_YNB1_UP_LFC",
          "F_YNB1_UP_LFC",
          "D_YNB2_UP_LFC",
          "E_YNB2_UP_LFC",
          "F_YNB2_UP_LFC",
          "D_YNB3_UP_LFC",
          "E_YNB3_UP_LFC",
          "F_YNB3_UP_LFC")
)

column_groups_dn <- list(
  FBS = c("A_FBS_DN_LFC",
          "B_FBS_DN_LFC",
          "C_FBS_DN_LFC",
          "D_FBS_DN_LFC",
          "E_FBS_DN_LFC",
          "F_FBS_DN_LFC"),
  Iron = c("A_FeFree_DN_LFC",
           "B_FeFree_DN_LFC",
           "C_FeFree_DN_LFC",
           "D_Iron_DN_LFC",
           "E_Iron_DN_LFC",
           "F_Iron_DN_LFC"),
  Temp37 = c("A_37_DN_LFC",
             "B_37_DN_LFC",
             "C_37_DN_LFC",
             "D_Temp37_DN_LFC",
             "E_Temp37_DN_LFC",
             "F_Temp37_DN_LFC"),
  YPD = c("A_YPD_DN_LFC",
          "B_YPD_DN_LFC",
          "C_YPD_DN_LFC",
          "D_YPD_DN_LFC",
          "E_YPD_DN_LFC",
          "F_YPD_DN_LFC"),
  NaCl = c("A_NaCl_DN_LFC",
           "B_NaCl_DN_LFC",
           "C_NaCl_DN_LFC",
           "NaCl_E_S2_L001_DOWN_LFC",
           "NaCl_D_S1_L001_DOWN_LFC",
           "NaCl_F_S3_L001_DOWN_LFC"),
  SDS = c("A_SDS_DN_LFC",
          "B_SDS_DN_LFC",
          "C_SDS_DN_LFC",
          "D_SDS_DN_LFC",
          "E_SDS_DN_LFC",
          "F_SDS_DN_LFC"),
  Sorbitol = c("A_Sorbitol_DN_LFC",
               "B_Sorbitol_DN_LFC",
               "C_Sorbitol_DN_LFC",
               "D_Sorbitol_DN_LFC",
               "E_Sorbitol_DN_LFC",
               "F_Sorbitol_DN_LFC"),
  YNB = c("A_CONT_DN_LFC",
          "B_CONT_DN_LFC",
          "C_CONT_DN_LFC",
          "A_YNB_BIO_DN_LFC",
          "B_YNB_BIO_DN_LFC",
          "C_YNB_BIO_DN_LFC",
          "D_YNB1_DN_LFC",
          "E_YNB1_DN_LFC",
          "F_YNB1_DN_LFC",
          "D_YNB2_DN_LFC",
          "E_YNB2_DN_LFC",
          "F_YNB2_DN_LFC",
          "D_YNB3_DN_LFC",
          "E_YNB3_DN_LFC",
          "F_YNB3_DN_LFC")
)

##############################################################
# Moderated t-test for single condition evaluation
##############################################################

# Conduct one sample moderated t-test on every condition, integrating UP and DN tags

# Set output and QC directory
output_dir = "../../data/output/output_scores/"
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
qc_output_dir = "../../data/output/qc/"
if (!dir.exists(qc_output_dir)) { dir.create(qc_output_dir, recursive = TRUE) }

# Create dataframe to summarize all results for global comparison
df_summary <- data.frame(plate = df_all$plate, row.names = df_all$orf19)

# Initialize lists to collect LFC dataframes for conditions (only needed for QC)
up_df_lfc_list <- c()
dn_df_lfc_list <- c()

# Iterate the scoring process through conditions
groups <- unique(df_meta$condition)
for (group in groups) {
  print(paste("Condition:",group))
  group_dir = paste0(output_dir, group, '/')
  if (!dir.exists(group_dir)) { dir.create(group_dir, recursive = TRUE) }
  
  df_meta_group <- df_meta[df_meta$condition==group, ]
  df_meta_group_dox <- df_meta_group[df_meta_group$DOX==TRUE, ]
  df_meta_group_nodox <- df_meta_group[df_meta_group$DOX==FALSE, ]
  conditions_dox_list <- c(df_meta_group_dox$UP, df_meta_group_dox$DN)
  conditions_no_dox_list <- c(df_meta_group_nodox$UP, df_meta_group_nodox$DN)
  
  # For each group (condition), process LFC values and score them
  df_lfc <- data_preprocess_get_lfc(df_all, conditions_dox_list, conditions_no_dox_list, no_dox_cutoff=nodox_qc_filter)
  
  # Plot saturation plot for baseline condition (YNB)
  if (group=="YNB") {
    replicates_UP <- group_list(column_groups_up$YNB)
    replicates_DN <- group_list(column_groups_dn$YNB)
    
    sorted_hits_rep_len_up <- saturation_plot_baseline_replicates(df_lfc, replicates_UP,
                                        output_direc=paste0(qc_output_dir, "YNB_UP_LFC_replicates_saturation_plot.pdf"),
                                        gene_col="orf19", lfc_cutoff=2, baseline_cond=group)
    sorted_hits_rep_len_dn <- saturation_plot_baseline_replicates(df_lfc, replicates_DN,
                                        output_direc=paste0(qc_output_dir, "YNB_DN_LFC_replicates_saturation_plot.pdf"),
                                        gene_col="orf19", lfc_cutoff=2, baseline_cond=group)
    iteration <- c(1:length(sorted_hits_rep_len_up))
    pdf(paste0(qc_output_dir, "YNB_LFC_replicates_saturation_plot.pdf"))  # Save the plot into a PDF file
    gap.plot(iteration, sorted_hits_rep_len_up, type = "b", pch = 1, 
             gap = c(5, 210), brw = 0.02,
             ytics=c(0, 220, 230, 240, 250), xtics=seq(0,max(iteration),by=1),
         #ylim = range(c(sorted_hits_rep_len_up, sorted_hits_rep_len_dn)),
         xlim = c(0, max(iteration)),
         ylim = c(0, max(c(sorted_hits_rep_len_up, sorted_hits_rep_len_dn))),
         xlab = paste0(group, " replicates"), ylab = "Number of total hits (LFC >= 2)",
         main = "Number of total hits given new replicates (add one at a time)")
    #lines(iteration, sorted_hits_rep_len_dn, type = "b", pch = 2)
    gap.plot(iteration, sorted_hits_rep_len_dn, type = "b", pch = 2, 
             gap = c(5, 210), brw = 0.02, add=TRUE,
             ytics=c(0, 220, 230, 240, 250), xtics=seq(0,max(iteration),by=1),
             #ylim = range(c(sorted_hits_rep_len_up, sorted_hits_rep_len_dn)),
             xlim = c(0, max(iteration)),
             ylim = c(0, max(c(sorted_hits_rep_len_up, sorted_hits_rep_len_dn))))
    axis.break(2, 5, breakcol="snow", style="gap")
    axis.break(2, 5.5, breakcol="black", style="slash")
    axis.break(4, 5.5, breakcol="black", style="slash")
    legend("topleft", legend = c("UP tag", "DN tag"), pch=c(1,2))
    dev.off() # Close the PDF device
  }
  
  # Generate QC plots (PCC heatmap & scatter plots)
  if (isTRUE(bool_qc)) {
    df_lfc_updn <- generate_qc_plots(df_lfc, lfc_up_list=paste0(df_meta_group_nodox$UP, "_LFC"),
                      lfc_dn_list=paste0(df_meta_group_nodox$DN, "_LFC"),
                      condition_name=group, output_direc=qc_output_dir, input_type="LFC", bk_min=0.5)
    up_df_lfc_list <- c(up_df_lfc_list, df_lfc_updn$lfc_up)
    dn_df_lfc_list <- c(dn_df_lfc_list, df_lfc_updn$lfc_dn)
  }

}

# Plot global LFC PCC heatmap
if (isTRUE(bool_qc)) {
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

##############################################################
# Moderated t-test for condition interaction effect compared to reference standard
##############################################################

# Conduct one sample moderated t-test on the difference between every condition and reference, integrating UP and DN tags

# Define output directory:
dLFC_output_dir <- "../../data/output/output_each_vs_REF/"
if (!dir.exists(dLFC_output_dir)) { dir.create(dLFC_output_dir, recursive = TRUE) }

# Create dataframe to summarize all results for global comparison
df_summary_dLFC <- data.frame(plate = df_all$plate, row.names = df_all$orf19)

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
  group_dir = paste0(dLFC_output_dir, group, '/')
  if (!dir.exists(group_dir)) { dir.create(group_dir, recursive = TRUE) }

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

# Plot global LFC PCC heatmap
if (isTRUE(bool_qc)) {
  df_dLFC_cond_up_all <- do.call(cbind, up_df_dlfc_list)
  df_dLFC_cond_dn_all <- do.call(cbind, dn_df_dlfc_list)
  plot_qc_correlation(df_dLFC_cond_up_all, paste0(qc_output_dir, "global_heatmap_UP"),
                      show_num=FALSE, font_size=6, bk_min=0.0, color_interval=100, input_type="dLFC")
  plot_qc_correlation(df_dLFC_cond_dn_all, paste0(qc_output_dir, "global_heatmap_DN"),
                      show_num=FALSE, font_size=6, bk_min=0.0, color_interval=100, input_type="dLFC")
}

