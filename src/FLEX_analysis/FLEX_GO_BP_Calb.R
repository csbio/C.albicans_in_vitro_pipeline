## ============================ ##
## profile similarity co-annotation enrichment analysis using the FLEX R Package
## ============================ ##

setwd('C.albicans_in_vitro_pipeline/src/FLEX_analysis')

require(devtools)

# Link to the FLEX R package: https://github.com/csbio/FLEX_R
load_all('replace_with_FLEX_local_install_directory/FLEX_R-master/R')

generate_stepwise_contribution_data<-function(complex_info, data_complex_info, output_folder, output_file_name_noext)
{
  pairs_in <- data.frame(true = complex_info$true, 
                         predicted = complex_info$predicted, 
                         ID =  as.character(complex_info$ID), 
                         stringsAsFactors = FALSE)
  stepwise_contrib <- GenerateDataForPerfCurve(value.predicted = complex_info$predicted, 
                                               value.true = complex_info$true, 
                                               x.axis = 'TP', y.axis = 'precision')
  
  precision_cutoffs <- c(stepwise_contrib$y[length(stepwise_contrib$y)], seq(0.01, max(stepwise_contrib$y), 0.01)) #.1 0.025
  precision_cutoffs[1] <- round(precision_cutoffs[1], 4)
  
  stepwise_contrib <- GetStepwiseContributionOfEntities(pairs_in, precision_cutoffs, data_complex_info)
  stepwise_contrib_file <- paste('Stepwise_cont_',output_file_name_noext,'.tsv', sep = "")
  write.table(stepwise_contrib, stepwise_contrib_file, sep="\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  print(stepwise_contrib_file)
  
  #cols<-brewer.pal(10, "Paired")
  stepwise_contrib <- stepwise_contrib[!duplicated(stepwise_contrib$Name), ]
  structure_plot_file <- paste(output_folder, 'Cont_struct_', output_file_name_noext, sep = "")
  PlotContributionStructure(stepwise_contrib, cutoff.all = precision_cutoffs, min.precision.cutoff = 0.01 , # 0.1
                            min.pairs = 10, outfile.name = structure_plot_file, outfile.type = 'pdf',
                            save.figure = TRUE)
  files <- c( 
    stepwise_contrib_file, paste(structure_plot_file, '.pdf', sep = "")
  )
  # for (f in files)
  # {
  #   file.copy(f, file.path(output_folder, f), overwrite = TRUE)
  #   file.remove(f)
  # }
}

generate_auprc_data<-function(complex_info, data_complex_info, corum_info, output_folder, output_file_name_noext)
{
  complex_info$ID=as.character(complex_info$ID)
  complex_df <- as.data.frame(complex_info, stringsAsFactors = FALSE)
  data_AUPRC <- GetAreaUnderPRCurveForEntities (data_complex_info, corum_info, complex_df)
  auprc_file <- paste('AUPRC_',output_file_name_noext,'.txt', sep = "")
  write.table(data_AUPRC, auprc_file, sep = '\t', row.names = FALSE, quote = FALSE)
  files <- c( 
    auprc_file    
  )
  for (f in files)
  {
    file.copy(f, file.path(output_folder, f), overwrite = TRUE)
    file.remove(f)
  }
}


# 1. Read in the data included with the package
# ---------------------------------------------
# GO BP terms with less than 5 or more than 200 genes annotated were filtered for more specific results
data_GO_BP = read.table("../../data/input/go_bp_standard_Calb_5_200terms.txt", sep='\t', header = TRUE,quote = "")

# 2. Create (or read once created) the co-annotation data
# -------------------------------------------------------
file_name <- "GO_BP_Calb_5_200.Rdata"
data.ca <- MakeCoAnnotationFromGeneSymbols(data_standard = data_GO_BP, 
                                           overlap_length = 1, 
                                           file_location = file_name)


# 3. Read in the interaction/dependency data
# ---------------------------------------------
file.int <- '../../data/output/summary_mod_t_test_compare_median_normvar_LFCgeq1_noNA_noYPD.txt'

# tab delimited
data.interaction <- GetInteractionData(file.int)
data.interaction <- na.omit(data.interaction)
# LFC>=1 cutoff
row_max <- apply(data.interaction, 1, max)
cutoff <- 1 # set your maximum cutoff value
data.interaction <- data.interaction[row_max >= cutoff, ]


# 4. Associate the pairwise scores to co-annotation
# ---------------------------------------------
GOBP_Calb_set1 <- CalculatePredictionAndTrueOnLibraryProfiles(data.ca, data.interaction)
save(GOBP_Calb_set1, file = '../../data/output/PR_curves/GOBP_Calb_set1_5_200_noNA_LFCgeq1_noYPD.Rdata')


# 5. Plot global PR Curves
# ---------------------------------------------
load('../../data/output/PR_curves/GOBP_Calb_set1_5_200_noNA_LFCgeq1_noYPD.Rdata')

pred.ca <- list(out_Calb_1 = list(true = GOBP_Calb_set1$true, 
                                  predicted = GOBP_Calb_set1$predicted))

PlotPRSimilarity (pred.ca, subsample = TRUE, type.plot = 'log',
                  fig.title = 'C.albicans GO_BP co-annotation', 
                  fig.labs = c('TP', 'Precision'), legend.names = c('C.albicans in vitro set (LFC>=1, noNA, noYPD)'), 
                  legend.color = c('#de2d26'), save.figure = TRUE, #, '#3182bd'
                  outfile.name = '../../data/output/PR_curves/Calbicans_GOBP_5_200_noNA_LFCgeq1_noYPD', outfile.type = 'pdf')


# 6. Individual AUPRC (Contribution Scatter)
# ---------------------------------------------
# entity.matrix <- as.data.frame(GOBP_Calb_set1, stringsAsFactors = FALSE)
# data.AUPRC <- GetAreaUnderPRCurveForEntities (summary.standard = data_GO_BP, data.standard = data.ca, entity.matrix = entity.matrix)
# write.table(data.AUPRC, 'GO_BP_AUPRC_Calb_set1.txt', sep = '\t', row.names = FALSE, quote = FALSE)
# 
# #data.AUPRC <- read.table('Pathway_AUPRC_DepMap_19Q2_set1.txt', stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
# plot.data <- data.AUPRC[!duplicated(data.AUPRC$Name), ]
# 
# # Don't advise using show.text = TRUE, looks like a mess
# PlotContributionScatter (plot.data, length.cutoff = 10, AUPRC.cutoff = 0.1, fig.labs = c('AUPRC', 'Pathway size'), show.text = FALSE, save.figure = TRUE, outfile.name = 'Pathway_Scatter_Plot')

# 7. Diversity plot
output_folder <- file.path("../../data/output/diversity_plot/")
# auprc_file_name_noext <- paste("test_output",sep = '')
# generate_auprc_data(ss, data_GO_BP, corum, output_folder, auprc_file_name_noext)
stepwise_contribution_file_name_noext <- paste("GO_BP_5_200_noNA_LFCgeq1_noYPD",sep = '')	
generate_stepwise_contribution_data(GOBP_Calb_set1, data_GO_BP, output_folder, stepwise_contribution_file_name_noext)
