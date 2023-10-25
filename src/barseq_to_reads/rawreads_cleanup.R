#-------------------------------------------------
# Raw sequencing reads data clean-up
#-------------------------------------------------
# Last updated October 24th, 2023 by Emily Xion

library(gplots) #for heatmap
library(ggplot2) # for plotting
library(ggrepel) # for text spaces in ggplot if there are too many texts
library(reshape2) # for table management
library(GGally) # for correlation plots
library(dplyr) # for table management
library(tidyr)
library(ggvenn)

setwd("/Volumes/EXIONG_T7/Cowen_BARSEQ/ALL") #set working directory to file with data 
#load data
og_data <- read.csv("merged_results.csv", header = TRUE)
orf19 <- og_data$orf19
feature <- og_data$Feature
common <- og_data$Common
description <- og_data$Description
plate <- og_data$GRACE.Plate
position <- og_data$GRACE.Position

#making regex for sequencing technical replicates
regex_FBSAUP <- "(FBS_A_S\\d+_L00\\d+_UP)"
regex_FBSADN <- "(FBS_A_S\\d+_L00\\d+_DOWN)"
regex_FBSBUP <- "(FBS_B_S\\d+_L00\\d+_UP)"
regex_FBSBDN <- "(FBS_B_S\\d+_L00\\d+_DOWN)"
regex_FBSCUP <- "(FBS_C_S\\d+_L00\\d+_UP)"
regex_FBSCDN <- "(FBS_C_S\\d+_L00\\d+_DOWN)"
#--
regex_FBSDOXAUP <- "(FBS_DOX_A_S\\d+_L00\\d+_UP)"
regex_FBSDOXADN <- "(FBS_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_FBSDOXBUP <- "(FBS_DOX_B_S\\d+_L00\\d+_UP)"
regex_FBSDOXBDN <- "(FBS_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_FBSDOXCUP <- "(FBS_DOX_C_S\\d+_L00\\d+_UP)"
regex_FBSDOXCDN <- "(FBS_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_37AUP <- "(37_A_S\\d+_L00\\d+_UP)"
regex_37ADN <- "(37_A_S\\d+_L00\\d+_DOWN)"
regex_37BUP <- "(37_B_S\\d+_L00\\d+_UP)"
regex_37BDN <- "(37_B_S\\d+_L00\\d+_DOWN)"
regex_37CUP <- "(37_C_S\\d+_L00\\d+_UP)"
regex_37CDN <- "(37_C_S\\d+_L00\\d+_DOWN)"
#--
regex_37DOXAUP <- "(37_DOX_A_S\\d+_L00\\d+_UP)"
regex_37DOXADN <- "(37_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_37DOXBUP <- "(37_DOX_B_S\\d+_L00\\d+_UP)"
regex_37DOXBDN <- "(37_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_37DOXCUP <- "(37_DOX_C_S\\d+_L00\\d+_UP)"
regex_37DOXCDN <- "(37_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_FeFreeAUP <- "(FeFree_A_S\\d+_L00\\d+_UP)"
regex_FeFreeADN <- "(FeFree_A_S\\d+_L00\\d+_DOWN)"
regex_FeFreeBUP <- "(FeFree_B_S\\d+_L00\\d+_UP)"
regex_FeFreeBDN <- "(FeFree_B_S\\d+_L00\\d+_DOWN)"
regex_FeFreeCUP <- "(FeFree_C_S\\d+_L00\\d+_UP)"
regex_FeFreeCDN <- "(FeFree_C_S\\d+_L00\\d+_DOWN)"
#--
regex_FeFreeDOXAUP <- "(FeFree_DOX_A_S\\d+_L00\\d+_UP)"
regex_FeFreeDOXADN <- "(FeFree_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_FeFreeDOXBUP <- "(FeFree_DOX_B_S\\d+_L00\\d+_UP)"
regex_FeFreeDOXBDN <- "(FeFree_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_FeFreeDOXCUP <- "(FeFree_DOX_C_S\\d+_L00\\d+_UP)"
regex_FeFreeDOXCDN <- "(FeFree_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_YPDAUP <- "(YPD_A_S\\d+_L00\\d+_UP)"
regex_YPDADN <- "(YPD_A_S\\d+_L00\\d+_DOWN)"
regex_YPDBUP <- "(YPD_B_S\\d+_L00\\d+_UP)"
regex_YPDBDN <- "(YPD_B_S\\d+_L00\\d+_DOWN)"
regex_YPDCUP <- "(YPD_C_S\\d+_L00\\d+_UP)"
regex_YPDCDN <- "(YPD_C_S\\d+_L00\\d+_DOWN)"
#--
regex_YPDDOXAUP <- "(YPD_DOX_A_S\\d+_L00\\d+_UP)"
regex_YPDDOXADN <- "(YPD_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_YPDDOXBUP <- "(YPD_DOX_B_S\\d+_L00\\d+_UP)"
regex_YPDDOXBDN <- "(YPD_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_YPDDOXCUP <- "(YPD_DOX_C_S\\d+_L00\\d+_UP)"
regex_YPDDOXCDN <- "(YPD_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_CONTAUP <- "(CONT_A_S\\d+_L00\\d+_UP)"
regex_CONTADN <- "(CONT_A_S\\d+_L00\\d+_DOWN)"
regex_CONTBUP <- "(CONT_B_S\\d+_L00\\d+_UP)"
regex_CONTBDN <- "(CONT_B_S\\d+_L00\\d+_DOWN)"
regex_CONTCUP <- "(CONT_C_S\\d+_L00\\d+_UP)"
regex_CONTCDN <- "(CONT_C_S\\d+_L00\\d+_DOWN)"
#--
regex_CONTDOXAUP <- "(CONT_DOX_A_S\\d+_L00\\d+_UP)"
regex_CONTDOXADN <- "(CONT_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_CONTDOXBUP <- "(CONT_DOX_B_S\\d+_L00\\d+_UP)"
regex_CONTDOXBDN <- "(CONT_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_CONTDOXCUP <- "(CONT_DOX_C_S\\d+_L00\\d+_UP)"
regex_CONTDOXCDN <- "(CONT_DOX_C_S\\d+_L00\\d+_DOWN)"

#--------------

regex_NaClAUP <- "(NaCl_A_S\\d+_L00\\d+_UP)"
regex_NaClADN <- "(NaCl_A_S\\d+_L00\\d+_DOWN)"
regex_NaClBUP <- "(NaCl_B_S\\d+_L00\\d+_UP)"
regex_NaClBDN <- "(NaCl_B_S\\d+_L00\\d+_DOWN)"
regex_NaClCUP <- "(NaCl_C_S\\d+_L00\\d+_UP)"
regex_NaClCDN <- "(NaCl_C_S\\d+_L00\\d+_DOWN)"
#--
regex_NaClDOXAUP <- "(NaCl_DOX_A_S\\d+_L00\\d+_UP)"
regex_NaClDOXADN <- "(NaCl_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_NaClDOXBUP <- "(NaCl_DOX_B_S\\d+_L00\\d+_UP)"
regex_NaClDOXBDN <- "(NaCl_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_NaClDOXCUP <- "(NaCl_DOX_C_S\\d+_L00\\d+_UP)"
regex_NaClDOXCDN <- "(NaCl_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_SDSAUP <- "(SDS_A_S\\d+_L00\\d+_UP)"
regex_SDSADN <- "(SDS_A_S\\d+_L00\\d+_DOWN)"
regex_SDSBUP <- "(SDS_B_S\\d+_L00\\d+_UP)"
regex_SDSBDN <- "(SDS_B_S\\d+_L00\\d+_DOWN)"
regex_SDSCUP <- "(SDS_C_S\\d+_L00\\d+_UP)"
regex_SDSCDN <- "(SDS_C_S\\d+_L00\\d+_DOWN)"
#--
regex_SDSDOXAUP <- "(SDS_DOX_A_S\\d+_L00\\d+_UP)"
regex_SDSDOXADN <- "(SDS_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_SDSDOXBUP <- "(SDS_DOX_B_S\\d+_L00\\d+_UP)"
regex_SDSDOXBDN <- "(SDS_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_SDSDOXCUP <- "(SDS_DOX_C_S\\d+_L00\\d+_UP)"
regex_SDSDOXCDN <- "(SDS_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_SorbitolAUP <- "(Sorbitol_A_S\\d+_L00\\d+_UP)"
regex_SorbitolADN <- "(Sorbitol_A_S\\d+_L00\\d+_DOWN)"
regex_SorbitolBUP <- "(Sorbitol_B_S\\d+_L00\\d+_UP)"
regex_SorbitolBDN <- "(Sorbitol_B_S\\d+_L00\\d+_DOWN)"
regex_SorbitolCUP <- "(Sorbitol_C_S\\d+_L00\\d+_UP)"
regex_SorbitolCDN <- "(Sorbitol_C_S\\d+_L00\\d+_DOWN)"
#--
regex_SorbitolDOXAUP <- "(Sorbitol_DOX_A_S\\d+_L00\\d+_UP)"
regex_SorbitolDOXADN <- "(Sorbitol_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_SorbitolDOXBUP <- "(Sorbitol_DOX_B_S\\d+_L00\\d+_UP)"
regex_SorbitolDOXBDN <- "(Sorbitol_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_SorbitolDOXCUP <- "(Sorbitol_DOX_C_S\\d+_L00\\d+_UP)"
regex_SorbitolDOXCDN <- "(Sorbitol_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_YNBAUP <- "(YNB_A_S\\d+_L00\\d+_UP)"
regex_YNBADN <- "(YNB_A_S\\d+_L00\\d+_DOWN)"
regex_YNBBUP <- "(YNB_B_S\\d+_L00\\d+_UP)"
regex_YNBBDN <- "(YNB_B_S\\d+_L00\\d+_DOWN)"
regex_YNBCUP <- "(YNB_C_S\\d+_L00\\d+_UP)"
regex_YNBCDN <- "(YNB_C_S\\d+_L00\\d+_DOWN)"
#--
regex_YNBDOXAUP <- "(YNB_DOX_A_S\\d+_L00\\d+_UP)"
regex_YNBDOXADN <- "(YNB_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_YNBDOXBUP <- "(YNB_DOX_B_S\\d+_L00\\d+_UP)"
regex_YNBDOXBDN <- "(YNB_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_YNBDOXCUP <- "(YNB_DOX_C_S\\d+_L00\\d+_UP)"
regex_YNBDOXCDN <- "(YNB_DOX_C_S\\d+_L00\\d+_DOWN)"

#average sequencing technical replicates
A_FBS_UP<- rowSums(og_data[,grep(regex_FBSAUP, colnames(og_data))])
A_FBS_DN <- rowSums(og_data[,grep(regex_FBSADN, colnames(og_data))])
B_FBS_UP <- rowSums(og_data[,grep(regex_FBSBUP, colnames(og_data))])
B_FBS_DN <- rowSums(og_data[,grep(regex_FBSBDN, colnames(og_data))])
C_FBS_UP <- rowSums(og_data[,grep(regex_FBSCUP, colnames(og_data))])
C_FBS_DN <- rowSums(og_data[,grep(regex_FBSCDN, colnames(og_data))])
#--
A_FBS_DOX_UP <- rowSums(og_data[,grep(regex_FBSDOXAUP, colnames(og_data))])
A_FBS_DOX_DN <- rowSums(og_data[,grep(regex_FBSDOXADN, colnames(og_data))])
B_FBS_DOX_UP <- rowSums(og_data[,grep(regex_FBSDOXBUP, colnames(og_data))])
B_FBS_DOX_DN <- rowSums(og_data[,grep(regex_FBSDOXBDN, colnames(og_data))])
C_FBS_DOX_UP <- rowSums(og_data[,grep(regex_FBSDOXCUP, colnames(og_data))])
C_FBS_DOX_DN <- rowSums(og_data[,grep(regex_FBSDOXCDN, colnames(og_data))])
#--
A_37_UP <- rowSums(og_data[,grep(regex_37AUP, colnames(og_data))])
A_37_DN <- rowSums(og_data[,grep(regex_37ADN, colnames(og_data))])
B_37_UP <- rowSums(og_data[,grep(regex_37BUP, colnames(og_data))])
B_37_DN <- rowSums(og_data[,grep(regex_37BDN, colnames(og_data))])
C_37_UP <- rowSums(og_data[,grep(regex_37CUP, colnames(og_data))])
C_37_DN <- rowSums(og_data[,grep(regex_37CDN, colnames(og_data))])
#--
A_37_DOX_UP <- rowSums(og_data[,grep(regex_37DOXAUP, colnames(og_data))])
A_37_DOX_DN <- rowSums(og_data[,grep(regex_37DOXADN, colnames(og_data))])
B_37_DOX_UP <- rowSums(og_data[,grep(regex_37DOXBUP, colnames(og_data))])
B_37_DOX_DN <- rowSums(og_data[,grep(regex_37DOXBDN, colnames(og_data))])
C_37_DOX_UP <- rowSums(og_data[,grep(regex_37DOXCUP, colnames(og_data))])
C_37_DOX_DN <- rowSums(og_data[,grep(regex_37DOXCDN, colnames(og_data))])
#--
A_FeFree_UP <- rowSums(og_data[,grep(regex_FeFreeAUP, colnames(og_data))])
A_FeFree_DN <- rowSums(og_data[,grep(regex_FeFreeADN, colnames(og_data))])
B_FeFree_UP <- rowSums(og_data[,grep(regex_FeFreeBUP, colnames(og_data))])
B_FeFree_DN <- rowSums(og_data[,grep(regex_FeFreeBDN, colnames(og_data))])
C_FeFree_UP <- rowSums(og_data[,grep(regex_FeFreeCUP, colnames(og_data))])
C_FeFree_DN <- rowSums(og_data[,grep(regex_FeFreeCDN, colnames(og_data))])
#--
A_FeFree_DOX_UP <- rowSums(og_data[,grep(regex_FeFreeDOXAUP, colnames(og_data))])
A_FeFree_DOX_DN <- rowSums(og_data[,grep(regex_FeFreeDOXADN, colnames(og_data))])
B_FeFree_DOX_UP <- rowSums(og_data[,grep(regex_FeFreeDOXBUP, colnames(og_data))])
B_FeFree_DOX_DN <- rowSums(og_data[,grep(regex_FeFreeDOXBDN, colnames(og_data))])
C_FeFree_DOX_UP <- rowSums(og_data[,grep(regex_FeFreeDOXCUP, colnames(og_data))])
C_FeFree_DOX_DN <- rowSums(og_data[,grep(regex_FeFreeDOXCDN, colnames(og_data))])
#--
A_YPD_UP <- rowSums(og_data[,grep(regex_YPDAUP, colnames(og_data))])
A_YPD_DN <- rowSums(og_data[,grep(regex_YPDADN, colnames(og_data))])
B_YPD_UP <- rowSums(og_data[,grep(regex_YPDBUP, colnames(og_data))])
B_YPD_DN <- rowSums(og_data[,grep(regex_YPDBDN, colnames(og_data))])
C_YPD_UP <- rowSums(og_data[,grep(regex_YPDCUP, colnames(og_data))])
C_YPD_DN <- rowSums(og_data[,grep(regex_YPDCDN, colnames(og_data))])
#--
A_YPD_DOX_UP <- rowSums(og_data[,grep(regex_YPDDOXAUP, colnames(og_data))])
A_YPD_DOX_DN <- rowSums(og_data[,grep(regex_YPDDOXADN, colnames(og_data))])
B_YPD_DOX_UP <- rowSums(og_data[,grep(regex_YPDDOXBUP, colnames(og_data))])
B_YPD_DOX_DN <- rowSums(og_data[,grep(regex_YPDDOXBDN, colnames(og_data))])
C_YPD_DOX_UP <- rowSums(og_data[,grep(regex_YPDDOXCUP, colnames(og_data))])
C_YPD_DOX_DN <- rowSums(og_data[,grep(regex_YPDDOXCDN, colnames(og_data))])
#--
A_CONT_UP <- rowSums(og_data[,grep(regex_CONTAUP, colnames(og_data))])
A_CONT_DN <- rowSums(og_data[,grep(regex_CONTADN, colnames(og_data))])
B_CONT_UP <- rowSums(og_data[,grep(regex_CONTBUP, colnames(og_data))])
B_CONT_DN <- rowSums(og_data[,grep(regex_CONTBDN, colnames(og_data))])
C_CONT_UP <- rowSums(og_data[,grep(regex_CONTCUP, colnames(og_data))])
C_CONT_DN <- rowSums(og_data[,grep(regex_CONTCDN, colnames(og_data))])
#--
A_CONT_DOX_UP <- rowSums(og_data[,grep(regex_CONTDOXAUP, colnames(og_data))])
A_CONT_DOX_DN <- rowSums(og_data[,grep(regex_CONTDOXADN, colnames(og_data))])
B_CONT_DOX_UP <- rowSums(og_data[,grep(regex_CONTDOXBUP, colnames(og_data))])
B_CONT_DOX_DN <- rowSums(og_data[,grep(regex_CONTDOXBDN, colnames(og_data))])
C_CONT_DOX_UP <- rowSums(og_data[,grep(regex_CONTDOXCUP, colnames(og_data))])
C_CONT_DOX_DN <- rowSums(og_data[,grep(regex_CONTDOXCDN, colnames(og_data))])
#--
A_NaCl_UP <- rowSums(og_data[,grep(regex_NaClAUP, colnames(og_data))])
A_NaCl_DN <- rowSums(og_data[,grep(regex_NaClADN, colnames(og_data))])
B_NaCl_UP <- rowSums(og_data[,grep(regex_NaClBUP, colnames(og_data))])
B_NaCl_DN <- rowSums(og_data[,grep(regex_NaClBDN, colnames(og_data))])
C_NaCl_UP <- rowSums(og_data[,grep(regex_NaClCUP, colnames(og_data))])
C_NaCl_DN <- rowSums(og_data[,grep(regex_NaClCDN, colnames(og_data))])
#--
A_NaCl_DOX_UP <- rowSums(og_data[,grep(regex_NaClDOXAUP, colnames(og_data))])
A_NaCl_DOX_DN <- rowSums(og_data[,grep(regex_NaClDOXADN, colnames(og_data))])
B_NaCl_DOX_UP <- rowSums(og_data[,grep(regex_NaClDOXBUP, colnames(og_data))])
B_NaCl_DOX_DN <- rowSums(og_data[,grep(regex_NaClDOXBDN, colnames(og_data))])
C_NaCl_DOX_UP <- rowSums(og_data[,grep(regex_NaClDOXCUP, colnames(og_data))])
C_NaCl_DOX_DN <- rowSums(og_data[,grep(regex_NaClDOXCDN, colnames(og_data))])
#--
A_SDS_UP <- rowSums(og_data[,grep(regex_SDSAUP, colnames(og_data))])
A_SDS_DN <- rowSums(og_data[,grep(regex_SDSADN, colnames(og_data))])
B_SDS_UP <- rowSums(og_data[,grep(regex_SDSBUP, colnames(og_data))])
B_SDS_DN <- rowSums(og_data[,grep(regex_SDSBDN, colnames(og_data))])
C_SDS_UP <- rowSums(og_data[,grep(regex_SDSCUP, colnames(og_data))])
C_SDS_DN <- rowSums(og_data[,grep(regex_SDSCDN, colnames(og_data))])
#--
A_SDS_DOX_UP <- rowSums(og_data[,grep(regex_SDSDOXAUP, colnames(og_data))])
A_SDS_DOX_DN <- rowSums(og_data[,grep(regex_SDSDOXADN, colnames(og_data))])
B_SDS_DOX_UP <- rowSums(og_data[,grep(regex_SDSDOXBUP, colnames(og_data))])
B_SDS_DOX_DN <- rowSums(og_data[,grep(regex_SDSDOXBDN, colnames(og_data))])
C_SDS_DOX_UP <- rowSums(og_data[,grep(regex_SDSDOXCUP, colnames(og_data))])
C_SDS_DOX_DN <- rowSums(og_data[,grep(regex_SDSDOXCDN, colnames(og_data))])
#--
A_Sorbitol_UP <- rowSums(og_data[,grep(regex_SorbitolAUP, colnames(og_data))])
A_Sorbitol_DN <- rowSums(og_data[,grep(regex_SorbitolADN, colnames(og_data))])
B_Sorbitol_UP <- rowSums(og_data[,grep(regex_SorbitolBUP, colnames(og_data))])
B_Sorbitol_DN <- rowSums(og_data[,grep(regex_SorbitolBDN, colnames(og_data))])
C_Sorbitol_UP <- rowSums(og_data[,grep(regex_SorbitolCUP, colnames(og_data))])
C_Sorbitol_DN <- rowSums(og_data[,grep(regex_SorbitolCDN, colnames(og_data))])
#--
A_Sorbitol_DOX_UP <- rowSums(og_data[,grep(regex_SorbitolDOXAUP, colnames(og_data))])
A_Sorbitol_DOX_DN <- rowSums(og_data[,grep(regex_SorbitolDOXADN, colnames(og_data))])
B_Sorbitol_DOX_UP <- rowSums(og_data[,grep(regex_SorbitolDOXBUP, colnames(og_data))])
B_Sorbitol_DOX_DN <- rowSums(og_data[,grep(regex_SorbitolDOXBDN, colnames(og_data))])
C_Sorbitol_DOX_UP <- rowSums(og_data[,grep(regex_SorbitolDOXCUP, colnames(og_data))])
C_Sorbitol_DOX_DN <- rowSums(og_data[,grep(regex_SorbitolDOXCDN, colnames(og_data))])
#--
A_YNB_BIO_UP <- rowSums(og_data[,grep(regex_YNBAUP, colnames(og_data))])
A_YNB_BIO_DN <- rowSums(og_data[,grep(regex_YNBADN, colnames(og_data))])
B_YNB_BIO_UP <- rowSums(og_data[,grep(regex_YNBBUP, colnames(og_data))])
B_YNB_BIO_DN <- rowSums(og_data[,grep(regex_YNBBDN, colnames(og_data))])
C_YNB_BIO_UP <- rowSums(og_data[,grep(regex_YNBCUP, colnames(og_data))])
C_YNB_BIO_DN <- rowSums(og_data[,grep(regex_YNBCDN, colnames(og_data))])
#--
A_YNB_BIO_DOX_UP <- rowSums(og_data[,grep(regex_YNBDOXAUP, colnames(og_data))])
A_YNB_BIO_DOX_DN <- rowSums(og_data[,grep(regex_YNBDOXADN, colnames(og_data))])
B_YNB_BIO_DOX_UP <- rowSums(og_data[,grep(regex_YNBDOXBUP, colnames(og_data))])
B_YNB_BIO_DOX_DN <- rowSums(og_data[,grep(regex_YNBDOXBDN, colnames(og_data))])
C_YNB_BIO_DOX_UP <- rowSums(og_data[,grep(regex_YNBDOXCUP, colnames(og_data))])
C_YNB_BIO_DOX_DN <- rowSums(og_data[,grep(regex_YNBDOXCDN, colnames(og_data))])

#generate dataframe with all summed barcode counts 
sum_data <- cbind(orf19, feature, common, description, plate, position,
                  A_FBS_UP, A_FBS_DN,
                  B_FBS_UP, B_FBS_DN,
                  C_FBS_UP, C_FBS_DN,
                  A_FeFree_UP, A_FeFree_DN,
                  B_FeFree_UP, B_FeFree_DN,
                  C_FeFree_UP, C_FeFree_DN,
                  A_37_UP, A_37_DN,
                  B_37_UP, B_37_DN,
                  C_37_UP, C_37_DN,
                  A_YPD_UP, A_YPD_DN,
                  B_YPD_UP, B_YPD_DN,
                  C_YPD_UP, C_YPD_DN,
                  A_CONT_UP, A_CONT_DN,
                  B_CONT_UP, B_CONT_DN,
                  C_CONT_UP, C_CONT_DN,
                  A_NaCl_UP, A_NaCl_DN,
                  B_NaCl_UP, B_NaCl_DN,
                  C_NaCl_UP, C_NaCl_DN,
                  A_SDS_UP, A_SDS_DN,
                  B_SDS_UP, B_SDS_DN,
                  C_SDS_UP, C_SDS_DN,
                  A_Sorbitol_UP, A_Sorbitol_DN,
                  B_Sorbitol_UP, B_Sorbitol_DN,
                  C_Sorbitol_UP, C_Sorbitol_DN,
                  A_YNB_BIO_UP, A_YNB_BIO_DN,
                  B_YNB_BIO_UP, B_YNB_BIO_DN,
                  C_YNB_BIO_UP, C_YNB_BIO_DN,
                  A_FBS_DOX_UP, A_FBS_DOX_DN,
                  B_FBS_DOX_UP, B_FBS_DOX_DN,
                  C_FBS_DOX_UP, C_FBS_DOX_DN,
                  A_FeFree_DOX_UP, A_FeFree_DOX_DN,
                  B_FeFree_DOX_UP, B_FeFree_DOX_DN,
                  C_FeFree_DOX_UP, C_FeFree_DOX_DN,
                  A_37_DOX_UP, A_37_DOX_DN,
                  B_37_DOX_UP, B_37_DOX_DN,
                  C_37_DOX_UP, C_37_DOX_DN,
                  A_YPD_DOX_UP, A_YPD_DOX_DN,
                  B_YPD_DOX_UP, B_YPD_DOX_DN,
                  C_YPD_DOX_UP, C_YPD_DOX_DN,
                  A_CONT_DOX_UP, A_CONT_DOX_DN,
                  B_CONT_DOX_UP, B_CONT_DOX_DN,
                  C_CONT_DOX_UP, C_CONT_DOX_DN,
                  A_NaCl_DOX_UP, A_NaCl_DOX_DN,
                  B_NaCl_DOX_UP, B_NaCl_DOX_DN,
                  C_NaCl_DOX_UP, C_NaCl_DOX_DN,
                  A_SDS_DOX_UP, A_SDS_DOX_DN,
                  B_SDS_DOX_UP, B_SDS_DOX_DN,
                  C_SDS_DOX_UP, C_SDS_DOX_DN,
                  A_Sorbitol_DOX_UP, A_Sorbitol_DOX_DN,
                  B_Sorbitol_DOX_UP, B_Sorbitol_DOX_DN,
                  C_Sorbitol_DOX_UP, C_Sorbitol_DOX_DN,
                  A_YNB_BIO_DOX_UP, A_YNB_BIO_DOX_DN,
                  B_YNB_BIO_DOX_UP, B_YNB_BIO_DOX_DN,
                  C_YNB_BIO_DOX_UP, C_YNB_BIO_DOX_DN)
write.csv(sum_data, file = "sum_data.csv")

#data pruning to remove strains with very low/inconsistent read counts
sum_data <- read.csv("sum_data.csv", header = T)
v1 <- c("A_FBS_UP", "A_FBS_DN",
        "B_FBS_UP", "B_FBS_DN",
        "C_FBS_UP", "C_FBS_DN",
        "A_FeFree_UP", "A_FeFree_DN",
        "B_FeFree_UP", "B_FeFree_DN",
        "C_FeFree_UP", "C_FeFree_DN",
        "A_37_UP", "A_37_DN",
        "B_37_UP", "B_37_DN",
        "C_37_UP", "C_37_DN",
        "A_YPD_UP", "A_YPD_DN",
        "B_YPD_UP", "B_YPD_DN",
        "C_YPD_UP", "C_YPD_DN",
        "A_CONT_UP", "A_CONT_DN",
        "B_CONT_UP", "B_CONT_DN",
        "C_CONT_UP", "C_CONT_DN",
        "A_NaCl_UP", "A_NaCl_DN",
        "B_NaCl_UP", "B_NaCl_DN",
        "C_NaCl_UP", "C_NaCl_DN",
        "A_SDS_UP", "A_SDS_DN",
        "B_SDS_UP", "B_SDS_DN",
        "C_SDS_UP", "C_SDS_DN",
        "A_Sorbitol_UP", "A_Sorbitol_DN",
        "B_Sorbitol_UP", "B_Sorbitol_DN",
        "C_Sorbitol_UP", "C_Sorbitol_DN",
        "A_YNB_BIO_UP", "A_YNB_BIO_DN",
        "B_YNB_BIO_UP", "B_YNB_BIO_DN",
        "C_YNB_BIO_UP", "C_YNB_BIO_DN")
z <- sum_data
z[v1][z[v1] < 50] <- NA
sum_data_purgelow <- z

sum_data$sum <- rowSums(sum_data[ , c(8:115)], na.rm=TRUE)
sum_data_purgelow <- sum_data[sum_data$sum > 150, ]

