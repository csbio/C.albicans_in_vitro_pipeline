#-------------------------------------------------
# Raw sequencing reads data clean-up
#-------------------------------------------------
# Last updated October 25th, 2024 by Emily Xiong

# Install required packages (only required upon initial setup).
install.packages("reshape2") 
install.packages("dplyr")
install.packages("tidyr")

# Load packages. 
library(reshape2) # for table management
library(dplyr) # for table management
library(tidyr)

#set working directory to file with data 
setwd("C.albicans_in_vitro_pipeline/STAR_Protocols_code_data
      /assign_example/output/MappingOutput-2024-10-25/")
### FOR WRITING: setwd("/Volumes/T7/assign_example/output/MappingOutput-2024-10-25")
#load data
data <- read.delim("merged_results.txt", header = TRUE)
orf19 <- data$orf19
feature <- data$Feature
common <- data$Common
description <- data$Description
plate <- data$GRACE.Plate
position <- data$GRACE.Position

#making regular expressions for sequencing technical replicates
regex_YNB1_AUP <- "(YNB1_A_S\\d+_L00\\d+_UP)"
regex_YNB1_ADN <- "(YNB1_A_S\\d+_L00\\d+_DOWN)"
regex_YNB1_BUP <- "(YNB1_B_S\\d+_L00\\d+_UP)"
regex_YNB1_BDN <- "(YNB1_B_S\\d+_L00\\d+_DOWN)"
regex_YNB1_CUP <- "(YNB1_C_S\\d+_L00\\d+_UP)"
regex_YNB1_CDN <- "(YNB1_C_S\\d+_L00\\d+_DOWN)"

regex_YNB1_DOX_AUP <- "(YNB1_DOX_A_S\\d+_L00\\d+_UP)"
regex_YNB1_DOX_ADN <- "(YNB1_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_YNB1_DOX_BUP <- "(YNB1_DOX_B_S\\d+_L00\\d+_UP)"
regex_YNB1_DOX_BDN <- "(YNB1_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_YNB1_DOX_CUP <- "(YNB1_DOX_C_S\\d+_L00\\d+_UP)"
regex_YNB1_DOX_CDN <- "(YNB1_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_YNB2_AUP <- "(YNB2_A_S\\d+_L00\\d+_UP)"
regex_YNB2_ADN <- "(YNB2_A_S\\d+_L00\\d+_DOWN)"
regex_YNB2_BUP <- "(YNB2_B_S\\d+_L00\\d+_UP)"
regex_YNB2_BDN <- "(YNB2_B_S\\d+_L00\\d+_DOWN)"
regex_YNB2_CUP <- "(YNB2_C_S\\d+_L00\\d+_UP)"
regex_YNB2_CDN <- "(YNB2_C_S\\d+_L00\\d+_DOWN)"

regex_YNB2_DOX_AUP <- "(YNB2_DOX_A_S\\d+_L00\\d+_UP)"
regex_YNB2_DOX_ADN <- "(YNB2_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_YNB2_DOX_BUP <- "(YNB2_DOX_B_S\\d+_L00\\d+_UP)"
regex_YNB2_DOX_BDN <- "(YNB2_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_YNB2_DOX_CUP <- "(YNB2_DOX_C_S\\d+_L00\\d+_UP)"
regex_YNB2_DOX_CDN <- "(YNB2_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_YNB3_AUP <- "(YNB3_A_S\\d+_L00\\d+_UP)"
regex_YNB3_ADN <- "(YNB3_A_S\\d+_L00\\d+_DOWN)"
regex_YNB3_BUP <- "(YNB3_B_S\\d+_L00\\d+_UP)"
regex_YNB3_BDN <- "(YNB3_B_S\\d+_L00\\d+_DOWN)"
regex_YNB3_CUP <- "(YNB3_C_S\\d+_L00\\d+_UP)"
regex_YNB3_CDN <- "(YNB3_C_S\\d+_L00\\d+_DOWN)"

regex_YNB3_DOX_AUP <- "(YNB3_DOX_A_S\\d+_L00\\d+_UP)"
regex_YNB3_DOX_ADN <- "(YNB3_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_YNB3_DOX_BUP <- "(YNB3_DOX_B_S\\d+_L00\\d+_UP)"
regex_YNB3_DOX_BDN <- "(YNB3_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_YNB3_DOX_CUP <- "(YNB3_DOX_C_S\\d+_L00\\d+_UP)"
regex_YNB3_DOX_CDN <- "(YNB3_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_YNB4_AUP <- "(YNB4_A_S\\d+_L00\\d+_UP)"
regex_YNB4_ADN <- "(YNB4_A_S\\d+_L00\\d+_DOWN)"
regex_YNB4_BUP <- "(YNB4_B_S\\d+_L00\\d+_UP)"
regex_YNB4_BDN <- "(YNB4_B_S\\d+_L00\\d+_DOWN)"
regex_YNB4_CUP <- "(YNB4_C_S\\d+_L00\\d+_UP)"
regex_YNB4_CDN <- "(YNB4_C_S\\d+_L00\\d+_DOWN)"

regex_YNB4_DOX_AUP <- "(YNB4_DOX_A_S\\d+_L00\\d+_UP)"
regex_YNB4_DOX_ADN <- "(YNB4_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_YNB4_DOX_BUP <- "(YNB4_DOX_B_S\\d+_L00\\d+_UP)"
regex_YNB4_DOX_BDN <- "(YNB4_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_YNB4_DOX_CUP <- "(YNB4_DOX_C_S\\d+_L00\\d+_UP)"
regex_YNB4_DOX_CDN <- "(YNB4_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_YNB5_AUP <- "(YNB5_A_S\\d+_L00\\d+_UP)"
regex_YNB5_ADN <- "(YNB5_A_S\\d+_L00\\d+_DOWN)"
regex_YNB5_BUP <- "(YNB5_B_S\\d+_L00\\d+_UP)"
regex_YNB5_BDN <- "(YNB5_B_S\\d+_L00\\d+_DOWN)"
regex_YNB5_CUP <- "(YNB5_C_S\\d+_L00\\d+_UP)"
regex_YNB5_CDN <- "(YNB5_C_S\\d+_L00\\d+_DOWN)"

regex_YNB5_DOX_AUP <- "(YNB5_DOX_A_S\\d+_L00\\d+_UP)"
regex_YNB5_DOX_ADN <- "(YNB5_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_YNB5_DOX_BUP <- "(YNB5_DOX_B_S\\d+_L00\\d+_UP)"
regex_YNB5_DOX_BDN <- "(YNB5_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_YNB5_DOX_CUP <- "(YNB5_DOX_C_S\\d+_L00\\d+_UP)"
regex_YNB5_DOX_CDN <- "(YNB5_DOX_C_S\\d+_L00\\d+_DOWN)"
#--an X is added before the name because R will automatically add that to any headers that start with a number
# this may differ for your samples, check your specific dataframe 
regex_Temp1_AUP <- "(X37_1_A_S\\d+_L00\\d+_UP)"
regex_Temp1_ADN <- "(X37_1_A_S\\d+_L00\\d+_DOWN)"
regex_Temp1_BUP <- "(X37_1_B_S\\d+_L00\\d+_UP)"
regex_Temp1_BDN <- "(X37_1_B_S\\d+_L00\\d+_DOWN)"
regex_Temp1_CUP <- "(X37_1_C_S\\d+_L00\\d+_UP)"
regex_Temp1_CDN <- "(X37_1_C_S\\d+_L00\\d+_DOWN)"

regex_Temp1_DOX_AUP <- "(X37_1_DOX_A_S\\d+_L00\\d+_UP)"
regex_Temp1_DOX_ADN <- "(X37_1_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_Temp1_DOX_BUP <- "(X37_1_DOX_B_S\\d+_L00\\d+_UP)"
regex_Temp1_DOX_BDN <- "(X37_1_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_Temp1_DOX_CUP <- "(X37_1_DOX_C_S\\d+_L00\\d+_UP)"
regex_Temp1_DOX_CDN <- "(X37_1_DOX_C_S\\d+_L00\\d+_DOWN)"
#--
regex_Temp2_AUP <- "(X37_2_A_S\\d+_L00\\d+_UP)"
regex_Temp2_ADN <- "(X37_2_A_S\\d+_L00\\d+_DOWN)"
regex_Temp2_BUP <- "(X37_2_B_S\\d+_L00\\d+_UP)"
regex_Temp2_BDN <- "(X37_2_B_S\\d+_L00\\d+_DOWN)"
regex_Temp2_CUP <- "(X37_2_C_S\\d+_L00\\d+_UP)"
regex_Temp2_CDN <- "(X37_2_C_S\\d+_L00\\d+_DOWN)"

regex_Temp2_DOX_AUP <- "(X37_2_DOX_A_S\\d+_L00\\d+_UP)"
regex_Temp2_DOX_ADN <- "(X37_2_DOX_A_S\\d+_L00\\d+_DOWN)"
regex_Temp2_DOX_BUP <- "(X37_2_DOX_B_S\\d+_L00\\d+_UP)"
regex_Temp2_DOX_BDN <- "(X37_2_DOX_B_S\\d+_L00\\d+_DOWN)"
regex_Temp2_DOX_CUP <- "(X37_2_DOX_C_S\\d+_L00\\d+_UP)"
regex_Temp2_DOX_CDN <- "(X37_2_DOX_C_S\\d+_L00\\d+_DOWN)"


#average sequencing technical replicates
A_YNB1_UP <- rowSums(data[,grep(regex_YNB1_AUP, colnames(data))])
A_YNB1_DN <- rowSums(data[,grep(regex_YNB1_ADN, colnames(data))])
B_YNB1_UP <- rowSums(data[,grep(regex_YNB1_BUP, colnames(data))])
B_YNB1_DN <- rowSums(data[,grep(regex_YNB1_BDN, colnames(data))])
C_YNB1_UP <- rowSums(data[,grep(regex_YNB1_CUP, colnames(data))])
C_YNB1_DN <- rowSums(data[,grep(regex_YNB1_CDN, colnames(data))])

A_YNB1_DOX_UP <- rowSums(data[,grep(regex_YNB1_DOX_AUP, colnames(data))])
A_YNB1_DOX_DN <- rowSums(data[,grep(regex_YNB1_DOX_ADN, colnames(data))])
B_YNB1_DOX_UP <- rowSums(data[,grep(regex_YNB1_DOX_BUP, colnames(data))])
B_YNB1_DOX_DN <- rowSums(data[,grep(regex_YNB1_DOX_BDN, colnames(data))])
C_YNB1_DOX_UP <- rowSums(data[,grep(regex_YNB1_DOX_CUP, colnames(data))])
C_YNB1_DOX_DN <- rowSums(data[,grep(regex_YNB1_DOX_CDN, colnames(data))])
#--
A_YNB2_UP <- rowSums(data[,grep(regex_YNB2_AUP, colnames(data))])
A_YNB2_DN <- rowSums(data[,grep(regex_YNB2_ADN, colnames(data))])
B_YNB2_UP <- rowSums(data[,grep(regex_YNB2_BUP, colnames(data))])
B_YNB2_DN <- rowSums(data[,grep(regex_YNB2_BDN, colnames(data))])
C_YNB2_UP <- rowSums(data[,grep(regex_YNB2_CUP, colnames(data))])
C_YNB2_DN <- rowSums(data[,grep(regex_YNB2_CDN, colnames(data))])

A_YNB2_DOX_UP <- rowSums(data[,grep(regex_YNB2_DOX_AUP, colnames(data))])
A_YNB2_DOX_DN <- rowSums(data[,grep(regex_YNB2_DOX_ADN, colnames(data))])
B_YNB2_DOX_UP <- rowSums(data[,grep(regex_YNB2_DOX_BUP, colnames(data))])
B_YNB2_DOX_DN <- rowSums(data[,grep(regex_YNB2_DOX_BDN, colnames(data))])
C_YNB2_DOX_UP <- rowSums(data[,grep(regex_YNB2_DOX_CUP, colnames(data))])
C_YNB2_DOX_DN <- rowSums(data[,grep(regex_YNB2_DOX_CDN, colnames(data))])
#--
A_YNB3_UP <- rowSums(data[,grep(regex_YNB3_AUP, colnames(data))])
A_YNB3_DN <- rowSums(data[,grep(regex_YNB3_ADN, colnames(data))])
B_YNB3_UP <- rowSums(data[,grep(regex_YNB3_BUP, colnames(data))])
B_YNB3_DN <- rowSums(data[,grep(regex_YNB3_BDN, colnames(data))])
C_YNB3_UP <- rowSums(data[,grep(regex_YNB3_CUP, colnames(data))])
C_YNB3_DN <- rowSums(data[,grep(regex_YNB3_CDN, colnames(data))])

A_YNB3_DOX_UP <- rowSums(data[,grep(regex_YNB3_DOX_AUP, colnames(data))])
A_YNB3_DOX_DN <- rowSums(data[,grep(regex_YNB3_DOX_ADN, colnames(data))])
B_YNB3_DOX_UP <- rowSums(data[,grep(regex_YNB3_DOX_BUP, colnames(data))])
B_YNB3_DOX_DN <- rowSums(data[,grep(regex_YNB3_DOX_BDN, colnames(data))])
C_YNB3_DOX_UP <- rowSums(data[,grep(regex_YNB3_DOX_CUP, colnames(data))])
C_YNB3_DOX_DN <- rowSums(data[,grep(regex_YNB3_DOX_CDN, colnames(data))])
#--
A_YNB4_UP <- rowSums(data[,grep(regex_YNB4_AUP, colnames(data))])
A_YNB4_DN <- rowSums(data[,grep(regex_YNB4_ADN, colnames(data))])
B_YNB4_UP <- rowSums(data[,grep(regex_YNB4_BUP, colnames(data))])
B_YNB4_DN <- rowSums(data[,grep(regex_YNB4_BDN, colnames(data))])
C_YNB4_UP <- rowSums(data[,grep(regex_YNB4_CUP, colnames(data))])
C_YNB4_DN <- rowSums(data[,grep(regex_YNB4_CDN, colnames(data))])

A_YNB4_DOX_UP <- rowSums(data[,grep(regex_YNB4_DOX_AUP, colnames(data))])
A_YNB4_DOX_DN <- rowSums(data[,grep(regex_YNB4_DOX_ADN, colnames(data))])
B_YNB4_DOX_UP <- rowSums(data[,grep(regex_YNB4_DOX_BUP, colnames(data))])
B_YNB4_DOX_DN <- rowSums(data[,grep(regex_YNB4_DOX_BDN, colnames(data))])
C_YNB4_DOX_UP <- rowSums(data[,grep(regex_YNB4_DOX_CUP, colnames(data))])
C_YNB4_DOX_DN <- rowSums(data[,grep(regex_YNB4_DOX_CDN, colnames(data))])
#--
A_YNB5_UP <- rowSums(data[,grep(regex_YNB5_AUP, colnames(data))])
A_YNB5_DN <- rowSums(data[,grep(regex_YNB5_ADN, colnames(data))])
B_YNB5_UP <- rowSums(data[,grep(regex_YNB5_BUP, colnames(data))])
B_YNB5_DN <- rowSums(data[,grep(regex_YNB5_BDN, colnames(data))])
C_YNB5_UP <- rowSums(data[,grep(regex_YNB5_CUP, colnames(data))])
C_YNB5_DN <- rowSums(data[,grep(regex_YNB5_CDN, colnames(data))])

A_YNB5_DOX_UP <- rowSums(data[,grep(regex_YNB5_DOX_AUP, colnames(data))])
A_YNB5_DOX_DN <- rowSums(data[,grep(regex_YNB5_DOX_ADN, colnames(data))])
B_YNB5_DOX_UP <- rowSums(data[,grep(regex_YNB5_DOX_BUP, colnames(data))])
B_YNB5_DOX_DN <- rowSums(data[,grep(regex_YNB5_DOX_BDN, colnames(data))])
C_YNB5_DOX_UP <- rowSums(data[,grep(regex_YNB5_DOX_CUP, colnames(data))])
C_YNB5_DOX_DN <- rowSums(data[,grep(regex_YNB5_DOX_CDN, colnames(data))])
#--
A_Temp1_UP <- rowSums(data[,grep(regex_Temp1_AUP, colnames(data))])
A_Temp1_DN <- rowSums(data[,grep(regex_Temp1_ADN, colnames(data))])
B_Temp1_UP <- rowSums(data[,grep(regex_Temp1_BUP, colnames(data))])
B_Temp1_DN <- rowSums(data[,grep(regex_Temp1_BDN, colnames(data))])
C_Temp1_UP <- rowSums(data[,grep(regex_Temp1_CUP, colnames(data))])
C_Temp1_DN <- rowSums(data[,grep(regex_Temp1_CDN, colnames(data))])

A_Temp1_DOX_UP <- rowSums(data[,grep(regex_Temp1_DOX_AUP, colnames(data))])
A_Temp1_DOX_DN <- rowSums(data[,grep(regex_Temp1_DOX_ADN, colnames(data))])
B_Temp1_DOX_UP <- rowSums(data[,grep(regex_Temp1_DOX_BUP, colnames(data))])
B_Temp1_DOX_DN <- rowSums(data[,grep(regex_Temp1_DOX_BDN, colnames(data))])
C_Temp1_DOX_UP <- rowSums(data[,grep(regex_Temp1_DOX_CUP, colnames(data))])
C_Temp1_DOX_DN <- rowSums(data[,grep(regex_Temp1_DOX_CDN, colnames(data))])
#--
A_Temp2_UP <- rowSums(data[,grep(regex_Temp2_AUP, colnames(data))])
A_Temp2_DN <- rowSums(data[,grep(regex_Temp2_ADN, colnames(data))])
B_Temp2_UP <- rowSums(data[,grep(regex_Temp2_BUP, colnames(data))])
B_Temp2_DN <- rowSums(data[,grep(regex_Temp2_BDN, colnames(data))])
C_Temp2_UP <- rowSums(data[,grep(regex_Temp2_CUP, colnames(data))])
C_Temp2_DN <- rowSums(data[,grep(regex_Temp2_CDN, colnames(data))])

A_Temp2_DOX_UP <- rowSums(data[,grep(regex_Temp2_DOX_AUP, colnames(data))])
A_Temp2_DOX_DN <- rowSums(data[,grep(regex_Temp2_DOX_ADN, colnames(data))])
B_Temp2_DOX_UP <- rowSums(data[,grep(regex_Temp2_DOX_BUP, colnames(data))])
B_Temp2_DOX_DN <- rowSums(data[,grep(regex_Temp2_DOX_BDN, colnames(data))])
C_Temp2_DOX_UP <- rowSums(data[,grep(regex_Temp2_DOX_CUP, colnames(data))])
C_Temp2_DOX_DN <- rowSums(data[,grep(regex_Temp2_DOX_CDN, colnames(data))])

#generate dataframe with all summed barcode counts 
sum_data <- cbind(orf19, feature, common, description, plate, position,
                  A_YNB1_UP, A_YNB1_DN,
                  B_YNB1_UP, B_YNB1_DN,
                  C_YNB1_UP, C_YNB1_DN,
                  A_YNB1_DOX_UP, A_YNB1_DOX_DN,
                  B_YNB1_DOX_UP, B_YNB1_DOX_DN,
                  C_YNB1_DOX_UP, C_YNB1_DOX_DN,
                  A_YNB2_UP, A_YNB2_DN,
                  B_YNB2_UP, B_YNB2_DN,
                  C_YNB2_UP, C_YNB2_DN,
                  A_YNB2_DOX_UP, A_YNB2_DOX_DN,
                  B_YNB2_DOX_UP, B_YNB2_DOX_DN,
                  C_YNB2_DOX_UP, C_YNB2_DOX_DN,
                  A_YNB3_UP, A_YNB3_DN,
                  B_YNB3_UP, B_YNB3_DN,
                  C_YNB3_UP, C_YNB3_DN,
                  A_YNB3_DOX_UP, A_YNB3_DOX_DN,
                  B_YNB3_DOX_UP, B_YNB3_DOX_DN,
                  C_YNB3_DOX_UP, C_YNB3_DOX_DN,
                  A_YNB4_UP, A_YNB4_DN,
                  B_YNB4_UP, B_YNB4_DN,
                  C_YNB4_UP, C_YNB4_DN,
                  A_YNB4_DOX_UP, A_YNB4_DOX_DN,
                  B_YNB4_DOX_UP, B_YNB4_DOX_DN,
                  C_YNB4_DOX_UP, C_YNB4_DOX_DN,
                  A_YNB5_UP, A_YNB5_DN,
                  B_YNB5_UP, B_YNB5_DN,
                  C_YNB5_UP, C_YNB5_DN,
                  A_YNB5_DOX_UP, A_YNB5_DOX_DN,
                  B_YNB5_DOX_UP, B_YNB5_DOX_DN,
                  C_YNB5_DOX_UP, C_YNB5_DOX_DN,
                  A_Temp1_UP, A_Temp1_DN,
                  B_Temp1_UP, B_Temp1_DN,
                  C_Temp1_UP, C_Temp1_DN,
                  A_Temp1_DOX_UP, A_Temp1_DOX_DN,
                  B_Temp1_DOX_UP, B_Temp1_DOX_DN,
                  C_Temp1_DOX_UP, C_Temp1_DOX_DN,
                  A_Temp2_UP, A_Temp2_DN,
                  B_Temp2_UP, B_Temp2_DN,
                  C_Temp2_UP, C_Temp2_DN,
                  A_Temp2_DOX_UP, A_Temp2_DOX_DN,
                  B_Temp2_DOX_UP, B_Temp2_DOX_DN,
                  C_Temp2_DOX_UP, C_Temp2_DOX_DN)
write.csv(sum_data, file = "sum_data_numeric.csv")

#data pruning to remove strains with very low/inconsistent read counts
sum_data_numeric <- read.csv("sum_data_numeric.csv", header = T)
sum_data_numeric$sum <- rowSums(sum_data_numeric[ , c(8:90)], na.rm=TRUE) 
# c(X:Y) replace X with column # of first condition and Y with column # of last
sum_data_purge <- sum_data_numeric[sum_data_numeric$sum > 150, ]

#export final csv file
write.csv(sum_data_purge, file = "sum_data.csv")
