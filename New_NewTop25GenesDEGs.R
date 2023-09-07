##########################################################################################################################################################################
######Select the top 25 genes with the highest Log2Fold change c--> ignore sing (absoulute values),all genes sould be expressed in the small intestine acording to (NCBI)#
##########################################################################################################################################################################

#clean workspace
rm(list=ls())

#set working directory
setwd("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput")

install.packages("dplyr")
library(dplyr)
library("writexl")

########################################
######Read TopTables and make DEGS######
########################################

#Read duodenum tables
toptable_SvsCinLowAbs_D <- read.table("toptable_duodenum_multi-level_Campesterol_chol_Stanol_VS_Control_LowAbsorber.csv", header=TRUE,sep=",")
toptable_SvsCinHighAbs_D <- read.table("toptable_duodenum_multi-level_Campesterol_chol_Stanol_VS_Control_HighAbsorber.csv", header=TRUE,sep=",")

#Make DEGS
#Determine DEG in Duodenum
DEG_SvsCinLowAbs_D <- toptable_SvsCinLowAbs_D[toptable_SvsCinLowAbs_D$P.Value< 0.05, c(1,2,5,8)]
DEG_SvsCinHighAbs_D <- toptable_SvsCinHighAbs_D[toptable_SvsCinHighAbs_D$P.Value< 0.05, c(1,2,5,8)]
# DEG_Diff_D <- toptable_Diff_D[toptable_Diff_D$P.Value< 0.05, c(1,2,5,8)] #E3E418FF
###1:ID,2:logFC,5:P.value,8:SYMBOL.

# specify if a gene is up or down regulated
DEG_SvsCinLowAbs_D$class <- ifelse(DEG_SvsCinLowAbs_D$logFC>0,"up","down")
DEG_SvsCinHighAbs_D$class <- ifelse(DEG_SvsCinHighAbs_D$logFC>0,"up","down")

#Read jejunum tables
toptable_SvsCinLowAbs_J <- read.table("toptable_jejunum_multi-level_Campesterol_chol_Stanol_VS_Control_LowAbsorber.csv", header=TRUE,sep=",")
toptable_SvsCinHighAbs_J <- read.table("toptable_jejunum_multi-level_Campesterol_chol_Stanol_VS_Control_HighAbsorber.csv", header=TRUE,sep=",")

#Make DEGs
#Determine DEG in Jejunum
DEG_SvsCinLowAbs_J <- toptable_SvsCinLowAbs_J[toptable_SvsCinLowAbs_J$P.Value< 0.05, c(1,2,5,8)]
DEG_SvsCinHighAbs_J <- toptable_SvsCinHighAbs_J[toptable_SvsCinHighAbs_J$P.Value< 0.05, c(1,2,5,8)]

# DEG_Diff_J <- toptable_Diff_J[toptable_Diff_J$P.Value< 0.05, c(1,2,5,8)] #E3E418FF

# specify if a gene is up or down regulated
DEG_SvsCinLowAbs_J$class <- ifelse(DEG_SvsCinLowAbs_J$logFC>0,"up","down")
DEG_SvsCinHighAbs_J$class <- ifelse(DEG_SvsCinHighAbs_J$logFC>0,"up","down")


################################################################################################
################################################################################################
#1-Selct the first 25 genes from DEG_SvsCinLowAbs_D# Top (absoulute values)
################################################################################################

# Filter by Top 25 absolute values in column "logFC"
Top25_Abso_DEGS_SvsCinLowAbs_D <- DEG_SvsCinLowAbs_D[order(abs(DEG_SvsCinLowAbs_D$logFC), decreasing = TRUE), ][1:25, ]

#Save in excel file#
write_xlsx(Top25_Abso_DEGS_SvsCinLowAbs_D,"/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/25TopGenes/Top25_Abso_DEGS_SvsCinLowAbs_D.xlsx")


################################################################################
#Now we would like to filter the same 25 genes from the high absorbers_duodenum#
################################################################################

#Remove  unnecessary columns
toptable_SvsCinHighAbs_D <- toptable_SvsCinHighAbs_D[c(1,2,5,8)]

#Add up or down
toptable_SvsCinHighAbs_D$class <- ifelse(toptable_SvsCinHighAbs_D$logFC>0,"up","down")

#Filter the 25 genes from the HighAbsorbers_D#
Filtered25genesOfHighAbsorbers_D <- subset(toptable_SvsCinHighAbs_D, ID %in% Top25_Abso_DEGS_SvsCinLowAbs_D$ID)

#Reorder the table to much the Low absorber table#
Filtered25genesOfHighAbsorbers_D <- Filtered25genesOfHighAbsorbers_D [match(Top25_Abso_DEGS_SvsCinLowAbs_D$ID, Filtered25genesOfHighAbsorbers_D$ID),]

#Save in excel file#
write_xlsx(Filtered25genesOfHighAbsorbers_D,"/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/25TopGenes/Filtered25genesOfHighAbsorbers_D.xlsx")

################################################################################################
################################################################################################
#2-Selct the first 25 genes from DEG_SvsCinHighAbs_D# Top (absoulute values)
################################################################################################

# Filter by Top 25 absolute values in column "logFC"
Top25_Abso_DEGS_SvsCinHighAbs_D <- DEG_SvsCinHighAbs_D[order(abs(DEG_SvsCinHighAbs_D$logFC), decreasing = TRUE), ][1:25, ]

#Save in excel file#
write_xlsx(Top25_Abso_DEGS_SvsCinHighAbs_D,"/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/25TopGenes/Top25_Abso_DEGS_SvsCinHighAbs_D.xlsx")

################################################################################
#Now we would like to filter the same genes from the Low absorbers_duodenum#
################################################################################

#Remove  unnecessary columns
toptable_SvsCinLowAbs_D <- toptable_SvsCinLowAbs_D[c(1,2,5,8)]

#Add up or down
toptable_SvsCinLowAbs_D$class <- ifelse(toptable_SvsCinLowAbs_D$logFC>0,"up","down")

#Filter the 25 genes from the HighAbsorbers_D#
Filtered25genesOfLowAbsorbers_D <- subset(toptable_SvsCinLowAbs_D, ID %in% Top25_Abso_DEGS_SvsCinHighAbs_D$ID)

#Reorder the table to much the Low absorber table#
Filtered25genesOfLowAbsorbers_D <- Filtered25genesOfLowAbsorbers_D [match(Top25_Abso_DEGS_SvsCinHighAbs_D$ID, Filtered25genesOfLowAbsorbers_D$ID),]

#Save in excel file#
write_xlsx(Filtered25genesOfLowAbsorbers_D,"/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/25TopGenes/Filtered25genesOfLowbsorbers_D.xlsx")

################################################################################################
################################################################################################
#3-Selct the first 25 genes from DEG_SvsCinLowAbs_J# Top (absoulute values)
################################################################################################

# Filter by Top 25 absolute values in column "logFC"
Top25_Abso_DEGS_SvsCinLowAbs_J <- DEG_SvsCinLowAbs_J[order(abs(DEG_SvsCinLowAbs_J$logFC), decreasing = TRUE), ][1:25, ]

#Save in excel file#
write_xlsx(Top25_Abso_DEGS_SvsCinLowAbs_J,"/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/25TopGenes/Top25_Abso_DEGS_SvsCinLowAbs_J.xlsx")

################################################################################
#Now we would like to filter the same 27 genes from the high absorbers_jejunum#
################################################################################

#Remove  unnecessary columns
toptable_SvsCinHighAbs_J <- toptable_SvsCinHighAbs_J[c(1,2,5,8)]

#Add up or down
toptable_SvsCinHighAbs_J$class <- ifelse(toptable_SvsCinHighAbs_J$logFC>0,"up","down")

#Filter the 25 genes from the HighAbsorbers_D#
Filtered25genesOfHighAbsorbers_J <- subset(toptable_SvsCinHighAbs_J, ID %in% Top25_Abso_DEGS_SvsCinLowAbs_J$ID)

#Reorder the table to much the Low absorber table#
Filtered25genesOfHighAbsorbers_J <- Filtered25genesOfHighAbsorbers_J [match(Top25_Abso_DEGS_SvsCinLowAbs_J$ID, Filtered25genesOfHighAbsorbers_J$ID),]

#Save in excel file#
write_xlsx(Filtered25genesOfHighAbsorbers_J,"/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/25TopGenes/Filtered25genesOfHighAbsorbers_J.xlsx")

################################################################################################
################################################################################################
#4-Selct the first 25 genes from DEG_SvsCinHighAbs_J# Top (absoulute values)
################################################################################################

# Filter by Top 25 absolute values in column "logFC"
Top25_Abso_DEGS_SvsCinHighAbs_J <- DEG_SvsCinHighAbs_J[order(abs(DEG_SvsCinHighAbs_J$logFC), decreasing = TRUE), ][1:25, ]

#Save in excel file#
write_xlsx(Top25_Abso_DEGS_SvsCinHighAbs_J,"/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/25TopGenes/Top25_Abso_DEGS_SvsCinHighAbs_J.xlsx")

################################################################################
#Now we would like to filter the same genes from the Low absorbers_duodenum#
################################################################################

#Remove  unnecessary columns
toptable_SvsCinLowAbs_J <- toptable_SvsCinLowAbs_J[c(1,2,5,8)]

#Add up or down
toptable_SvsCinLowAbs_J$class <- ifelse(toptable_SvsCinLowAbs_J$logFC>0,"up","down")

#Filter the 25 genes from the HighAbsorbers_D#
Filtered25genesOfLowAbsorbers_J <- subset(toptable_SvsCinLowAbs_J, ID %in% Top25_Abso_DEGS_SvsCinHighAbs_J$ID)

#Reorder the table to much the Low absorber table#
Filtered25genesOfLowAbsorbers_J <- Filtered25genesOfLowAbsorbers_J [match(Top25_Abso_DEGS_SvsCinHighAbs_J$ID, Filtered25genesOfLowAbsorbers_J$ID),]

#Save in excel file#
write_xlsx(Filtered25genesOfLowAbsorbers_J,"/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/25TopGenes/Filtered25genesOfLowbsorbers_J.xlsx")

