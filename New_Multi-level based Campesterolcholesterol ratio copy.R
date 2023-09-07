########################################################
######Starting functions for Multi-level analysis#######
########################################################

#clean workspace TEST
rm(list=ls())

#set working directory
setwd("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC")

#Load packages
#library(heatmap.plus)
#library(Heatplus)
library(dplyr)
library(abind)
library(Biobase)
library(limma)
library(tibble)
library(enrichplot)
library(VennDiagram)
library(EnhancedVolcano)
library(clusterProfiler)
library(plyr)
library(showtext)
library(tidyverse)
install.packages("stringr")    
library("stringr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


devtools::install_github("r-lib/xml2")



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu95av2.db")
BiocManager::install("S4Vectors")



########################################################
###Multi-level based on Campesterol/cholesterol ratio###
########################################################

#Read duodenum en jejunum tables
db_mLOG2 <- read.table("input/Duodenum_blauw_log2.txt", header=TRUE)
dg_mLOG2 <- read.table("input/Duodenum_geel_log2.txt", header=TRUE)
jb_mLOG2 <- read.table("input/Jejunum_blauw_log2.txt", header=TRUE)
jg_mLOG2 <- read.table("input/Jejunum_geel_log2.txt", header=TRUE)

#create a matrix out of the data frames
db_mLOG2 <- as.matrix(db_mLOG2)
dg_mLOG2 <- as.matrix(dg_mLOG2)
jb_mLOG2 <- as.matrix(jb_mLOG2)
jg_mLOG2 <- as.matrix(jg_mLOG2)

#Heatmap creation Duodenum and Jejunum separtly at control condition (Baseline) 
# Define colors used in the heatmap
my_palette <- colorRampPalette(c("blue",
                                 "white",
                                 "red"))

#1-heatmap creation Duodenum control condition
png('output/heatmap_Duodenum_control_synth.png')
heatmap.2(dg_mLOG2[,-1], 
          trace = "none",               # turn off trace lines from heatmap
          cexCol = 1.0,
          cexRow = 0.7,
          col = my_palette(25),          # use my colour scheme with 25 levels
          main = "Duodenum control condition")    # add title
dev.off()

#2-heatmap creation Jejunum control condition
png('output/heatmap_Jejunum_control.png')
heatmap.2(jg_mLOG2[,-1], 
          trace = "none",               # turn off trace lines from heatmap
          cexCol = 1.0,
          cexRow = 0.7,
          col = my_palette(25),          # use my colour scheme with 25 levels
          main = "Jejunum control condition")    # add title
dev.off()

############################################
#create experimental group factor variable
#Make a list of three times High absorber and three times Low absorber
pd_camp_c <- cbind(Type=c(rep("High_absorber",3),rep("Low_absorber",3),rep("High_absorber",3),rep("Low_absorber",3)))
pj_camp_c <- cbind(Type=c(rep("High_absorber",3),rep("Low_absorber",3),rep("High_absorber",3),rep("Low_absorber",3)))

#Used samples duodenum: High absorbers (p13, p1, p14), Low absorbers (p20, p8, p15)
#Used samples jejunum: High absorbers (p13, p1, p3), Low absorbers (p8, p15, p11)

#Change rownames to the desired samples (samples selected based on campe/chol)
rownames(pd_camp_c) =c("blauw_p13","blauw_p1","blauw_p14","blauw_p20","blauw_p8","blauw_p15","geel_p13","geel_p1","geel_p14","geel_p20","geel_p8","geel_p15")
rownames(pj_camp_c) =c("blauw_p13","blauw_p1","blauw_p3","blauw_p8","blauw_p15","blauw_p11","geel_p13","geel_p1","geel_p3","geel_p8","geel_p15","geel_p11")

#Add corresponding treatment type
# Geel = Control
# Blauw = Plant sterols

Treatment <- c(rep("Stanol",6),rep("Control",6))
pd_camp_c<- cbind(pd_camp_c, Treatment)
pj_camp_c<- cbind(pj_camp_c, Treatment)

#Add corresponding SibShip type 
SibShip <- c(1,2,3,4,5,6,1,2,3,4,5,6)
pd_camp_c<- cbind(pd_camp_c, SibShip)
pj_camp_c<- cbind(pj_camp_c, SibShip)

#Change to DF --> by this we will have for each treatment the corresponding control (They get the same number)
pd_camp_c <- as.data.frame(pd_camp_c)
pj_camp_c <- as.data.frame(pj_camp_c)

#Select gene expression data based on selected samples
# Type_db_camp_c <- db_mLOG2[,c(11,3,12,2,6,10)]
# Type_dg_camp_c <- dg_mLOG2[,c(11,4,12,2,6,10)]
# Type_jb_camp_c <- jb_mLOG2[,c(11,2,4,5,9,8)]
# Type_jg_camp_c <- jg_mLOG2[,c(11,3,6,7,13,10)]

Type_db_camp_c <- db_mLOG2[,c("blauw_p13","blauw_p1","blauw_p14","blauw_p20","blauw_p8","blauw_p15")]
Type_dg_camp_c <- dg_mLOG2[,c("geel_p13","geel_p1","geel_p14","geel_p20","geel_p8","geel_p15")]
Type_jb_camp_c <- jb_mLOG2[,c("blauw_p13","blauw_p1","blauw_p3","blauw_p8","blauw_p15","blauw_p11")]
Type_jg_camp_c <- jg_mLOG2[,c("geel_p13","geel_p1","geel_p3","geel_p8","geel_p15","geel_p11")]

#Combine RNA expression files (of the treatment and control)<- volgorde verandert niet
Type_d_camp_c <- abind(Type_db_camp_c, Type_dg_camp_c)
Type_j_camp_c <- abind(Type_jb_camp_c, Type_jg_camp_c)

#ADD Genes IDs to the list of combined treatment and control:
rownames(Type_d_camp_c) <- db_mLOG2[,"ID"]
rownames(Type_j_camp_c) <- jb_mLOG2[,"ID"]

#Check class
class(Type_d_camp_c) #matrix and array
class(pd_camp_c) #data.frame

#Duodenum design based on camp/c ratio
InterD <- factor(paste(pd_camp_c$Treatment,pd_camp_c$Type,sep="."))
designD <- model.matrix(~0+InterD)
colnames(designD) <- levels(InterD)
#Jejunum design based on camp/c ratio
InterJ <- factor(paste(pj_camp_c$Treatment,pj_camp_c$Type,sep="."))
designJ <- model.matrix(~0+InterJ)
colnames(designJ) <- levels(InterJ)

#Duodenum correlation based on camp/c ratio
corfitD <- duplicateCorrelation(Type_d_camp_c,designD,block=pd_camp_c$SibShip)
corfitD$consensus #0.3712276
#Jejunum correlation based on camp/c ratio
corfitJ <- duplicateCorrelation(Type_j_camp_c,designJ,block=pj_camp_c$SibShip)
corfitJ$consensus #0.2924498

#Duodenum fit based on camp/c ratio
fitD <- lmFit(Type_d_camp_c,designD,block=pd_camp_c$SibShip,correlation=corfitD$consensus)
#Jejunum fit based on camp/c ratio
fitJ <- lmFit(Type_j_camp_c,designJ,block=pj_camp_c$SibShip,correlation=corfitJ$consensus)

#Duodenum cont.matrix based on camp/c ratio
cont.matrix_D <- makeContrasts(
  SvsCinHighAbs=Stanol.High_absorber-Control.High_absorber,
  SvsCinLowAbs=Stanol.Low_absorber-Control.Low_absorber,
  Diff=(Stanol.High_absorber-Control.High_absorber)-(Stanol.Low_absorber-Control.Low_absorber),
  levels=designD)
#Jejunum cont.matrix based on camp/c ratio
cont.matrix_J <- makeContrasts(
  SvsCinHighAbs=Stanol.High_absorber-Control.High_absorber,
  SvsCinLowAbs=Stanol.Low_absorber-Control.Low_absorber,
  Diff=(Stanol.High_absorber-Control.High_absorber)-(Stanol.Low_absorber-Control.Low_absorber),
  levels=designJ)

#Duodenum fit2 based on camp/c ratio
fitD2 <- contrasts.fit(fitD, cont.matrix_D)
fitD2 <- eBayes(fitD2)
#Jejunum fit2 based on camp/c ratio
fitJ2 <- contrasts.fit(fitJ, cont.matrix_J)
fitJ2 <- eBayes(fitJ2)

#Toptable function for Duodenum
toptable_SvsCinHighAbs_D <- topTable(fitD2, adjust.method="BH", coef="SvsCinHighAbs", number=Inf)
toptable_SvsCinLowAbs_D <- topTable(fitD2, adjust.method="BH", coef="SvsCinLowAbs", number=Inf)

#Toptable function for Jejunum
toptable_SvsCinHighAbs_J <- topTable(fitJ2, adjust.method="BH", coef="SvsCinHighAbs", number=Inf)
toptable_SvsCinLowAbs_J <- topTable(fitJ2, adjust.method="BH", coef="SvsCinLowAbs", number=Inf)

#Add rownames as column
toptable_SvsCinHighAbs_D <- rownames_to_column(toptable_SvsCinHighAbs_D, var = "ID")
toptable_SvsCinLowAbs_D <- rownames_to_column(toptable_SvsCinLowAbs_D, var = "ID")

toptable_SvsCinHighAbs_J <- rownames_to_column(toptable_SvsCinHighAbs_J, var = "ID")
toptable_SvsCinLowAbs_J <- rownames_to_column(toptable_SvsCinLowAbs_J, var = "ID")

# Reannotated the EntrezGeneIDs
ids_D <- bitr(toptable_SvsCinHighAbs_D$ID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
ids_J <- bitr(toptable_SvsCinHighAbs_J$ID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")

# Combine toptable and gene symbol
toptable_SvsCinHighAbs_D <- merge(toptable_SvsCinHighAbs_D, ids_D, by.x="ID", by.y="ENTREZID", all=T)
toptable_SvsCinLowAbs_D <- merge(toptable_SvsCinLowAbs_D, ids_D, by.x="ID", by.y="ENTREZID", all=T)

toptable_SvsCinHighAbs_J <- merge(toptable_SvsCinHighAbs_J, ids_J, by.x="ID", by.y="ENTREZID", all=T)
toptable_SvsCinLowAbs_J <- merge(toptable_SvsCinLowAbs_J, ids_J, by.x="ID", by.y="ENTREZID", all=T)

# To determine the number of NAs in the Symbol column
#1-HighAbs_D
Duodenum_High_NA <- sum(is.na(toptable_SvsCinHighAbs_D$SYMBOL))
Duodenum_High_NA
toptable_SvsCinHighAbs_D_NA <- toptable_SvsCinHighAbs_D %>% filter(is.na(SYMBOL))
dim(toptable_SvsCinHighAbs_D_NA)

#2-LowAbs_D
Duodenum_Low_NA <- sum(is.na(toptable_SvsCinLowAbs_D$SYMBOL))
Duodenum_Low_NA
toptable_SvsCinLowAbs_D_NA <- toptable_SvsCinLowAbs_D %>% filter(is.na(SYMBOL))
dim(toptable_SvsCinLowAbs_D_NA)

#3-HighAbs_J
Jejunum_High_NA <- sum(is.na(toptable_SvsCinHighAbs_J$SYMBOL))
Jejunum_High_NA
toptable_SvsCinHighAbs_J_NA <- toptable_SvsCinHighAbs_J %>% filter(is.na(SYMBOL))
dim(toptable_SvsCinHighAbs_J_NA)

#4-LowAbs_J
Jejunum_Low_NA <- sum(is.na(toptable_SvsCinLowAbs_J$SYMBOL))
Jejunum_Low_NA
toptable_SvsCinLowAbs_J_NA <- toptable_SvsCinLowAbs_J %>% filter(is.na(SYMBOL))
dim(toptable_SvsCinLowAbs_J_NA)

#Write NA of all toptable
write.table(toptable_SvsCinHighAbs_D_NA, file = "output/NA/toptable_SvsCinHighAbs_D_NA.csv", sep=",", row.names=FALSE)
write.table(toptable_SvsCinLowAbs_D_NA, file = "output/NA/toptable_SvsCinLowAbs_D_NA.csv", sep=",", row.names=FALSE)

write.table(toptable_SvsCinHighAbs_J_NA, file = "output/NA/toptable_toptable_SvsCinHighAbs_J_NA.csv", sep=",", row.names=FALSE)
write.table(toptable_SvsCinLowAbs_J_NA, file = "output/NA/toptable_toptable_SvsCinLowAbs_J_NA.csv", sep=",", row.names=FALSE)

#Remove NA and replace with the new ID and if the ID = 0 remove it from the list

#1-High absorbers duodenum
#For genes that are replaced new ID (Current ID) replace that in the toptable and genes that were removed from NCBI have a 0, so remove all the zeros:
New_IDs_HighAbs_D <- read.delim("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NA/gene_result_HighAbsorbers_D.txt", header = TRUE, sep = "\t", dec = ".")

# Replace NA IDs with new IDs and remove rows with new ID of 0
merged_IDs_HighAbs_D <- merge(toptable_SvsCinHighAbs_D, New_IDs_HighAbs_D, by = "ID", all.x = TRUE)
merged_IDs_HighAbs_D$CurrentID <- ifelse(is.na(merged_IDs_HighAbs_D$CurrentID), merged_IDs_HighAbs_D$ID, merged_IDs_HighAbs_D$CurrentID)
merged_IDs_HighAbs_D <- merged_IDs_HighAbs_D[merged_IDs_HighAbs_D$CurrentID != 0,]

#Remove ID and Symbol coloumns, as we will use the new ID and call the new Symbols
merged_IDs_HighAbs_D$ID <- NULL
merged_IDs_HighAbs_D$SYMBOL <- NULL

#Change the header (CurrentID) to ID
colnames(merged_IDs_HighAbs_D)[7] <- "ID"

# Reannotated the EntrezGeneIDs
ids_D <- bitr(merged_IDs_HighAbs_D$ID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")

# Combine toptable and gene symbol
toptable_SvsCinHighAbs_D <- merge(merged_IDs_HighAbs_D, ids_D, by.x="ID", by.y="ENTREZID", all=T)

#2-Low absorbers duodenum
#For genes that are replaced new ID (Current ID) replace that in the toptable and genes that were removed from NCBI have a 0, so remove all the zeros:
New_IDs_LowAbs_D <- read.delim("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NA/gene_result_LowAbsorbers_D.txt", header = TRUE, sep = "\t", dec = ".")

# Replace NA IDs with new IDs and remove rows with new ID of 0
merged_IDs_LowAbs_D <- merge(toptable_SvsCinLowAbs_D, New_IDs_LowAbs_D, by = "ID", all.x = TRUE)
merged_IDs_LowAbs_D$CurrentID <- ifelse(is.na(merged_IDs_LowAbs_D$CurrentID), merged_IDs_LowAbs_D$ID, merged_IDs_LowAbs_D$CurrentID)
merged_IDs_LowAbs_D <- merged_IDs_LowAbs_D[merged_IDs_LowAbs_D$CurrentID != 0,]

#Remove ID and Symbol coloumns, as we will use the new ID and call the new Symbols
merged_IDs_LowAbs_D$ID <- NULL
merged_IDs_LowAbs_D$SYMBOL <- NULL

#Change the header (CurrentID) to ID
colnames(merged_IDs_LowAbs_D)[7] <- "ID"

# Reannotated the EntrezGeneIDs
ids_D <- bitr(merged_IDs_LowAbs_D$ID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")

# Combine toptable and gene symbol
toptable_SvsCinLowAbs_D <- merge(merged_IDs_LowAbs_D, ids_D, by.x="ID", by.y="ENTREZID", all=T)

#3-High absorbers Jejunum
#For genes that are replaced new ID (Current ID) replace that in the toptable and genes that were removed from NCBI have a 0, so remove all the zeros:
New_IDs_HighAbs_J <- read.delim("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NA/gene_result_HighAbsorbers_J.txt", header = TRUE, sep = "\t", dec = ".")

# Replace NA IDs with new IDs and remove rows with new ID of 0
merged_IDs_HighAbs_J <- merge(toptable_SvsCinHighAbs_J, New_IDs_HighAbs_J, by = "ID", all.x = TRUE)
merged_IDs_HighAbs_J$CurrentID <- ifelse(is.na(merged_IDs_HighAbs_J$CurrentID), merged_IDs_HighAbs_J$ID, merged_IDs_HighAbs_J$CurrentID)
merged_IDs_HighAbs_J <- merged_IDs_HighAbs_J[merged_IDs_HighAbs_J$CurrentID != 0,]

#Remove ID and Symbol coloumns, as we will use the new ID and call the new Symbols
merged_IDs_HighAbs_J$ID <- NULL
merged_IDs_HighAbs_J$SYMBOL <- NULL

#Change the header (CurrentID) to ID
colnames(merged_IDs_HighAbs_J)[7] <- "ID"

# Reannotated the EntrezGeneIDs
ids_J <- bitr(merged_IDs_HighAbs_J$ID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")

# Combine toptable and gene symbol
toptable_SvsCinHighAbs_J <- merge(merged_IDs_HighAbs_J, ids_J, by.x="ID", by.y="ENTREZID", all=T)

#4-Low absorbers Jejunum
#For genes that are replaced new ID (Current ID) replace that in the toptable and genes that were removed from NCBI have a 0, so remove all the zeros:
New_IDs_LowAbs_J <- read.delim("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NA/gene_result_LowAbsorbers_J.txt", header = TRUE, sep = "\t", dec = ".")

# Replace NA IDs with new IDs and remove rows with new ID of 0
merged_IDs_LowAbs_J <- merge(toptable_SvsCinLowAbs_J, New_IDs_LowAbs_J, by = "ID", all.x = TRUE)
merged_IDs_LowAbs_J$CurrentID <- ifelse(is.na(merged_IDs_LowAbs_J$CurrentID), merged_IDs_LowAbs_J$ID, merged_IDs_LowAbs_J$CurrentID)
merged_IDs_LowAbs_J <- merged_IDs_LowAbs_J[merged_IDs_LowAbs_J$CurrentID != 0,]

#Remove ID and Symbol coloumns, as we will use the new ID and call the new Symbols
merged_IDs_LowAbs_J$ID <- NULL
merged_IDs_LowAbs_J$SYMBOL <- NULL

#Change the header (CurrentID) to ID
colnames(merged_IDs_LowAbs_J)[7] <- "ID"

# Reannotated the EntrezGeneIDs
ids_J <- bitr(merged_IDs_LowAbs_J$ID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")

# Combine toptable and gene symbol
toptable_SvsCinLowAbs_J <- merge(merged_IDs_LowAbs_J, ids_J, by.x="ID", by.y="ENTREZID", all=T)

#Write toptable for Duodenum
write.table(toptable_SvsCinHighAbs_D, file = "output/NewOutput/toptable_duodenum_multi-level_Campesterol_chol_Stanol_VS_Control_HighAbsorber.csv", sep=",", row.names=FALSE)
write.table(toptable_SvsCinLowAbs_D, file = "output/NewOutput/toptable_duodenum_multi-level_Campesterol_chol_Stanol_VS_Control_LowAbsorber.csv", sep=",", row.names=FALSE)

#Write toptable for Jejunum
write.table(toptable_SvsCinHighAbs_J, file = "output/NewOutput/toptable_jejunum_multi-level_Campesterol_chol_Stanol_VS_Control_HighAbsorber.csv", sep=",", row.names=FALSE)
write.table(toptable_SvsCinLowAbs_J, file = "output/NewOutput/toptable_jejunum_multi-level_Campesterol_chol_Stanol_VS_Control_LowAbsorber.csv", sep=",", row.names=FALSE)

########################################################
#####Venn diagrams on Campesterol/cholesterol ratio#####
########################################################

#Determine DEG in Duodenum
DEG_SvsCinHighAbs_D <- toptable_SvsCinHighAbs_D[toptable_SvsCinHighAbs_D$P.Value< 0.05, c(1,2,5,8)]
DEG_SvsCinLowAbs_D <- toptable_SvsCinLowAbs_D[toptable_SvsCinLowAbs_D$P.Value< 0.05, c(1,2,5,8)]
# DEG_Diff_D <- toptable_Diff_D[toptable_Diff_D$P.Value< 0.05, c(1,2,5,8)] #E3E418FF
###1:ID,2:logFC,5:P.value,8:SYMBOL.

# specify if a gene is up or down regulated
DEG_SvsCinHighAbs_D$class <- ifelse(DEG_SvsCinHighAbs_D$logFC>0,"up","down")
DEG_SvsCinLowAbs_D$class <- ifelse(DEG_SvsCinLowAbs_D$logFC>0,"up","down")

write.table(DEG_SvsCinHighAbs_D, file = "output/NewOutput/DEG_SvsCinHighAbs_D_camp_chol.csv", sep=",", row.names=FALSE)
write.table(DEG_SvsCinLowAbs_D, file = "output/NewOutput/DEG_SvsCinLowAbs_D_camp_chol.csv", sep=",", row.names=FALSE)

#Determine DEG in Jejunum
DEG_SvsCinHighAbs_J <- toptable_SvsCinHighAbs_J[toptable_SvsCinHighAbs_J$P.Value< 0.05, c(1,2,5,8)]
DEG_SvsCinLowAbs_J <- toptable_SvsCinLowAbs_J[toptable_SvsCinLowAbs_J$P.Value< 0.05, c(1,2,5,8)]
# DEG_Diff_J <- toptable_Diff_J[toptable_Diff_J$P.Value< 0.05, c(1,2,5,8)] #E3E418FF

# specify if a gene is up or down regulated
DEG_SvsCinHighAbs_J$class <- ifelse(DEG_SvsCinHighAbs_J$logFC>0,"up","down")
DEG_SvsCinLowAbs_J$class <- ifelse(DEG_SvsCinLowAbs_J$logFC>0,"up","down")

write.table(DEG_SvsCinHighAbs_J, file = "output/NewOutput/DEG_SvsCinHighAbs_J_camp_chol.csv", sep=",", row.names=FALSE)
write.table(DEG_SvsCinLowAbs_J, file = "output/NewOutput/DEG_SvsCinLowAbs_J_camp_chol.csv", sep=",", row.names=FALSE)
###########################
#Venn diagram of DEG 
##We would like to have the font type in "calibri"
#Find the path to all fonts
font_paths()
font.files()
font_files() %>% tibble()
font_files() %>% tibble() %>% filter(str_detect(family, "Calibri"))

#Venn diagram of DEG in Duodenum
venn.diagram(x = list(DEG_SvsCinHighAbs_D$ID, DEG_SvsCinLowAbs_D$ID),
             category.names = c("High-cholesterol absorbers","Low-cholesterol absorbers"),
             cat.cex = 1.2,
             fontfamily = "sans",
             cex = 1.5,
             cat.fontfamily ="sans",
             filename = 'output/NewOutput/venn_diagram_duodenum_camp_chol.png',
             output=FALSE,
             col=c("#472D7BFF","#1F9A8AFF"),
             fill = c(alpha("#472D7BFF",0.3),alpha("#1F9A8AFF",0.3)),
             cat.pos = 4,
             main = "Duodenum: plant stanol versus control",
             main.cex = 1.5,
             main.pos = c (0.39, 1.05),
             main.fontfamily = "sans",
)

#Venn diagram of DEG in Jejunum
venn.diagram(x = list(DEG_SvsCinHighAbs_J$ID, DEG_SvsCinLowAbs_J$ID),
             category.names = c("High-cholesterol absorbers","Low-cholesterol absorbers"),
             cat.cex = 1.2,
             fontfamily = "sans",
             cex = 1.5,
             cat.fontfamily ="sans",
             filename = 'output/NewOutput/venn_diagram_jejunum_camp_chol.png',
             output=FALSE,
             col=c("#472D7BFF","#1F9A8AFF"),
             fill = c(alpha("#472D7BFF",0.3),alpha("#1F9A8AFF",0.3)),
             cat.pos = 4,
             main =  "Jejunum: plant stanol versus control",
             main.cex = 1.5,
             main.pos = c (0.39, 1.05),
             main.fontfamily = "sans",
)

#Groups genes based on Venn diagram in Duodenum
VD_D <- calculate.overlap(x = list(DEG_SvsCinHighAbs_D$SYMBOL, DEG_SvsCinLowAbs_D$SYMBOL))
D_a1 <- as.data.frame(VD_D$a1)
D_a2 <- as.data.frame(VD_D$a2)
D_a3 <- as.data.frame(VD_D$a3)
# D_a4 <- as.data.frame(VD_D$a4)
# D_a5 <- as.data.frame(VD_D$a5)
# D_a6 <- as.data.frame(VD_D$a6)
# D_a7 <- as.data.frame(VD_D$a7)

#Groups genes based on Venn diagram in Jejunum
VD_J <- calculate.overlap(x = list(DEG_SvsCinHighAbs_J$SYMBOL, DEG_SvsCinLowAbs_J$SYMBOL))
J_a1 <- as.data.frame(VD_J$a1)
J_a2 <- as.data.frame(VD_J$a2)
J_a3 <- as.data.frame(VD_J$a3)
# J_a4 <- as.data.frame(VD_J$a4)
# J_a5 <- as.data.frame(VD_J$a5)
# J_a6 <- as.data.frame(VD_J$a6)
# J_a7 <- as.data.frame(VD_J$a7)

#Combine list of genes from Venn diagram in Duodenun
VD_D_genes <- rbind.fill(D_a1, D_a2, D_a3)
#Combine list of genes from Venn diagram in Jejunum
VD_J_genes <- rbind.fill(J_a1, J_a2, J_a3)

#Write Venn diagram genes for Duodenum
write.table(VD_D_genes, file = "output/NewOutput/Venn diagram genelist_duodenum_camp_chol.csv", sep=",", row.names=FALSE)
#Write Venn diagram genes for Jejunum
write.table(VD_J_genes, file = "output/NewOutput/Venn diagram genelist_jejunum_camp_chol.csv", sep=",", row.names=FALSE)

########################################################
#####Volcano plots on Campesterol/cholesterol ratio#####
########################################################

#Volcano plot of DEG in SvsCinHighAbs at Duodenum
png('output/NewOutput/volcanoplot_duodenum_camp_chol_Stanol_VS_Control_HighAbs.png')
EnhancedVolcano(toptable_SvsCinHighAbs_D, title = "Duodenum: stanol versus control/high-cholesterol absorbers", lab = toptable_SvsCinHighAbs_D$SYMBOL, titleLabSize = 14,
                labSize = 3, x = 'logFC', xlim = c(-2.5,2.5), y = 'P.Value', ylim = c(0,4), pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

#Volcano plot of DEG in SvsCinLowAbs at Duodenum
png('output/NewOutput/volcanoplot_duodenum_camp_chol_Stanol_VS_Control_LowAbs.png')
EnhancedVolcano(toptable_SvsCinLowAbs_D, title = "Duodenum: stanol versus control/low-cholesterol absorbers", lab = toptable_SvsCinLowAbs_D$SYMBOL, titleLabSize = 14,
                labSize = 3, x = 'logFC', xlim = c(-2.5,2.5), y = 'P.Value', ylim = c(0,4), pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

#Volcano plot of DEG in SvsCinHighAbs at Jejunum
png('output/NewOutput/volcanoplot_jejunum_camp_chol_Stanol_VS_Control_HighAbsorber.png')
EnhancedVolcano(toptable_SvsCinHighAbs_J, title = "Jejunum: stanol versus control/high-cholesterol absorbers", lab = toptable_SvsCinHighAbs_J$SYMBOL, titleLabSize = 14,
                labSize = 3, x = 'logFC', xlim = c(-2.5,2.5), y = 'P.Value', ylim = c(0,4), pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

#Volcano plot of DEG in SvsCinLowAbs at Jejunum
png('output/NewOutput/volcanoplot_jejunum_camp_chol_Stanol_VS_Control_LowAbsorber.png')
EnhancedVolcano(toptable_SvsCinLowAbs_J, title = "Jejunum: stanol versus control/low-cholesterol absorbers", lab = toptable_SvsCinLowAbs_J$SYMBOL, titleLabSize = 14,
                labSize = 3, x = 'logFC', xlim = c(-2.5,2.5), y = 'P.Value', ylim = c(0,4), pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

