########################################################
######Enrichment analysis pathway analysis#######
########################################################

#clean workspace TEST
rm(list=ls())

library(ggplot2)
library(magrittr)

#set working directory
setwd("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC")

##Duodenum##

#1# High_Absorber_Duodenum

###KEGG###
#############################################################
#1# KEGG_Camp_chol_duodenum_Stanol_vs_Control_in_HighAbsorber

##Read table
Pathway_KEGG_HighAbs_d <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_KEGG_Camp_chol_duodenum_Stanol_vs_Control_in_HighAbsorber", header = TRUE, sep = "\t", dec = ".")
   
##Select the first 15 pathways#
Pathway_KEGG_HighAbs_d_15Top <- head(Pathway_KEGG_HighAbs_d, 15)

#Make a ggplot (barplot) of the top 15 pathways and save as PDF
ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/KEGG_High_d.pdf", plot = ggplot(Pathway_KEGG_HighAbs_d_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
        geom_bar(stat='identity', position=position_dodge()) +
         labs(x='KEGG Pathway', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(3))+
         coord_flip() +
         theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
               axis.title = element_text(size = 8)))  # Increase font size of axis titles) 

###Reactome###
###################################################
#2# Reactome_Camp_chol_duodenum_Stanol_vs_Control_in_HighAbsorber

##Read table
Pathway_Reactome_HighAbs_d <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_Reactome_Camp_chol_duodenum_Stanol_vs_Control_in_HighAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_Reactome_HighAbs_d_15Top <- head(Pathway_Reactome_HighAbs_d, 15)

#Make a ggplot (barplot) of the top 15 pathways and save as PDF
ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/Reactome_High_d.pdf", plot = ggplot(Pathway_Reactome_HighAbs_d_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
         geom_bar(stat='identity', position=position_dodge()) +
         labs(x='Reactome Pathway', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(5))+
         coord_flip() +
         theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
               axis.title = element_text(size = 8)))  # Increase font size of axis titles) 

###Wikipathways###
###################################################
#3# Wikipathways_Camp_chol_duodenum_Stanol_vs_Control_in_HighAbsorber

##Read table
Pathway_Wikipathways_HighAbs_d <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_enrichResults_Camp_chol_duodenum_Stanol_vs_Control_in_HighAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_Wikipathways_HighAbs_d_15Top <- head(Pathway_Wikipathways_HighAbs_d, 16)

#Change name of WP5333
Pathway_Wikipathways_HighAbs_d_15Top$Description[10] <- "Enterocyte cholesterol metabolism"

##Do delete WP5122##
# Remove the row with ID = 3
Pathway_Wikipathways_HighAbs_d_15Top <- Pathway_Wikipathways_HighAbs_d_15Top[-3, ]


#Make a ggplot (barplot) of the top 15 pathways and save as PDF
ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/Wikipathways_High_d.pdf", plot = ggplot(Pathway_Wikipathways_HighAbs_d_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
         geom_bar(stat='identity', position=position_dodge()) +
         labs(x='WikiPathways', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(4))+
         coord_flip() + 
       theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
             axis.title = element_text(size = 8)))  # Increase font size of axis titles) 


##########################
#2# Low_Absorber_Duodenum

###KEGG###
#############################################################
#1# KEGG_Camp_chol_duodenum_Stanol_vs_Control_in_LowAbsorber

##Read table
Pathway_KEGG_LowAbs_d <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_KEGG_Camp_chol_duodenum_Stanol_vs_Control_in_LowAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_KEGG_LowAbs_d_15Top <- head(Pathway_KEGG_LowAbs_d, 15)

#Make a ggplot (barplot) of the top 15 pathways and save as PDF
ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/KEGG_Low_d.pdf", plot = ggplot(Pathway_KEGG_LowAbs_d_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
         geom_bar(stat='identity', position=position_dodge()) +
         labs(x='KEGG Pathway', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(7))+
         coord_flip() +
         theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
               axis.title = element_text(size = 8)))  # Increase font size of axis titles)

###Reactome###
###################################################
#2# Reactome_Camp_chol_duodenum_Stanol_vs_Control_in_LowAbsorber

##Read table
Pathway_Reactome_LowAbs_d <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_Reactome_Camp_chol_duodenum_Stanol_vs_Control_in_LowAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_Reactome_LowAbs_d_15Top <- head(Pathway_Reactome_LowAbs_d, 15)

#Make a ggplot (barplot) of the top 15 pathways, and save it as a pdf.

ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/Reactome_Low_d.pdf", plot = ggplot(Pathway_Reactome_LowAbs_d_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
         geom_bar(stat='identity', position=position_dodge()) +
         labs(x='Reactome Pathway', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(10))+
         coord_flip() +
         theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
               axis.title = element_text(size = 8)))  # Increase font size of axis titles)


###Wikipathways###
###################################################
#3# Wikipathways_Camp_chol_duodenum_Stanol_vs_Control_in_LowAbsorber

##Read table
Pathway_Wikipathways_LowAbs_d <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_enrichResults_Camp_chol_duodenum_Stanol_vs_Control_in_LowAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_Wikipathways_LowAbs_d_15Top <- head(Pathway_Wikipathways_LowAbs_d, 15)

#Make a ggplot (barplot) of the top 15 pathways, and save it as a pdf.

ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/Wikipathways_Low_d.pdf", plot = ggplot(Pathway_Wikipathways_LowAbs_d_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
         geom_bar(stat='identity', position=position_dodge()) +
         labs(x='WikiPathways', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(5))+
         coord_flip() +
theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
      axis.title = element_text(size = 8)))  # Increase font size of axis titles) 





##Jejunum##

#1# High_Absorber_Jejunum

###KEGG###
#############################################################
#1# KEGG_Camp_chol_jejunum_Stanol_vs_Control_in_HighAbsorber

##Read table
Pathway_KEGG_HighAbs_J <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_KEGG_Camp_chol_jejunum_Stanol_vs_Control_in_HighAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_KEGG_HighAbs_J_15Top <- head(Pathway_KEGG_HighAbs_J, 15)

#Make a ggplot (barplot) of the top 15 pathways, and save it as a pdf.
ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/KEGG_high_J.pdf", plot = ggplot(Pathway_KEGG_HighAbs_J_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
         geom_bar(stat='identity', position=position_dodge()) +
         labs(x='KEGG Pathway', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(2))+
         coord_flip() +
         theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
               axis.title = element_text(size = 8)))  # Increase font size of axis titles) 

###Reactome###
###################################################
#2# Reactome_Camp_chol_jejunum_Stanol_vs_Control_in_HighAbsorber

##Read table
Pathway_Reactome_HighAbs_J <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_Reactome_Camp_chol_jejunum_Stanol_vs_Control_in_HighAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_Reactome_HighAbs_J_15Top <- head(Pathway_Reactome_HighAbs_J, 15)

#Make a ggplot (barplot) of the top 15 pathways, and save it as a pdf.
ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/Reactome_high_J.pdf", plot = ggplot(Pathway_Reactome_HighAbs_J_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
         geom_bar(stat='identity', position=position_dodge()) +
         labs(x='Reactome Pathway', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(5))+
         coord_flip() +
         theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
               axis.title = element_text(size = 8)))  # Increase font size of axis titles) 

###Wikipathways###
###################################################
#3# Wikipathways_Camp_chol_jejunum_Stanol_vs_Control_in_HighAbsorber

##Read table
Pathway_Wikipathways_HighAbs_J <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_enrichResults_Camp_chol_jejunum_Stanol_vs_Control_in_HighAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_Wikipathways_HighAbs_J_15Top <- head(Pathway_Wikipathways_HighAbs_J, 15)

#Make a ggplot (barplot) of the top 15 pathways, and save it as a pdf.
ggsave("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/output/NewOutput/Enrichment/p1q1/Wikipathways_high_J.pdf", plot = ggplot(Pathway_Wikipathways_HighAbs_J_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
         geom_bar(stat='identity', position=position_dodge()) +
         labs(x='WikiPathways', y='Gene number') +
         scale_y_continuous(breaks = scales::pretty_breaks(2))+
         coord_flip() +
         theme(axis.text = element_text(size = 8),  # Increase font size of axis labels
               axis.title = element_text(size = 8)))  # Increase font size of axis titles)

##########################
#2# Low_Absorber_Jejunum

###KEGG###
#############################################################
#1# KEGG_Camp_chol_jejunum_Stanol_vs_Control_in_LowAbsorber

##Read table
Pathway_KEGG_LowAbs_J <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_KEGG_Camp_chol_jejunum_Stanol_vs_Control_in_LowAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_KEGG_LowAbs_J_15Top <- head(Pathway_KEGG_LowAbs_J, 15)

#Make a ggplot (barplot) of the top 15 pathways#
ggplot(Pathway_KEGG_LowAbs_J_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
  geom_bar(stat='identity', position=position_dodge()) +
  labs(x='KEGG Pathway', y='Gene number') +
  scale_y_continuous(breaks = scales::pretty_breaks(8))+
  coord_flip()

###Reactome###
###################################################
#2# Reactome_Camp_chol_Jejunum_Stanol_vs_Control_in_LowAbsorber

##Read table
Pathway_Reactome_LowAbs_J <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_Reactome_Camp_chol_jejunum_Stanol_vs_Control_in_LowAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_Reactome_LowAbs_J_15Top <- head(Pathway_Reactome_LowAbs_J, 15)

#Make a ggplot (barplot) of the top 15 pathways#
ggplot(Pathway_Reactome_LowAbs_J_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
  geom_bar(stat='identity', position=position_dodge()) +
  labs(x='Reactome Pathway', y='Gene number') +
  scale_y_continuous(breaks = scales::pretty_breaks(12))+
  coord_flip()

###Wikipathways###
###################################################
#3# Wikipathways_Camp_chol_Jejunum_Stanol_vs_Control_in_LowAbsorber

##Read table
Pathway_Wikipathways_LowAbs_J <- read.delim("output/NewOutput/Enrichment/p1q1/p=1_q=1_enrichResults_Camp_chol_jejunum_Stanol_vs_Control_in_LowAbsorber", header = TRUE, sep = "\t", dec = ".")

##Select the first 15 pathways#
Pathway_Wikipathways_LowAbs_J_15Top <- head(Pathway_Wikipathways_LowAbs_J, 15)

#Make a ggplot (barplot) of the top 15 pathways#
ggplot(Pathway_Wikipathways_LowAbs_J_15Top, aes(x=reorder(Description,Count), y=Count, fill=qvalue)) +
  geom_bar(stat='identity', position=position_dodge()) +
  labs(x='WikiPathways', y='Gene number') +
  scale_y_continuous(breaks = scales::pretty_breaks(9))+
  coord_flip()






