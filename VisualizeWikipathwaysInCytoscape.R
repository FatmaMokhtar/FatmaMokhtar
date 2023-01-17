######################################
#Visualize Wikipathways with Cytoscape#
######################################

#clean workspace
rm(list=ls())

#set working directory
setwd("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/input/Wikipathway/VisualizeWikiCyto")

if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways")
}
#In order to create networks in Cytoscape from R you need RCy3#
library(rWikiPathways)
if(!"RCy3" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("RCy3")
}
#install.packages("brewer.pal")

library(RCy3)
library(brewer.pal)

##Install WikipathwayApp##
installApp("Wikipathways")

###############################################################################
###############################################################################
#1- Use one toptable of the Low absorbers in the Duodenum and view in Cytscape#
###############################################################################
###############################################################################

###########Prepare Tables#################
#####READ the needed Tables#######

##Read toptable##
toptable_SvsCinLowAbs_D <- read.table("toptable_duodenum_multi-level_Campesterol_chol_Stanol_VS_Control_LowAbsorber.csv", header=TRUE,sep=",")
toptable_SvsCinLowAbs_D <- toptable_SvsCinLowAbs_D[,c("ID","logFC","P.Value","SYMBOL")]


###########Open Cytoscape and load pathway of interest and the toptable data#################

#The whole point of RCy3 is to connect with Cytoscape, Make sure that Cytoscape is open#
cytoscapePing()

#Confirm that Cytoscape is installed and opened#
cytoscapeVersionInfo ()

#Import wikipathway of interest ## T-cell activation SARS-CoV-2##WP5098##
wp.cmd = 'wikipathways import-as-pathway id="WP5098"'
commandsGET(wp.cmd)

#Import the toptable into Cytoscape#
#load data to the imported pathway in cytoscape by key column as SYMBOL#
loadTableData(table = "node", data = toptable_SvsCinLowAbs_D, data.key.column = "SYMBOL", table.key.column = "name")

###########Visual style#################
RCy3::copyVisualStyle("default","ppi")#Create a new visual style (ppi) by copying a specified style (default)
RCy3::setNodeLabelMapping("name", style.name="ppi")
RCy3::lockNodeDimensions(TRUE, style.name="ppi")#Set a boolean value to have node width and height fixed to a single size value.

#Adjust the dimensions of the nodes#
RCy3::setNodeWidthDefault(100, style.name = NULL)
RCy3::setNodeHeightDefault(25, style.name = NULL)
RCy3::setNodeBorderWidthDefault(0.001, style.name = NULL)
RCy3::setNodeBorderColorDefault('#89D0F5')#

#threshold is set based of differential expressed gene criteria
data.values<-c(-1,0,1)
#red-blue color schema chosen
node.colors <- c(brewer.pal(length(data.values), "RdBu"))
#nodes are split to show both log2fc values for both diseases
RCy3::setNodeCustomHeatMapChart(c("logFC"), slot = 1, style.name = "ppi", colors = c("#CC3300","#FFFFFF","#6699FF","#CCCCCC"))
RCy3::setNodeColorMapping("logFC", colors=paletteColorBrewerRdBu,style.name = "ppi")
RCy3::setVisualStyle("ppi")

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #Set node border width and color based on p-value
  #First we need to get all p-values from node table
  pvalues <- getTableColumns(table = 'node', columns = 'P.Value')
pvalues <- na.omit(pvalues)
#Create a range for all sign. p-values, and one for all not significant.
significant_pvalues <- pvalues[(pvalues < 0.05)]
not.significant_pvalues <- pvalues[(pvalues >= 0.05)]
significant_pvalues.colors <- rep("#2e9d1d", length(significant_pvalues))
not.significant_pvalues.colors <- rep("#FFFFFF", length(not.significant_pvalues))
RCy3::setNodeBorderWidthMapping('P.Value', table.column.values = NULL , c(6,6) , mapping.type = "c", style.name = "ppi")
RCy3::setNodeBorderColorMapping('P.Value', c(significant_pvalues,not.significant_pvalues), c(significant_pvalues.colors, not.significant_pvalues.colors), default.color = "#AAAAAA", mapping.type = "d", style.name = "ppi")
##Update relevant interactions to directional ones:
#Save output
setwd("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/input/Wikipathway/VisualizeWikiCyto")
exportImage('Pathwayvisualisation_Duodenum_LowAbsorbers_LogFC_2_T-cellactivationSARS-CoV-2_WP5098', 'PDF')

###############################################################################
###############################################################################
#2- Use one toptable of the High absorbers in the Duodenum and view in Cytscape#
###############################################################################
###############################################################################

###########Prepare Tables#################
#####READ the needed Tables#######

##Read toptable##

toptable_SvsCinHighAbs_D <- read.table("toptable_duodenum_multi-level_Campesterol_chol_Stanol_VS_Control_HighAbsorber.csv", header=TRUE,sep=",")
toptable_SvsCinHighAbs_D <- toptable_SvsCinHighAbs_D[,c("ID","logFC","P.Value","SYMBOL")]

###########Open Cytoscape and load pathway of interest and the toptable data#################

#The whole point of RCy3 is to connect with Cytoscape, Make sure that Cytoscape is open#
cytoscapePing()

#Confirm that Cytoscape is installed and opened#
cytoscapeVersionInfo ()

#Import wikipathway of interest ## T-cell activation SARS-CoV-2##WP5098##
wp.cmd = 'wikipathways import-as-pathway id="WP5098"'
commandsGET(wp.cmd)

#Import the toptable into Cytoscape#
#load data to the imported pathway in cytoscape by key column as SYMBOL#
loadTableData(table = "node", data = toptable_SvsCinHighAbs_D, data.key.column = "SYMBOL", table.key.column = "name")

###########Visual style#################
RCy3::copyVisualStyle("default","ppi")#Create a new visual style (ppi) by copying a specified style (default)
RCy3::setNodeLabelMapping("name", style.name="ppi")
RCy3::lockNodeDimensions(TRUE, style.name="ppi")#Set a boolean value to have node width and height fixed to a single size value.

#Adjust the dimensions of the nodes#
RCy3::setNodeWidthDefault(100, style.name = NULL)
RCy3::setNodeHeightDefault(25, style.name = NULL)
RCy3::setNodeBorderWidthDefault(0.001, style.name = NULL)
RCy3::setNodeBorderColorDefault('#89D0F5')#

#threshold is set based of differential expressed gene criteria
data.values<-c(-1,0,1)
#red-blue color schema chosen
node.colors <- c(brewer.pal(length(data.values), "RdBu"))
#nodes are split to show both log2fc values for both diseases
RCy3::setNodeCustomHeatMapChart(c("logFC"), slot = 1, style.name = "ppi", colors = c("#CC3300","#FFFFFF","#6699FF","#CCCCCC"))
RCy3::setNodeColorMapping("logFC", colors=paletteColorBrewerRdBu,style.name = "ppi")
RCy3::setVisualStyle("ppi")

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Set node border width and color based on p-value
#First we need to get all p-values from node table
pvalues <- getTableColumns(table = 'node', columns = 'P.Value')
pvalues <- na.omit(pvalues)
#Create a range for all sign. p-values, and one for all not significant.
significant_pvalues <- pvalues[(pvalues < 0.05)]
not.significant_pvalues <- pvalues[(pvalues >= 0.05)]
significant_pvalues.colors <- rep("#2e9d1d", length(significant_pvalues))
not.significant_pvalues.colors <- rep("#FFFFFF", length(not.significant_pvalues))
RCy3::setNodeBorderWidthMapping('P.Value', table.column.values = NULL , c(6,6) , mapping.type = "c", style.name = "ppi")
RCy3::setNodeBorderColorMapping('P.Value', c(significant_pvalues,not.significant_pvalues), c(significant_pvalues.colors, not.significant_pvalues.colors), default.color = "#AAAAAA", mapping.type = "d", style.name = "ppi")
##Update relevant interactions to directional ones:
#Save output
setwd("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/input/Wikipathway/VisualizeWikiCyto")
exportImage('Pathwayvisualisation_Duodenum_Highbsorbers_LogFC_2_T-cellactivationSARS-CoV-2_WP5098', 'PDF')



############################################################################
############################################################################
#3-Merge the two toptables of the High and Low absorbers, to view in Cytscape#
############################################################################
############################################################################

###########Prepare Tables#################

#READ needed Tables#
toptable_SvsCinLowAbs_D <- read.table("toptable_duodenum_multi-level_Campesterol_chol_Stanol_VS_Control_LowAbsorber.csv", header=TRUE,sep=",")
toptable_SvsCinHighAbs_D <- read.table("toptable_duodenum_multi-level_Campesterol_chol_Stanol_VS_Control_HighAbsorber.csv", header=TRUE,sep=",")

#Keep the needed columns#
toptable_SvsCinLowAbs_D <- toptable_SvsCinLowAbs_D[,c("ID","logFC","P.Value","SYMBOL")]
toptable_SvsCinHighAbs_D <- toptable_SvsCinHighAbs_D[,c("ID","logFC","P.Value","SYMBOL")]

#Add High or low to the headers of the tables#
colnames(toptable_SvsCinLowAbs_D)[-4]  <- paste0(colnames(toptable_SvsCinLowAbs_D)[-4], "_Low")
colnames(toptable_SvsCinHighAbs_D)[-4] <- paste0(colnames(toptable_SvsCinHighAbs_D)[-4], "_High")

#Merge two data frames by SYMBOL#
toptable_SvsC_D_Merged_High_Low <- merge(toptable_SvsCinLowAbs_D, toptable_SvsCinHighAbs_D, by.x = "SYMBOL", by.y = "SYMBOL" )

###########Open Cytoscape and load pathway of interest and the merged data#################

#The whole point of RCy3 is to connect with Cytoscape, Make sure that Cytoscape is open#
cytoscapePing()

#Confirm that Cytoscape is installed and opened#
cytoscapeVersionInfo ()

#Import wikipathway of interest ## T-cell activation SARS-CoV-2##WP5098##
wp.cmd = 'wikipathways import-as-pathway id="WP5098"'
commandsGET(wp.cmd)

#Import the merged toptable into Cytoscape#
#load data to the imported pathway in cytoscape by key column as SYMBOL#
loadTableData(table = "node", data = toptable_SvsC_D_Merged_High_Low, data.key.column = "SYMBOL", table.key.column = "name")

###########Visual style#################
RCy3::copyVisualStyle("default","ppi")#Create a new visual style (ppi) by copying a specified style (default)
RCy3::setNodeLabelMapping("name", style.name="ppi")
RCy3::lockNodeDimensions(TRUE, style.name="ppi")#Set a boolean value to have node width and height fixed to a single size value.


#threshold is set based of differential expressed gene criteria
data.values<-c(-1,0,1)
#red-blue color schema chosen
node.colors <- c(brewer.pal(length(data.values), "RdBu"))
#nodes are split to show both log2fc values for both groups
RCy3::setNodeCustomHeatMapChart(c("logFC_Low","logFC_High"), slot = 2, style.name = "ppi", colors = c("#CC3300","#FFFFFF","#6699FF","#CCCCCC"))
RCy3::setVisualStyle("ppi")

#Save pathway
setwd("/Users/fatma.mokhtar/Documents/ProjectWithSusan/New_Campesterol_StandardisedTC/input/Wikipathway/VisualizeWikiCyto")
exportImage('Pathwayvisualisation_Duodenum_Merged_High_Low_LogFC_2_T-cellactivationSARS-CoV-2_WP5098', 'PDF')
