#############################################################
#
# R script for Analysing Annotations data and generating summary data both across all genomes and summarized for phylo oenetic clusters:
# 
# Dombrowski et al., 2020
# 
# Finalized: February 2002
#
#############################################################
#rm(list=ls())
sessionInfo()

################################################
#0.1 setting working directory (!adjust wdir accordingly!)
################################################
wdir <- "C:/Users/phr001/OneDrive - University of Bergen/PHD/Methods/ANME_SRB/ANME_Annotation/Dombrowski_R"
setwd(wdir)

################################################
#0.2 load required packages 
################################################
library("ggplot2")
library("plyr") 
library("grid")
library("gplots")
library(tidyr)
library("gridExtra")
library("multcomp")
library("seqinr")
library("reshape2")
library("RColorBrewer")
library("dplyr")
library('tidyr')
library('tidyverse')
library(data.table)
library(ggtree)
library("kableExtra")

################################################################################################
################################################################################################
#1. load and clean tables
################################################################################################
################################################################################################

######################################################################
#1.1. mapping file that links binIDs to different taxonomy levels
######################################################################

#read in the mapping file
#Notice: We use the GroupID to order things later, this ID was manually added using the order of taxa in the archaeal phylogenies
design <- read.table("1_Input/mapping_binID_to_taxa.txt", sep="\t", header=T, fill=TRUE, quote = "")
kable(head(design), format='markdown')

# grouping
design <- read.table("1_Input/mapping_binID_to_taxa_grouping.txt", sep="\t", header=T, fill=TRUE, quote = "")
kable(head(design), format='markdown')

#add a new column that links the BinID and the taxonomic cluster
design$NewName2 <-paste(design$BinID, design$Cluster, sep = "-")
kable(head(design), format='markdown')

#transform mapping file and summarize how many genomes we have/cluster
Number_of_taxa <- ddply(design, .(Cluster), summarize, NrGenomes = length(Cluster))

#order our data
Number_of_taxa <- Number_of_taxa[order(Number_of_taxa$NrGenomes, decreasing = TRUE),] 

#view data
kable((Number_of_taxa), format='markdown')

#print the summary of nr. of genomes per cluster for the records
write.table(Number_of_taxa, "2_Output/Number_of_clusters.txt",  sep = "\t", quote = F, row.names = F, na = "")

#add a new column with the cluster count info into mapping file (then it can be added as a label into figures)
design_2 <- merge(design, Number_of_taxa, by = "Cluster" )

#view whether all went fine
kable(head(design_2), format='markdown')

#generate a new column into the mapping file that links the cluster name with the number of genomes of clusters in each cluster
design_2$ClusterName <- paste(design_2$Cluster, " (", design_2$NrGenomes, ")", sep = "")

#view data
kable(head(design_2), format='markdown')

#generate some vectors to order the data frame if needed
#Notice at this moment the order is defined as by the order of the original mapping file. The mapping file is order by the group ID in the mapping file (see line 84), but this can be changed as wished
#check the old data
design$Cluster

#check the new data
design_2$Cluster

#If we check the individual factors, we can see that the new dataframe is order by alphabet. Since this is not what we want, lets correct this
#reorder our old dataframe by the original one
design_2 <- design_2[ order(match(design_2$BinID, design$BinID)), ]

#check, if the basal clade is now in the correct order
design_2$Cluster

#Now, that the basal clade is in the first position, lets make vectors of the bins and clusters in order this we use:
#unique() = we use unique on all ClusterNames in our mapping file. That way instead of having repeated names, we only have the unique ones
#as.character() = we make sure that our R object we generated is a character (so has words).
#make a list to order our bins
Bin_order <- as.character(unique(design_2$BinID))
Bin_order

#make cluster order
Cluster_order <- as.character(unique(design_2$ClusterName))
Cluster_order

Bin_order2 <- as.character(unique(design_2$NewName2))
Bin_order2




######################################################################
#1.2a. read in mapping files to subset the annotations for pathways of interest (Table S12)
######################################################################

###
#General comment:
#These files list either the arcogs or KOs ordered based on pathways of interest.
#For example, all genes of gluconeogenesis are listed in order.
#This was used to generate Table S12
###

#general mapping file for arcog IDs
Arcog_mapping <- read.table("1_Input/ar14_arCOGdef19.txt", sep="\t", header=T, fill=TRUE, quote = "")
kable((head(Arcog_mapping)), format='markdown')

#pathway mapping file
Metabolism_file_KEGG <- read.table("1_Input/Metabolism_Table_KO_Apr2020.txt", sep="\t", header=T, fill=TRUE, quote = "")
kable((head(Metabolism_file_KEGG)), format='markdown')

#1.2b. read in genes of interest that we want to plot in a heatmap
################################################
#list of metabolic genes
#load the genes of interest
Genes_of_interest <- read.table("1_Input/Genes_of_interest_metabolism_LithoHsensingCfix.txt", sep="\t", header=T, fill=TRUE, quote = "")
kable((head(Genes_of_interest)), format='markdown')

# all hydrogenases
Genes_of_interest <- read.table("1_Input/Genes_of_interest_metabolism_hydrogenases.txt", sep="\t", header=T, fill=TRUE, quote = "")
kable((head(Genes_of_interest)), format='markdown')

#AOM / Vulcano et al, 2022
Genes_of_interest <- read.table("1_Input/Genes_of_interest_metabolism_Vulcano_etal2022_1.txt", sep="\t", header=T, fill=TRUE, quote = "")
kable((head(Genes_of_interest)), format='markdown')



#order metabolic genes
Genes_Metabolism_order_temp <- Genes_of_interest %>% arrange(Order, Order2)
#make a unique vector for our genes of interest
Genes_Metabolism_order <- as.character(unique(Genes_Metabolism_order_temp$Gene))
Genes_Metabolism_order

#define a order for metabolic pathways
Pathway_order <- as.character(unique(Genes_of_interest$pathway_2))
Pathway_order

#read in a mapping file for informational genes of interest

# Input of annotations file!
Input <- read.table("1_Input/ANME_AMOR_all_Annotations.txt", sep="\t", header=T, fill=TRUE, quote = "")
kable((head(Input)), format='markdown')

#generate Description table for all DBs of interest
Arcogs_Description <- unique(Input[,c("arcogs","arcogs_Description" )])
colnames(Arcogs_Description) <- c("Gene", "Description")
kable((head(Arcogs_Description)), format='markdown')

KOs_Description <- Input[,c("KO_hmm","Definition" )]
colnames(KOs_Description) <- c("Gene", "Description")
kable((head(KOs_Description)), format='markdown')

Pfam_Description <- unique(Input[,c("PFAM_hmm","PFAM_description" )])
colnames(Pfam_Description) <- c("Gene", "Description")
kable((head(Pfam_Description)), format='markdown')

TIGR_Description <- unique(Input[,c("TIRGR","TIGR_description" )])
colnames(TIGR_Description) <- c("Gene", "Description")
kable((head(TIGR_Description)), format='markdown')

Cazy_Description <- unique(Input[,c("CAZy","Description" )])
colnames(Cazy_Description) <- c("Gene", "Description")
kable((head(Cazy_Description)), format='markdown')

HydDB_Description <- unique(Input[,c("Description.1","Description.1" )])
colnames(HydDB_Description) <- c("Gene", "Description")
kable((head(HydDB_Description)), format='markdown')

#make a file with a description of all the ids for each search
All_Genes_Description <- rbind(Arcogs_Description,KOs_Description,Pfam_Description,TIGR_Description, Cazy_Description, HydDB_Description)

#print the column names to subset our datatable
colnames(Input)

#only keep the columns we actually want to work with
Input_subset = Input[,c('BinID','accession','arcogs','KO_hmm','PFAM_hmm','TIRGR','CAZy','Description.1' )]
kable((head(Input_subset)), format='markdown')

#convert dataframe from wide to long
Input_long <- reshape2::melt(Input_subset,  id=c("accession","BinID"))

#give informative headers
colnames(Input_long) <- c("accession", "BinID", "DB", "gene")

#clean factors, to remove issues when counting
Input_long$gene <- as.factor(Input_long$gene)
kable((head(Input_long)), format='markdown')

##Make count tables
#Generate a count table for our genomes of interest
#Now we want to count, how often does a genome (i.e. NIOZ134_mb_b41_2) have a gene. I.e. how often do we want arCOG00570, arCOG01358, …
#Do this via a loop (not executed, just an example)
#Notice:
#Since we run this chunk with , eval = FALSE we can still see the code but it is not executed. This is done because some computations take some time, which we do not want to spend, but I still want to show the code to give some alternative examples.
#count the number of proteins for each genome of interest
y <- c()
for (i in Bin_order) {
  x <-  table(subset(Input_long, BinID %in% paste(i))$gene)
  y <- cbind (y,x)
}

#clean-up the table
Counts_Table_loop <- y
colnames(Counts_Table_loop) <- Bin_order
Counts_Table_loop <- as.data.frame(Counts_Table_loop)
kable((head(Counts_Table_loop)), format='markdown')

#the '-' (=not identified genes) is also counted and listed in the first column and removed at this step
Counts_Table_loop <- Counts_Table_loop[-1,]
kable((head(Counts_Table_loop)), format='markdown')

#count data and clean header
Counts_Table_long <- ddply(Input_long, .(BinID, gene), summarize, GeneCount = length(gene))
colnames(Counts_Table_long) <- c("BinID", "geneID", "count")
kable((head(Counts_Table_long)), format='markdown')

#transform to wide format, with fill = 0 instead of a NA we add a 0
Counts_Table_wide <- spread(Counts_Table_long, BinID, count, fill = 0 )

#view data
kable((head(Counts_Table_wide)), format='markdown')

#count data and clean header
Counts_Table_long <- Input_long[,c('accession', 'BinID','gene')] %>% count(BinID, gene, sort = FALSE)
colnames(Counts_Table_long) <- c("BinID", "geneID", "count")
kable((head(Counts_Table_long)), format='markdown')

#delete rows with a minus symbol
Counts_Table_long <- Counts_Table_long[Counts_Table_long$geneID!= "-", ]

#clean factors
Counts_Table_long$geneID <- factor(Counts_Table_long$geneID)

#view data
kable((head(Counts_Table_long)), format='markdown')

#transform to wide format, with fill = 0 instead of a NA we add a 0
Counts_Table_wide <- spread(Counts_Table_long, BinID, count, fill = 0 )
kable((head(Counts_Table_wide)), format='markdown')

#change the rownames
rownames(Counts_Table_wide) <- Counts_Table_wide$geneID

#view data
kable((head(Counts_Table_wide)), format='markdown')

#delete the first column
Counts_Table_wide <- Counts_Table_wide[,-1]
kable((head(Counts_Table_wide)), format='markdown')

#order our data so that the bins start first with the bins from the basal group
Counts_Table_wide <- Counts_Table_wide[,Bin_order]

#view data
kable((head(Counts_Table_wide)), format='markdown')

##Generate a count table for our clusters of interest
#Same as above, but now we want to know for our 4 aquifer genomes, how many have Gene Xx and show this as percent. I.e. if 1/4 genomes have a gene, then 25% have it.
#First, lets merge in our taxa info into our count table, we need this to summarize our data by clusters.
#merge the count table with mapping file to add in the taxa info (might take a while depending on size)
Counts_Table_long_Tax <- merge(Counts_Table_long, design_2[,c("BinID", "ClusterName", "NrGenomes")], by = "BinID")
kable((head(Counts_Table_long_Tax)), format='markdown')

#convert counts to presence/absence matrix (just using 0/1) (this is needed to calculate the percentage across clusters)
Counts_Table_long_Tax$count[Counts_Table_long_Tax$count > 1] <- 1
kable((head(Counts_Table_long_Tax)), format='markdown')

#Now, we can use tidyr to count of how many genomes in a cluster have a gene; count data and clean header
Counts_Table_long_Tax_sum <- Counts_Table_long_Tax[,c('ClusterName', 'geneID','NrGenomes', 'count')] %>% count(ClusterName, geneID, NrGenomes, sort = FALSE)
colnames(Counts_Table_long_Tax_sum) <- c("ClusterName", "geneID", "NrGenomes", "quantity")
kable((head(Counts_Table_long_Tax_sum)), format='markdown')

#calculate of the percentage to answer of the total genomes per cluster how many have a certain gene
#notice: if running for your own data check here that your percentage makes sense. I.e. we do not want values above 100
Counts_Table_long_Tax_sum$percentage <- round(Counts_Table_long_Tax_sum$quantity/Counts_Table_long_Tax_sum$NrGenomes*100, digits = 0)
kable((head(Counts_Table_long_Tax_sum)), format='markdown')

#convert long to wide format and clean table (i.e. place the rownames)
Counts_Table_long_Tax_sum_wide <- spread(Counts_Table_long_Tax_sum[,c("geneID", "ClusterName", "percentage")], ClusterName, percentage)

#change the rownames
rownames(Counts_Table_long_Tax_sum_wide) <- Counts_Table_long_Tax_sum_wide$geneID

#view data
kable((head(Counts_Table_long_Tax_sum_wide)), format='markdown')

#delete the first column
Counts_Table_long_Tax_sum_wide <- Counts_Table_long_Tax_sum_wide[,-1]
kable((head(Counts_Table_long_Tax_sum_wide)), format='markdown')

#replace NAs with 0
Counts_Table_long_Tax_sum_wide[is.na(Counts_Table_long_Tax_sum_wide)] <- 0
kable((head(Counts_Table_long_Tax_sum_wide)), format='markdown')

#sort by cluster order (defined by the order of the mapping file)
Counts_Table_long_Tax_sum_wide <- Counts_Table_long_Tax_sum_wide[,Cluster_order]
kable((head(Counts_Table_long_Tax_sum_wide)), format='markdown')

#Merge our tables with the mapping data we have
#Now, that we have our count tables both for the bins as well as for all the clusters, we now want to add some gene description and subset the data based on different categores.
#For the bins
#Add gene descriptions
#Remember above, we made a list of descriptions that links all geneIDs with what is behind all the gene IDs? Now we want to add this info back in in order to print all the counts.
#merge
Counts_Table_final <- merge(All_Genes_Description, Counts_Table_wide, by.x="Gene", by.y="row.names", all.x = T, sort = F)
kable((head(Counts_Table_final)), format='markdown')

#print (and beautify elsewhere)
write.table(Counts_Table_final, "2_Output/Counts_Table_final.txt",  sep = "\t", quote = F, row.names = T, na = "")

#merge, Merging with the arcog_table
Arcog_Data <- merge(Arcog_mapping, Counts_Table_wide, by.x="arcog", by.y="row.names", all.x = T, sort = F)
kable((head(Arcog_Data)), format='markdown')

#print (and beautify elsewhere)
write.table(Arcog_Data, "2_Output/ArCOG_Data.txt",  sep = "\t", quote = F, row.names = T, na = "")

#merge, Merging with the metabolism metadata file
KEGG_Metabolism <- merge(Metabolism_file_KEGG, Counts_Table_wide, by.x="KO", by.y="row.names", all.x = T, sort = F)
kable((head(KEGG_Metabolism)), format='markdown')

#print
write.table(KEGG_Metabolism, "2_Output/KEGG_Metabolism.txt",  sep = "\t", quote = F, row.names = T, na = "")

##############################################################################################
#5. pull out genes of interest and make stats and plots --> METABOLISM genes of interest 
##############################################################################################
#subset genes of interest and clean factors
#Genes_metabolism <- subset(Genes_of_interest, pathway0 %in% "Metabolism")
#Genes_metabolism$Gene <- factor(Genes_metabolism$Gene)
#Genes_metabolism$GeneID <- factor(Genes_metabolism$GeneID)

#check how many genes we have
#dim(Genes_metabolism)

kable((head(Genes_of_interest)), format='markdown')
#kable((head(Genes_metabolism)), format='markdown')

#define an order (we arange the dataframe based on two columsn, Order and Order2)
Genes_metabolism_order_temp <- Genes_of_interest %>% arrange(Order, Order2)
Genes_metabolism_order <- as.character(unique(Genes_metabolism_order_temp$Gene))
length(Genes_metabolism_order)

Genes_metabolism_order
#The metabolism genes belong to different pathways, these pathways we want to show in two separate heatmaps
Metabolism_Pathway_order <- as.character(unique(Genes_of_interest$pathway_2))
Metabolism_Pathway_order

#Now that we know what genes we are interested in, lets subset our original count table.
#subset our original count table for genes of interest and clean factors
Genes_metabolism_counts <- subset(Counts_Table_long_Tax_sum, geneID %in% as.character(Genes_of_interest$GeneID))
Genes_metabolism_counts$geneID <- factor(Genes_metabolism_counts$geneID)

#control that all went fine
length(unique(Genes_metabolism_counts$geneID))

dim(Genes_metabolism_counts)

kable((head(Genes_metabolism_counts)), format='markdown')

#With length we see that now we just have 15 genes. How can we find out what gene is missing?
#setdiff() = we compare two vectors and print the elements that differ.

setdiff(Genes_of_interest$GeneID, Genes_metabolism_counts$geneID)

#Now, since our count data does not know that we categorize our different genes into different pathways, lets add this info in with merge()
#add in metadata (from the pathway info)
Key_metabolism_genes_cluster <- merge(Genes_metabolism_counts, Genes_of_interest, by.x ="geneID", by.y = 'GeneID' , all.x = T )
kable((head(Key_metabolism_genes_cluster)), format='markdown')

#define color code (not used for the current figure, but can be changed)
#here , we define 3 color levels, which sometimes is useful to show very clear cutoffs
Key_metabolism_genes_cluster$category <- ifelse(Key_metabolism_genes_cluster$percentage == 100, "1",
                                                ifelse(Key_metabolism_genes_cluster$percentage >= 75, "0.75",
                                                       ifelse(Key_metabolism_genes_cluster$percentage >= 33, "0.33", "0")))

#Remember the annoying thing that R sorts alphabetically? Let’s make sure we ahve the order we want.
#define order for the plot
Key_metabolism_genes_cluster$ClusterName2 <-  factor(Key_metabolism_genes_cluster$ClusterName, levels = rev(Cluster_order))
Key_metabolism_genes_cluster$Gene2 <-  factor(Key_metabolism_genes_cluster$Gene, levels = Genes_metabolism_order)
Key_metabolism_genes_cluster$pathway_2b <-  factor(Key_metabolism_genes_cluster$pathway_2, levels = Metabolism_Pathway_order)

#In the example here, we use a gradual scale. If we would want to use our 4 categories we can use this code #scale_fill_manual(values= c("white", "blue", "blue", "dodgerblue")) and replacing the fill = percentage with fill = category.

#plot
p1_metabolism <- 
  ggplot(Key_metabolism_genes_cluster, aes(x=Gene2, y=(ClusterName2))) + 
  geom_tile(aes(fill = percentage)) +
  facet_wrap( ~ pathway_2b, nrow = 1, scales='free_x') +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme_bw() +
  #scale_fill_manual(values= c("white", "blue", "blue", "dodgerblue")) +
  labs(x="", y="", fill="Percentage") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(legend.position="left",
        axis.text.x=element_text(angle=90,vjust = 1, hjust=1, size=6),
        #axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        plot.margin=unit(c(0, 0, 0, 0), "mm"))

p1_metabolism

#If we plot with facets we can sometimes have the problem that different genes have different widths. We can correct this behaviour with ggplotGrob.
# convert ggplot object to grob object (used to rescale plot)
gp_metabolism <- ggplotGrob(p1_metabolism)

# optional: take a look at the grob object's layout
gtable::gtable_show_layout(gp_metabolism)

# get gtable columns corresponding to the facets (5 & 9, in this case)
facet.columns <- gp_metabolism$layout$l[grepl("panel", gp_metabolism$layout$name)]

# get the number of unique x-axis values per facet (1 & 3, in this case)
x.var <- sapply(ggplot_build(p1_metabolism)$layout$panel_scales_x,
                function(l) length(l$range$range))

# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
gp_metabolism$widths[facet.columns] <- gp_metabolism$widths[facet.columns] * x.var

# plot result
#print
pdf("2_Output/ANME_AMOR_Annotation_heatmap_AOM_Vulcano_1.pdf", paper="special", family="sans",width=15, height=7, useDingbats=FALSE)
grid::grid.draw(gp_metabolism)
dev.off() 

while (!is.null(dev.list()))  dev.off()
