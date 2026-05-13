# Relative abundance stacked barplot

phy_rel_abund <- transform_sample_counts(phy, function(x) {x / sum(x)})

# Relative abundance Faavne

setwd("C:/Users/phr001/OneDrive - University of Bergen/PHD/Methods/Gallionella/Which_samples")
MAGs_abundance_wideformat <- read.table('Gallionella_Zeta_RelativeAbundance.csv', header=T, dec = '.', sep = ';')
MAGs_abundance_wideformat <- as.data.frame(MAGs_abundance_wideformat)

#Adjust the data from a wide format to a long format
#using tidyr: https://r-crash-course.github.io/14-tidyr/

library("tidyr")
library("dplyr")

str(MAGs_abundance_wideformat)

colnames(MAGs_abundance_wideformat)
rownames(MAGs_abundance_wideformat) <- c('Gallionella (f. Gallionellaceae)','Ghiorsea (c. Zetaproteobacteria)','Mariprofundus (c. Zetaproteobacteria)','other Gallionellaceae','other Zetaproteobacteria')
MAGs_abundance_wideformat_long <- data.frame(rows = row.names(MAGs_abundance_wideformat), stack(MAGs_abundance_wideformat))
View(MAGs_abundance_wideformat_long)

write.csv(MAGs_abundance_wideformat_long, file = "MAGs_abundance_wideformat_long.csv")

MAGs_abundance_wideformat_long2 <- read.table('MAGs_abundance_wideformat_long2.csv', header=T, dec = '.', sep = ';')
MAGs_abundance_wideformat_long2 <- as.data.frame(MAGs_abundance_wideformat_long2)
MAGs_abundance_wideformat_long2
View(MAGs_abundance_wideformat_long2)

# Plotting on the Domain level

library(ggplot2)

MAGs <- ggplot(data=MAGs_abundance_wideformat_long2, aes(x=ind, y=values, fill=rows)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand=c(0, 0)) +
  theme_classic() +
  labs (x=NULL, y="Mean Relative Abundance (MAGs)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Mean Relative Abundance (MAGs)") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  guides(fill = guide_legend(title = "Taxa")) + 
  scale_fill_manual(values=c("#d95f0e","#fec44f","#a6bddb","#ece7f2","#2b8cbe"),
                    limits = c('Gallionella (f. Gallionellaceae)','other Gallionellaceae','Mariprofundus (c. Zetaproteobacteria)','Ghiorsea (c. Zetaproteobacteria)','other Zetaproteobacteria'))

MAGs

#Plot abundances (frames) for each separate MAG as nested category (fill=rows, color/alpha=MAGs)
MAGs_abundance_wideformat_long2_nested <- read.table('MAGs_abundance_wideformat_long2_nested.csv', header=T, dec = '.', sep = ';')
MAGs_abundance_wideformat_long2_nested <- as.data.frame(MAGs_abundance_wideformat_long2_nested)
MAGs_abundance_wideformat_long2_nested
View(MAGs_abundance_wideformat_long2_nested)

#to Increase number of axis ticks
#following this answer https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
dat <- data.frame(x = rnorm(100), y = rnorm(100))

##USE THIS ONE (COLOR) - frame around the nested category
MAGs2 <- ggplot(data=MAGs_abundance_wideformat_long2_nested, aes(x=ind, y=values, fill=rows, color=MAGs)) +
  geom_bar(stat = "identity", color="black", width = 0.9,size=0.1) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_classic() +
  labs (x=NULL, y="Relative Abundance (%) - MAGs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Relative Abundance (%) - MAGs") +
  scale_y_continuous(labels=function(x)x*100, expand=c(0,0), limits = c(0,0.61), breaks = round(seq(min(dat$y), max(dat$y), by = 0.1),1)) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  labs(tag = "A", size=20,face="bold") +
  theme(plot.tag = element_text(size=20,face="bold")) +
  #scale_x_discrete(name ="Faavne         Jan Mayen") +
  #theme(axis.title.x = element_text(size=15,face="italic")) +
  guides(fill = guide_legend(title = "Taxa")) + 
  scale_fill_manual(values=c("#d95f0e","#fec44f","#a6bddb","#ece7f2","#2b8cbe"),
                    limits = c('Gallionella (f. Gallionellaceae)','other Gallionellaceae','Mariprofundus (c. Zetaproteobacteria)','Ghiorsea (c. Zetaproteobacteria)','other Zetaproteobacteria'))
#guides(fill="none")


MAGs2

