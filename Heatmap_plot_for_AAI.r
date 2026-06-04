#Heatmap plotting for average nucleotide identity (ANI), average aminoacid identity (AAI), and alignment fraction (AF)

#Plotting AAI values from EzAAI tool output

setwd("C:/Users/phr001/OneDrive - University of Bergen/PHD/Supervision/2025_Emilie/LC_904_LAC/")

##Making a matrix from EzAAI results
library(tidyr)
aai_tomatrix <- read.delim("C:/Users/phr001/OneDrive - University of Bergen/PHD/Supervision/2025_Emilie/LC_904_LAC/aai_2cut.txt", header=TRUE)

##change header of column 2
colnames(aai_tomatrix)[2] <- "MAGs"
pivot_wider(aai_tomatrix, names_from = Label.1, values_from = AAI)
write.csv(pivot_wider(aai_tomatrix, names_from = Label.1, values_from = AAI), file = "AAI_EzAAI_matrix.csv", row.names = FALSE)

#Plotting AAI in a heatmap
#Remove V2 in first row and column
aai_matrixfile_EzAAI <- read.delim("C:/Users/phr001/OneDrive - University of Bergen/PHD/Rwork/AAI_EzAAI_matrix.csv", header=TRUE, sep = ",", row.names=1)

View(aai_matrixfile_EzAAI)

sapply(aai_matrixfile_EzAAI, class)

# https://stackoverflow.com/questions/37707060/converting-data-frame-column-from-character-to-numeric/37707117
AAI_EzAAI <- as.data.frame(lapply(aai_matrixfile_EzAAI, function(x) as.numeric(as.character(x))))
rownames(AAI_EzAAI) <- rownames(aai_matrixfile_EzAAI)

# And if you check again this time you get numericals
sapply(AAI_EzAAI, class)
# Now make you matrix
y <- data.matrix(AAI_EzAAI)
View(y)

# And your plot!

#set colour key cut-off: yellow for under 65% AAI, orange for over 65 AAI (genus threshold), red for over 95 AAI (species)
my_palette <- colorRampPalette(c("yellow", "orange", "#E4001A"))(n = 6)
col_breaks=c(seq(0,64.999,length=1)
              ,seq(65,94.99,length=3)
              ,seq(95,100,length=3))

# Heatmap with a small colour key without count, larger font size (cexRow, cexCol)
library(gplots)
library(viridis)
heatmap.2(y
          ,dendrogram="none"
          ,Rowv=FALSE
          ,Colv=FALSE
          ,trace='none'
          ,margins = c(15,15)
          ,cexRow = 1
          ,cexCol = 1
          ,symm = FALSE
          ,key = TRUE
          ,key.title = "AAI percentage identity"
          ,main = "AAI percentage identity"
          ,density.info = "none"
          ,symkey = FALSE
          ,densadj = 0
          ,tracecol = NA
          #,col=heat.colors(6,rev = TRUE)
          ,keysize = 0.5
          ,lhei = c(1,20)
          ,lwid = c(1,10)
          ,key.ylab = "Genome pairs"
          ,key.xtickfun = NULL
          ,key.ytickfun = NULL
          ,key.xlab = "AAI percentage identity"
          ,key.par=list(mar=c(1,1,1,1), cex=1)
          #,lwid = c(1,15)
          #,lhei = c(1,15)
          ,col=my_palette
          ,breaks=col_breaks,
          na.color = "Green")

dev.off()

#clustering
#dendrogram = c("both","row","column","none")
heatmap.2(y
          ,dendrogram="both"
          ,trace='none'
          ,margins = c(14.5,14.5)
          ,cexRow = 1
          ,cexCol = 1
          ,symm = FALSE
          #,key = TRUE
          #,key.title = "AAI percentage identity"
          #,main = "AAI percentage identity"
          ,density.info = "none"
          #,symkey = FALSE
          #,densadj = 0
          #,tracecol = NA
          #,col=heat.colors(6
          #  ,rev = TRUE)
          ,keysize = 0.5
          ,lhei = c(1,20)
          ,lwid = c(1,10)
          #,key.ylab = "Genome pairs"
          ,key.xtickfun = NULL
          ,key.ytickfun = NULL
          ,key.xlab = "AAI percentage identity"
          ,key.par=list(mar=c(1,1,1,1), cex=1)
          #,lwid = c(1,15)
          #,lhei = c(1,15)
          ,col=my_palette
          ,breaks=col_breaks)

#plot values on the heatmap
#Change Number of Digits in Global R Options
# Modify global options
options(digits = 2)     

dev.off()
heatmap.2(y
          ,dendrogram="both"
          ,trace='none'
          ,margins = c(15,15)
          ,cexRow = 1
          ,cexCol = 1
          ,symm = FALSE
          #,key = TRUE
          #,key.title = "AAI percentage identity"
          #,main = "AAI percentage identity"
          ,density.info = "none"
          #,symkey = FALSE
          #,densadj = 0
          #,tracecol = NA
          #,col=heat.colors(6
          #  ,rev = TRUE)
          ,keysize = 0.5
          ,lhei = c(1,20)
          ,lwid = c(1,10)
          #,key.ylab = "Genome pairs"
          ,key.xtickfun = NULL
          ,key.ytickfun = NULL
          ,key.xlab = "AAI percentage identity"
          ,key.par=list(mar=c(1,1,1,1), cex=1)
          #,lwid = c(1,15)
          #,lhei = c(1,15)
          ,col=my_palette
          ,breaks=col_breaks
          ,cellnote = y
          ,notecex=0.8
          ,notecol = "black")

