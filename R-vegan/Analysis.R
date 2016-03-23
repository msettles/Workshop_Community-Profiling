## install and load libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("Biobase","RColorBrewer","vegan","WGCNA","vegetarian","gplots","plotrix"))
library(Biobase)
library(RColorBrewer)
library(vegan)
library(WGCNA)
library(vegetarian)
library(gplots)
library(plotrix)


####################################################################
## Modify for the abundance files that come out of dbcAmplicons
####################################################################
abundancefile <- "16sV1V3.abundance.txt"
# contains at least the SAMPLEID column from the dbcAmplicons samplesheet 
metadatatable <- "slashpilesamplesheet.txt"
output = "16sv1v3" # prefex for output
source_level <- "genus" # chosen in dbcAmplicons abundance
####################################################################

## hacked version of heatmap
source("heatmap.microbe.R")
heatcol <- rainbow(256,start=0.2,end=1)

## SET up data
## Read in the abudance table, remove unnecessary columns, etc.
abund <- read.table(abundancefile,sep="\t",header=T,as.is=T)
rownames(abund) <- paste(abund$Taxon_Name,abund$Level) 
da <- dim(abund)
da
taxon_name <- abund$Taxon_Name
level <- abund$Level
abund$Taxon_Name <- NULL
abund$Level <- NULL

## read in the metadata
meta <- read.table(metadatatable,sep="\t",as.is=T, header=T)
abund <- abund[,match(meta$SampleID,colnames(abund))] # resort column on meta

#### check rows
rem <- rowSums(abund) > 0
table(rem)
if (length(rem) > 0){
  abund <- abund[-rem,]
  level <- level[-rem]
  taxon_name <- taxon_name[-rem]
}

## SAVE ORINGINAL DATA
abund.orig <- abund
level.orig <- level
taxon_name.orig <- taxon_name

### filter Samples, dumb filter, minimum number of reads per sample
remSample <- colSums(abund) >= 1000 & !is.na(colSums(abund))
table(remSample)

abund <- abund[,remSample]
meta <- meta[remSample,]

## Histogram of Read Counts per Sample
pdf(paste0(output,".readCountsHist.pdf"))
hist(colSums(abund),breaks=250,xlab="Number of Reads",main="Histgram of Sample Read Counts")
dev.off()

### compute frequencies (proportions)
freqAbund <- sweep(abund,2,colSums(abund),"/")

## filter out genus for table
### THE RULE: keep a taxon if it present (> 0) in more than 1/2 the samples or in 1 sample at a level of at leat 0.05%
#keepRows <- ((apply(freqAbund > 0.00,1,sum,na.rm=TRUE) >= ncol(freqAbund)/2) | (apply(freqAbund >= 0.05, 1, sum,na.rm=TRUE) > 0))
### THE RULE: keep a taxon if it present in more than 2 sample at a level of at least 0.01% or in 1 sample at a level of at leat 0.05%
keepRows <- ((apply(freqAbund >= 0.01,1,sum,na.rm=TRUE) > 2) | (apply(freqAbund >= 0.05, 1, sum,na.rm=TRUE) > 0))
keepRows["Bacteria domain"] <- FALSE
table(keepRows)


# Combine cut genera into 'other'
Other <- colSums(abund[!keepRows,])
abund <- abund[keepRows,]
level <- level[keepRows]
taxon_name <- taxon_name[keepRows]
abund <- rbind(abund,"Other Bacteria domain"=Other)
level <- c(level,"Bacteria")
taxon_name <- c(taxon_name,"root")

freqAbund <- sweep(abund,2,colSums(abund),"/")

#########################################################################################################################
#### WRITE OUT NEW DATA
### compute frequencies and write out
write.table(data.frame(Taxon_Name=taxon_name,Level=level,freqAbund),paste0(output,".reduced.proportions"),sep="\t",row.names=FALSE,col.names=TRUE,quote=F)
write.table(data.frame(Taxon_Name=taxon_name,Level=level,abund),paste0(output,".reduced.abundance"),sep="\t",row.names=FALSE,col.names=TRUE,quote=F)
#########################################################################################################################

######### Descriptive plots
# meta data has 3 metadata columns, Slash_pile_number, Dist_from_edge, Depth_cm all of which are factors
pile <- c("pink","red","grey")[as.numeric(as.factor(meta$Slash_pile_number))]
dist <- labels2colors(meta$Dist_from_edge)
depth <- c("salmon","orange")[as.numeric(as.factor(meta$Depth_cm))]
labelColors <- cbind(pile, dist, depth)

colnames(labelColors) <-   c("slashpile","dist","depth")

#### Standardize the abundance table using vegan, decostand, hellinger method, WHY? Because a biostatistician told me once it was the best
abund.standardized <- decostand(t(abund), method = "hellinger")
#### Compute the euclidean distances
abdist <- vegdist(t(abund))
#### Cluster using wards.D
hcl <- hclust(abdist,method="ward.D")

#cor(abdist, cophenetic(hcl))

plot(hcl)
rect.hclust(hcl,5)
grp <- cutree(hcl, 5)


####### Produce a dedrogram plot with annotation
pdf(file=paste0(output,".cluster-analysis-all.pdf"),width=11,height=5.5,pointsize=6)
  lColors <- labelColors
  subjectstring <- paste(meta$SampleID,meta$Slash_pile_number,meta$Dist_from_edge,meta$Depth_cm,sep="-")
  maintext <- "Cluster Analysis, all samples"
  plotDendroAndColors(hcl, 
                    colors = cbind(labels2colors(grp),lColors),groupLabels=c("group",colnames(labelColors)), rowText=cbind(grp,meta[,c("Slash_pile_number","Dist_from_edge","Depth_cm")]),
                    dendroLabels = subjectstring, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,cex.dendroLabels=0.8,
                    main = paste("veg distances and average linkage clustering","[",maintext,"]"))
  rect.hclust(hcl,5)
dev.off()


### ordered barplot of taxa representation
pdf(paste0(output,".taxaRepresentation.pdf"),width=7,height=10,pointsize=8)
  par(mar=c(3,15,2,1))
  maintext <- paste("Taxa represention across samples",sep="")
  boxplot(t(freqAbund[order(rowSums(freqAbund),decreasing=F),]),names=rownames(abund)[order(rowSums(freqAbund),decreasing=F)],horizontal=T,las=1,cex=0.3,cex.axis=0.75,
          main=maintext,cex.main=0.9)
dev.off()

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

### heatmap of the microbes
pdf(file=paste0(output,".microbeHeatMap.ordered.pdf"),bg="white",width=9,height=6.5,pointsize=10)
  par(oma=c(2,2,2,10))
  heatmap.microbe(as.matrix(freqAbund)[,hcl$order],
                col=heatcol,
                distfun=vegdist,
                labCol=subjectstring[hcl$order],
                ColSideColors=labels2colors(grp[hcl$order]),
                ColSideLabels="group",
                cexRow=0.7,
                cexCol=0.8,
                mar=c(7,5),
                dendrogram="none",
                Rowv=T,
                Colv=F,
                sepcolor="white",
                keysize=0.5,
                trace="none",
                key=T,
                density.info="none",family="sans")
dev.off()

## classical multidimentional scaling
cmdfit <- cmdscale(abdist,eig=TRUE, k=2) # k is the number of dim

# plot solution
pdf(file=paste0(output,"mds-AllData.pdf"),bg="white",width=10,height=7,pointsize=4)
  x <- cmdfit$points[,1]
  y <- cmdfit$points[,2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric MDS, Bray's Distance (All Data)", type="n")
  text(x, y, labels = subjectstring,col= labelColors[,"slashpile"])
dev.off()


### plot ordination, using meta MDS
ord <- metaMDS(t(abund))

pdf(file=paste0(output,"metricmdsLARGE-AllData.pdf"),bg="white",width=30,height=30,pointsize=8)
  plot(ord,type="n")
  text(ord,display="sites",cex=4,pch=21,col=labelColors[,"depth"],bg=labelColors[,"depth"])
  text(ord,display="spec",cex=2,col="black")
dev.off()

pdf(file=paste0(output,"metricmdsSmall-AllData.pdf"),bg="white",width=10,height=7,pointsize=12)
  plot(ord)
dev.off()

# fit environmental vairables onto the ordination
ord.fit <- envfit(ord ~ Slash_pile_number + Dist_from_edge + Depth_cm, data=meta, perm=999)
pdf(paste0(output,"envfit.MDSordination.pdf"),bg="white",width=10,height=10,pointsize=8)
  plot(ord, dis="site")
  plot(ord.fit)
  plot(ord, dis="spec")
  plot(ord.fit)
dev.off()

# constrained cononical ordination
ord.cca <- cca(t(abund) ~ Slash_pile_number + Dist_from_edge + Depth_cm, data=meta)
pdf(paste0(output,"envfit.CCAordination.pdf"),bg="white",width=10,height=10,pointsize=8)
  plot(ord.cca)
dev.off()

anova(ord.cca,by='term')
##################################################
############# DIVERSITY MEASURES IN VEGAN
# Shannon diversity index
shannon <- diversity(t(abund.orig))
shannon

summary(lm(shannon ~ as.factor(meta$Slash_pile_number) + as.factor(meta$Depth_cm) + as.factor(meta$Dist_from_edge)))

# Pielou’s evenness
pielou <- shannon/log(specnumber(t(abund.orig)))
pielou
summary(lm(pielou ~ as.factor(meta$Slash_pile_number) + as.factor(meta$Depth_cm) + as.factor(meta$Dist_from_edge)))

# R´enyi diversities
renyi <- renyi(t(abund.orig))
pdf(file=paste0(output,"renyi-diversities.pdf"),bg="white",width=7,height=7,pointsize=6)
  plot(renyi,cex=0.5)
dev.off()

# fisher's alpha
alpha <- fisher.alpha(t(abund.orig))
summary(lm(alpha ~ as.factor(meta$Slash_pile_number) + as.factor(meta$Depth_cm) + as.factor(meta$Dist_from_edge)))

# beta diversity
beta <- vegdist(t(abund.orig), binary=TRUE)

# species accumulation
sac <- specaccum(t(abund.orig))
pdf(file=paste0(output,"speciesAccumulation.pdf"),bg="white",width=7,height=7,pointsize=6)
  plot(sac, ci.type="polygon", ci.col="yellow")
dev.off()

