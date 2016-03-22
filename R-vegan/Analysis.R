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

abundancefile <- "16sV1V3.abundance.txt"
metadatatable <- "slashpilesamplesheet.txt"
output = "16sv1v3"

source("heatmap.microbe.R")
heatcol <- rainbow(256,start=0.2,end=1)

## SET up data
## Read in, remove unnecessary columns, etc.

source_level <- "genus"
abund <- read.table(abundancefile,sep="\t",header=T,as.is=T)

rownames(abund) <- paste(abund$Taxon_Name,abund$Level) 
da <- dim(abund)
taxon_name <- abund$Taxon_Name
level <- abund$Level
abund$Taxon_Name <- NULL
abund$Level <- NULL

meta <- read.table(metadatatable,sep="\t",as.is=T, header=T)

abund <- abund[,match(meta$SampleID,colnames(abund))] # resort column on meta

#### check rows
rem <- which(rowSums(abund)==0)
if (length(rem) > 0){
  abund <- abund[-rem,]
  level <- level[-rem]
  taxon_name <- taxon_name[-rem]
}

## SAVE ORINGINAL DATA
abund.orig <- abund
level.orig <- level
taxon_name.orig <- taxon_name

### filter Samples (NOT YET)
remSample <- colSums(abund) >= 1000 & !is.na(colSums(abund))
table(remSample)

abund <- abund[,remSample]
meta <- meta[remSample,]

## Histogram of Read Counts per Sample
pdf(paste0(output,".readCountsHist.pdf"))
hist(colSums(abund),breaks=50,xlab="Number of Reads",main="Histgram of Sample Read Counts")
dev.off()

### comput frequencies
freqAbund <- sweep(abund,2,colSums(abund),"/")

## filter out genus for table
### THE RULE: keep a taxon is it present in more than 1 sample at a level of at least 0.01% or in 1 sample at a level of at leat 0.05%
keepRows <- ((apply(freqAbund >= 0.01,1,sum,na.rm=TRUE) > 1) | (apply(freqAbund >= 0.05, 1, sum,na.rm=TRUE) > 0)) ## remove taxa with < 0.05 and only found in 1
keepRows["Bacteria domain"] <- FALSE
table(keepRows)


# Combine non species
Other <- colSums(abund[!keepRows,])
abund <- abund[keepRows,]
level <- level[keepRows]
taxon_name <- taxon_name[keepRows]

abund <- rbind(abund,"Other Bacteria"=Other)
level <- c(level,"Bacteria")
taxon_name <- c(taxon_name,"root")
freqAbund <- sweep(abund,2,colSums(abund),"/")

#########################################################################################################################
#### WRITE OUT NEW DATA
### compute frequencies and write out
write.table(data.frame(Names=rownames(freqAbund),Level=level,freqAbund),paste0(output,".reduced.proportions"),sep="\t",row.names=FALSE,col.names=TRUE,quote=F)
write.table(data.frame(Names=rownames(abund),Level=level,abund),paste0(output,".reduced.abundance"),sep="\t",row.names=FALSE,col.names=TRUE,quote=F)
#########################################################################################################################



##################################################################################################
### Paper 1 based Figures visit 1 to 3 only
##################################################################################################

Organic <- c("pink","red","grey")[as.numeric(as.factor(meta$Organic))]
SolidLiquid <- c("salmon","orange")[as.numeric(as.factor(meta$SolidLiquid))]
Type <- labels2colors(meta$TYPE)
labelColors <- cbind(Type, SolidLiquid, Organic)

colnames(labelColors) <-   c("Type","Solid/Liquid","Organic")

abund.standardized <- decostand(t(splitAbund), method = "hellinger")
abdist <- dist(abund.standardized, method="euc")

maintext <- "Cluster Analysis, all samples"

hcl <- hclust(abdist,method="ward.D")

lColors <- labelColors
subjectstring <- paste(meta$SAMPLE_ID,meta$TYPE,meta$SolidLiquid,meta$Organic,sep="-")


pdf(file="cluster-analysis-all.pdf",width=11,height=5.5,pointsize=6)
plotDendroAndColors(hcl, 
                    colors = lColors,groupLabels=names(labelColors), rowText=meta[,c("TYPE","SolidLiquid","Organic")],
                    dendroLabels = subjectstring, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,cex.dendroLabels=0.8,
                    main = paste("Hellinger standardization with euclidean distances and average linkage clustering","[",maintext,"]"))
dev.off()



renamed <- rownames(abund)

freqa <- sweep(abund,2,colSums(abund),"/")
pdf("TaxaRepresentation.pdf",width=7,height=10,pointsize=8)
par(mar=c(3,15,2,1))
  maintext <- paste("Taxa represention across samples",sep="")
  boxplot(t(freqa[order(rowSums(freqa),decreasing=F),]),names=renamed[order(rowSums(freqa),decreasing=F)],horizontal=T,las=1,cex=0.3,cex.axis=0.75,
          main=maintext,cex.main=0.9)

dev.off()

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

rsum <- rowSums(freqa)

pdf(file="MicrobeHeatMap.ordered.pdf",bg="white",width=9,height=6.5,pointsize=10)
par(oma=c(2,2,2,10))
heatmap.microbe(as.matrix(freqa)[,hcl$order],
                col=heatcol,
                distfun=vegdist,
                labCol=subjectstring[hcl$order],
                ColSideColors=labelColors[,"Type"],
                ColSideLabels="Type",
                cexRow=0.7,
                cexCol=0.8,
                mar=c(7,5),
                dendrogram="none",
                Rowv=F,
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
pdf(file="mds-AllData.pdf",bg="white",width=10,height=7,pointsize=4)
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS, Bray's Distance (All Data)", type="n")
text(x, y, labels = subjectstring,col= labelColors[,"Type"])
dev.off()


### plot ordination, using meta MDS
ord <- metaMDS(t(abund))

pdf(file="metricmdsLARGE-AllData.pdf",bg="white",width=30,height=30,pointsize=8)
plot(ord,type="n")
points(ord,display="sites",cex=2,pch=21,col=labelColors[,"Type"],bg=labelColors[,"Type"])
text(ord,display="spec",cex=1,col="black")
dev.off()

pdf(file="metricmdsSmall-AllData.pdf",bg="white",width=10,height=7,pointsize=12)
plot(ord)
dev.off()


############# DIVERSITY MEASURES IN VEGAN

# Shannon diversity index
shannon <- diversity(t(abund))

# Pielou’s evenness
pielou <- shannon/log(specnumber(t(abund)))

# R´enyi diversities
renyi <- renyi(t(abund))
pdf(file="renyi-diversities.pdf",bg="white",width=7,height=7,pointsize=6)
plot(renyi,cex=0.5)
dev.off()

# fisher's alpha
alpha <- fisher.alpha(t(abund))

quantile(rowSums((t(abund))))

# species accumulation
sac <- specaccum(t(abund))
pdf(file="speciesAccumulation.pdf",bg="white",width=7,height=7,pointsize=6)
plot(sac, ci.type="polygon", ci.col="yellow")
dev.off()

# beta diversity
beta <- vegdist(t(abund), binary=TRUE)
