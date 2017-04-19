#All the libraries you will need
library(dplyr)
library(reshape)
library(ggplot2)
library(latticeExtra)
library(pheatmap)
library(phyloseq)
library(ggpmisc)
library(Hmisc)
library(vegan)
library(Kendall)
library(corrplot)

##Heatmap of hmms using pheatmap, similar code used for Fig.4 and Fig.5
#Fig4
counts.bin <- read.table("counts_plot.bins.norm.txt", row.names="bin", check.names=FALSE, header=TRUE)
matrix.bin <- data.matrix(counts.bin, rownames.force=TRUE)
#have NA, do log10(n+1)
log10.matrix <-log10(matrix.bin+1)
plot.bin <- pheatmap(matrixlog.bin, cluster_cols= FALSE, cluster_rows = FALSE,color=colorRampPalette(c("white", "lightblue","darkblue"))(50), border_color="black",fontsize_row=8, fontsize_col=8, cellwidth=20, cellheight=8, fontsize=8, annotation_legend=NA)

#Fig5
counts.bin <-read.table("hmm_heatmap_paper.v6.txt", row.names="function", check.names=FALSE, header=TRUE)
matrix.bin <- data.matrix(counts.bin, rownames.force=TRUE)
plot.bin <- pheatmap(matrix.bin, cluster_cols= FALSE, cluster_rows = FALSE,color=colorRampPalette(c("honeydew2","green4"))(50), border_color="black",fontsize_row=16, fontsize_col=16, cellwidth=20, cellheight=20, fontsize=16, annotation_legend=NA)

#Subsetting different colors to make 0/NA one color and use gradient for values above 0/NA
#Set breaks, set length to number of gradients you want
matrix.min=min(sqrt.matrix[sqrt.matrix >0])
matrix.max=max(sqrt.matrix[sqrt.matrix >0])
pairs.breaks <-c(seq(0, length=1), seq(matrix.min, matrix.max, length.out=26))
#assign colors, set gradient to one more than length of breaks (else it will repeat colors)
col1 = rep("black")
col2 = colorRampPalette(brewer.pal(9,"Reds"))(27)
colors <-c(col1,col2)

##line plot Depth distributions of key hmms, Fig.6
#Put all the hmms counts into a table (bin,pathway,gene,count1,count2,count3,count4,count5)
hmm <-read.table("hmm_function_combo.txt", header=TRUE, check.names=FALSE)
hmmmelt <-melt(hmm, id.vars= c("bin", "path", "gene"))
hmmmelt$variable <- as.numeric(hmm, melt$variable)
#manually assign colors for Rubisco and N
hmmcolors <-c("green4", "red", "yellow2", "purple4", "lawngreen", "orange", "magenta", "cyan", "coral","wheat3", "chocolate4", "plum4","goldenrod4","mistyrose1", "sienna3", "turquoise4", "mediumspringgreen", "orangered3", "tomato4", "darkblue", "darkmagenta", "black", "grey40", "yellowgreen", "hotpink","darkolivegreen","honeydew3")
hmmplot <- ggplot(hmmmelt, aes(x=variable, y=value, color=bin, shape=gene)) + geom_line() + geom_point(size=2) + coord_flip() + scale_x_reverse() + facet_wrap(~path)
hmmplot + scale_color_manual(values=hmmcolors) + theme_bw() + scale_shape_manual(values=c(16,17,19,20,16)) +scale_y_log10()

##line plot depth distribution of 16S sequences, Fig. 7AB
#Use the mothur classification file.Format with depth rows (with 1 or 4 appended and bins as cols)
#Add gray boxes in illustrator because ain't nobody gots time for ggplot aestetics
#calculate mean and SD from melted
phyla <-read.table("LamMg.2.mat.v2.tax.phyla.txt", header=TRUE, check.names=FALSE)
meltphyla <-melt(phyla, id.vars=c("depth","time"))
dplyphyla <- ddply(meltphyla, c("depth","time","variable"), summarise, mean= mean(value), sd= sd(value))
#scale color for 0 0.1%, may need to alter
phylumcolors <-c("lightpink2", "brown", "yellow", "goldenrod2", "goldenrod4", "purple4", "lawngreen","green4", "darkviolet","slategray1","magenta","darkolivegreen","red","blue4","orange","bisque3", "cyan", "deeppink4", "black","grey30","deepskyblue1")
#plot, need to add breaks to yaxis to match figure in paper
pp <- ggplot(dplyphyla, aes(x=depth, y=mean, color=variable)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) + geom_line() + geom_point()+scale_color_manual(values=phylumcolors)
pp + facet_wrap(~time) + theme_bw() + theme(legend.position="none", strip.text=element_text(size=8), axis.text.x=element_text()) + coord_flip() + xlim(limits=c(6.5,0))

##NMDS plot done in phyloseq for figure 7C
#read in a biom file created in mothur (can do this with qiime too)
#load biome table
biome <- import_biom("LamMg.2.mat.0.03.subsample.0.03.biom", parseFunction=parse_taxonomy_default)
#counts to percent abundance
Lmg.abund.5mm <- transform_sample_counts(biome, function(x) x / sum(x))
#only those OTU who's percent abundance is greater than 0.01% across all samples, remove those from dataset(TRUE)
Lmg.hi.abund.5mm <- filter_taxa(Lmg.abund.5mm, function(x) mean (x) > 1e-4, TRUE)
#omit specific samples, in this case unmatched day and night
above3.75 <- subset_samples(Lmg.hi.abund.01, !depth%in%c("4","4.25","4.5","4.75","5","5.25","5.5","5.75","6","6.25","6.5"))
#Ordination
nmds <- ordinate(above3.75, method='NMDS', distance="bray", trymax=100)
print(nmds) #give stress and how good
stressplot(nmds) #shepardsplot
#NMDS plot
getPalette = colorRampPalette(c("red","yellow", "green4","blue"))
plot_ordination(above3.75, nmds, type="samples", color="depth", shape="time") + geom_point(size=3) + scale_color_manual(values=getPalette(24))


##Fig S1 was generated by Hans and Jim

##Plotting KEGG pathways, Fig S2
KO <-read.table("CA.KO.path.abund.txt", header=TRUE, check.names=FALSE)
KOmelt<-melt(KO, id.vars=c("path3"))
KOplot <-ggplot(KOmelt,aes(y=value, x=path3, fill=variable)) + geom_bar(stat="identity", position=position_dodge(-.9)) +theme_bw()
KOplot + scale_fill_manual(values=KOcolor, breaks=c("5","4","3","2","1")) + theme(axis.text.y= element_text(size="8") , panel.grid.major.x=element_blank(),legend.position="none") + scale_y_continuous(expand=c(0,0), limits=c(0,0.2), breaks=c(.05,.1,.15,.2)) + coord_flip()


#Regression lines plots that are dubious, Figure S3
hmm <-read.table("hmm_matches_TPM.unique.2.txt", header=TRUE, check.names=FALSE)
hmmmelt <-melt(hmm, id.vars=c("bin", "taxonomy", "complete","model","gene", "path"))
#replace 0 values with nothing so they won't plot but keep in dataset so all depths seen in plots
hmmmelt[,8][hmmmelt[,8] == 0] <- NA
#Sum the counts per gene model for cyano vs. non-cyano
#Add a column to denote dominant cyanos from others
bin.values <-c("C","C","C")
bin.index <-c("2","13","14")
hmmmelt.cyano <-hmmmelt
hmmmelt.cyano$type <- bin.values[match(hmmmelt.cyano$bin, bin.index)]
hmmmelt.cyano$type[is.na(hmmmelt.cyano$type)] <- "N"
hmmmelt.cyano$type <-as.factor(hmmmelt.cyano$type)
#remove NA since sum can't use those
hmmmelt.cyano <-na.omit(hmmmelt.cyano)
#Sum
melt.cyano.sum<-ddply(hmmmelt.cyano,.(gene,path,variable,type), summarize, sum=sum(value))
#plot this (have SE, next plot without )
hmmplot <- ggplot(melt.cyano.sum, aes(x=variable, y=sum, color=type)) + geom_point(size=1) + geom_smooth(aes(group=type, fill=type), method="lm", size=0.5) + facet_wrap(path~gene, scale="free") + theme_bw()
hmmplot + scale_color_manual(values=c("green4","red")) + scale_fill_manual(name="type", values=c("lightgreen","lightpink"))
#wo SE since sum
hmmplot <- ggplot(melt.cyano.sum, aes(x=variable, y=sum, color=type)) + geom_point(size=1) + geom_smooth(aes(group=type), method="lm", size=0.5, se = FALSE) + facet_wrap(path~gene, scale="free") + theme_bw()
hmmplot + scale_color_manual(values=c("green4","red"))
#Pearson correlations on the same data
#plot correlations for all the plots displayed using ddply, omit all NA
hmm.omitna <-na.omit(hmmmelt)
#Setup cor.test as a function,use pearson
pcorfun = function(x, y) {
  corr=(cor.test(x, y,method="pearson"))
}
#translate p-values to symbol notation
formatPvalues =function(pvalue) {
  ra<-""
  if(pvalue <= 0.1) ra<-"."
  if(pvalue <= 0.05) ra<-"*"
  if(pvalue <= 0.01) ra<-"**"
  if(pvalue <= 0.001) ra<-"***"
  return(ra)
}
#ddply for the magics
p.cors <-ddply(hmm.omitna,.(path, gene), summarize, r=round(pcorfun(depth, value)$estimate,2), pval=formatPvalues(pcorfun(depth, value)$p.value))

#Add pearson correlations (using code above) since we did a linear regression
melt.cyano.sum$variable = as.numeric(melt.cyano.sum$variable)
p.cors <-ddply(melt.cyano.sum,.(path, gene,type), summarize, r=round(pcorfun(variable, sum)$estimate,2), pval=formatPvalues(pcorfun(variable, sum)$p.value))
hmmplot + geom_text(data=p.cors, aes(x=-Inf,y=Inf, label=paste("r=", r, pval, sep=""), hjust=-0.2, vjust=1.2)) +scale_color_manual(values=c("green4","red")) + scale_fill_manual(name="type", values=c("lightgreen","lightpink"))

##Coabundance of genomes and biogeochem measurements, Fig. S4
#Create a file just containing the genome counts (rounded to nearest whole number) for the 5 different sections (omit the names of these)
bin <-read.table("complete.genome.coverage.txt", header=TRUE, check.names=FALSE)
#Spearman correlation for plotting
bin.cor <-cor(bin, method="spearman")
#Spearman correlation to get rho and significance
bin.rcoor = rcorr(as.matrix(bin), type="spearman")
#function to flatten coor matrix into columns
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
flat.rcoor = flattenCorrMatrix(bin.rcoor$r, bin.rcoor$P)
write.csv(flat.rcoor, "spearman.complete.genomes.rcoor.txt")
#Let's play with plotting
rcoor.P = bin.rcoor$P
col <-colorRampPalette(c("blue","white","red"))(20)
corrplot(bin.cor, method="circle", type="upper", order="hclust", col=col, bg="white")
#Sorted order
corrplot(bin.cor, method="square", type="upper", l=col, bg="white", p.mat=rcoor.P, sig.level=0.05, insig="blank", tl.col="black", tl.cex=0.5)


##Bar graphs of depth distribution hmms, Fig.S5
suppmelt <-melt(supp, id.vars= c("bin", "path", "gene"))
suppmelt$variable <- as.numeric(suppmelt$variable)
supplot <- ggplot(suppmelt, aes(x=variable, y=value, fill=bin)) + geom_bar(stat="identity")  + facet_wrap(~path+gene, scale="free")
supplot + theme_bw() + scale_fill_manual("Legend", values=bins.col2)

#Want to subset out facets but keep the bins the same color
#Grab unique number of bins from dataframe
unq.bins <- unique(suppmelt$bin)
#Select your palette
getPalette=colorRampPalette(brewer.pal(8,"Dark2"))
#Make a dataframe connecting colors to the number of bins and bin names
bins.col2 <-getPalette(length(unq.bins))
names(bins.col2) <- unq.bins
#to assign specific colors to specific bins, need to do for each then link to names
bins.col2 <-str_replace(bins.col2,"#666666","black")
bins.col2<-str_replace(bins.col2,"#1B9E77","green1")
bins.col2<-str_replace(bins.col2,"#3C9262","darkblue")
bins.col2<-str_replace(bins.col2,"#5D874E","deepskyblue")
bins.col2<-str_replace(bins.col2,"#4C8D58","darkred")
bins.col2<-str_replace(bins.col2,"#7E7C39","chartreuse3")
bins.col2<-str_replace(bins.col2,"#A3751F","cyan4")
bins.col2<-str_replace(bins.col2,"#966A77","red")
bins.col2<-str_replace(bins.col2,"#8F772F","gray87")

##Plotting the GH percent of proteins in bins, Figure S6
GH <- read.table("CA.GH.polysac.abund.txt", header=TRUE, check.names=FALSE)
GHmelt <-melt(GH, id.vars=c("bins"))
GHplot <-ggplot(GHmelt,aes(x=bins, y=value, fill=variable)) + geom_bar(stat="identity") +theme_bw()
GHcolors <-c("forestgreen", "black","yellow","cyan","brown","orange","purple","darkblue","red","chartreuse1")
GHplot + scale_fill_manual(values=GHcolors) + theme(axis.text.x= element_text(angle=90, size="7"), panel.grid.major.x=element_blank(), legend.position="none") + scale_y_continuous(expand=c(0,0), limits=c(0,9), breaks=c(0,1,2,3,4,5,6,7,8,9))

#Plotting distribution of 16S mapped to CA, Figure S7A and S7B
#Be sure to treat 0 and NA correctly since summarizing
data <-read.table("CA.v2.binsmatch.txt", header=TRUE, check.names=FALSE)
data[data ==0] <- NA
meltdata <- melt(data, id.vars=c("depth", "time"))
meltdata <- na.omit(meltdata)
dplydata <- ddply(meltdata, c("depth","time","variable"), summarise, mean= mean(value), sd= sd(value))
psubbin <- ggplot(dplydata, aes(x=depth, y=mean, color=time)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) + geom_line() + geom_point() + scale_color_manual(values=c("red","blue"))
psubbin + facet_wrap(~variable, scale="free", ncol=6) + theme_bw() + theme(legend.position="none", strip.text=element_text(size=8), axis.text.x=element_text()) + scale_x_continuous(limits=c(0,4))


