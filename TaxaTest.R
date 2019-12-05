#SETUP-----------
#install.packages("lookupTable")
library(lookupTable)
#install.packages("plyr")
library(plyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("magrittr")
library(magrittr)

data=read.csv("feature-table-6.csv")
#Remove controls as log(0)=Inf
datasub<-subset(data, auc_continuous!=0)
taxakey=read.csv("Level6TaxaMetadata1.csv")

#Plotting AUC vs Observed OTUs-----------

#Scatter plot of relative abundance of Taxa952 in relation to the auc/infection severity
#install.packages("tidyverse")
library(tidyverse)

pdf('Taxa.pdf')
ggplotRegression <- function(dat, xvar, yvar) {
  fml <- paste(yvar, "~", xvar)
  fit<-lm(fml,dat)
  ggplot(data=datasub, aes(x=observed_otus,y=auc)) + geom_point() +
    labs(x = "Observed OTUs") +
    labs(y = "log AUC") +
    labs(title = "Analyzing the relationship between parasite infection \n severity and skin bacterial microbiome diversity \n in P. reticulata", face = "bold") +
    scale_y_log10(breaks = trans_breaks("log10",function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, formula = y ~ x) +
    geom_text(aes(120, 3, label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                                           " P =",signif(summary(fit)$coef[2,4], 5)), size = 6)) +
    geom_text(aes(120, 2, label = paste ("Intercept =",signif(fit$coef[[1]],5 ),
                                           " Slope =",signif(fit$coef[[2]], 5)), size = 6))+
    theme(axis.text = element_text(size = 12))+
    theme(axis.title = element_text(size = 15),
          legend.position='none')+
    theme(plot.title = element_text(hjust = 0.5, size = 15))
}
ggplotRegression(datasub, "auc", "observed_otus")
dev.off()

help(ggtitle)

#FULL-DATASET---------------
#Store the column names of all the taxa into a variable (-1-26) excluding metadata columns
Taxa<-colnames(datasub[-c(1:26)])
#Have to organize the row names by the order the ANOVA outputs, otherwise values will not match their name
x<-data.frame(row.names= c("sex","line","batch","behav","TaxaN","seqbatch","sex:line"))

#For every bacterial taxa, conduct a linearized model to determine whether the relative abundance of that specific taxa has
#a significant affect on the log of the area under the curve (auc, infection severity), while taking into account
#factors such as sex, line, batch, etc.
for(i in Taxa){
  tempcolindex<-as.numeric(substr(i,5,nchar(i)))
  tempcolindex<- tempcolindex+26
  tempdatasub<-subset(datasub, select = c(1:26,tempcolindex))
  names(tempdatasub)[27]<-"TaxaN"
  templm<-lm(log(auc_continuous)~sex+line+sex:line+batch+behav+TaxaN+seqbatch, data=tempdatasub)
  tempanova<-anova(templm)
  x<-cbind(x,tempanova$`Pr(>F)`[1:7])
}

colnames(x)<-Taxa

#Number of taxa that should be removed in subsequent removedx df
rowSums(is.na(x))

#Filtering out the NA taxa should still result in removing the singularity-taxa.. since the NA result due to that
#The NAs are in the sex:line row because the ANOVA moves the value for that up where the TaxaN was
removedx<-na.omit(t(x))
removedx<-t(removedx)

#Num. sig taxa before correction
uncorrected<-removedx["TaxaN",]<=0.05
table(uncorrected)["TRUE"]
#In this case, 102 significant prior to correction

#Subset and extract only row 5 (anova p values for the Taxa)
removedxtaxa<-removedx[c(5), , drop=FALSE]
taxacol<-colnames(removedxtaxa)
#conduct FDR correction for multiple comparisons
#change ncol based on removedxtaxa column num
fdrcorrection<-(matrix(p.adjust(as.vector(as.matrix(removedxtaxa)), method='fdr'),ncol=833))

dfr<-as.data.frame(t(fdrcorrection))
dtaxa<-as.list(t(taxacol))

dfr$taxa<-dtaxa

#Add Taxa by taxonomic name
dfr$name<-taxakey$Name[match(dfr$taxa,taxakey$ID)]

#How many taxa have a p value less than or equal to 0.05 after FDR correction?
tm<-dfr[,"V1"]<=0.05
table(tm)["TRUE"]
#print(tm)
#Which taxa are those?
tempsigdata<-which(dfr[,"V1"] <=0.05, arr.ind=TRUE)
#returns row numbers as well as Taxa numbers

#Subset dfr by including only rows in tempsigdata, export
dfrsub<-dfr[tempsigdata, , drop=FALSE]
dfrsub <- apply(dfrsub,2,as.character)
write.csv(dfrsub, "full-analysis-fdr-corrected.csv")
