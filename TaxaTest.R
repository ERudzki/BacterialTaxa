#PLOTS-----------
#Plotting AUC vs Observed OTUs
#install.packages("ggplot2")
library(ggplot2)

#Scatter plot of relative abundance of Taxa952 in relation to the auc/infection severity
data=read.csv("feature-table-6.csv")
datasub<-subset(data, auc_continuous!=0)

#install.packages("tidyverse")
library(tidyverse)

ggplotRegression <- function(dat, xvar, yvar) {
  fml <- paste(yvar, "~", xvar)
  fit<-lm(fml,dat)
  ggplot(data=readssub, aes(x=observed_otus,y=auc)) + geom_point() +
    scale_y_log10(breaks = trans_breaks("log10",function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, formula = y ~ x) +
    geom_text(aes(120, 1000, label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                                           "Intercept =",signif(fit$coef[[1]],5 ),
                                           " Slope =",signif(fit$coef[[2]], 5),
                                           " P =",signif(summary(fit)$coef[2,4], 5))))
}
ggplotRegression(datasub, "auc", "observed_otus")

#SETUP-----------
#install.packages("lookupTable")
library(lookupTable)
#install.packages("plyr")
library(plyr)

data=read.csv("feature-table-6.csv")
#Remove controls as log(0)=Inf
datasub<-subset(data, auc_continuous!=0)
taxakey=read.csv("Level6TaxaMetadata1.csv")

#FULL-DATASET---------------
#Store the column names of all the taxa into a variable (-1-26) excluding metadata columns
Taxa<-colnames(datasub[-c(1:26)])
#Things tested as row names
#Have to organize the row names by the order the ANOVA outputs, otherwise values will not match their name
x<-data.frame(row.names= c("sex","line","batch","behav","TaxaN","seqbatch","sex:line"))
for(i in Taxa){
  tempcolindex<-as.numeric(substr(i,5,nchar(i)))
  tempcolindex<- tempcolindex+26
  tempdatasub<-subset(datasub, select = c(1:26,tempcolindex))
  names(tempdatasub)[27]<-"TaxaN"
  templm<-lm(log(auc_continuous)~sex+line+sex:line+batch+behav+TaxaN+seqbatch, data=tempdatasub)
  tempanova<-anova(templm)
  x<-cbind(x,tempanova$`Pr(>F)`[1:7])
}
#The NAs are in the sex:line row because the ANOVA moves the value for that up where the TaxaN was, when there were singularities in the Taxa rel abund
colnames(x)<-Taxa

#Number of taxa that should be removed in removedx
rowSums(is.na(x))

#Filtering out the NA taxa should still result in removing the singularity-taxa.. since the NA result due to that?
removedx<-na.omit(t(x))
removedx<-t(removedx)
#write.table(removedx, "pvalR.txt", sep = "\t")

#Num. sig taxa before correction
uncorrected<-removedx["TaxaN",]<=0.05
table(uncorrected)["TRUE"]
#In this case, 102 significant prior to correction

#Calculate Bonferroni p value by # of taxa left in removedx
#a<-0.05/833

#Subset and extract only row 5 (p values for the Taxa)
removedxtaxa<-removedx[c(5), , drop=FALSE]
taxacol<-colnames(removedxtaxa)
library(magrittr)
#ncol change
fdrcorrection<-(matrix(p.adjust(as.vector(as.matrix(removedxtaxa)), method='fdr'),ncol=833))
#class(fdrcorrection)
#class(subsetsigcol)
#class(subsetsig)
dfr<-as.data.frame(t(fdrcorrection))
dtaxa<-as.list(t(taxacol))

dfr$taxa<-dtaxa
View(dfr)

#Add Taxa by taxa name
dfr$name<-taxakey$Name[match(dfr$taxa,taxakey$ID)]

#How many taxa have a p value less than or equal to 0.05?
tm<-dfr[,"V1"]<=0.05
table(tm)["TRUE"]
#print(tm)
#Which taxa are those?
tempsigdata<-which(dfr[,"V1"] <=0.05, arr.ind=TRUE)
#list of row numbers where this is true
#print(tempsigdata)
#returns row numbers as well as Taxa numbers

#Subset dfr by including only rows in tempsigdata
dfrsub<-dfr[tempsigdata, , drop=FALSE]

write.table(dfr, "full-analysis-fdr-corrected.txt", sep="\t")