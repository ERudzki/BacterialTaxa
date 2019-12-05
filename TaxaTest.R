#-----------
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

#-----------
#Evaluating the bacterial taxa whose relative abundance on the fish skin significantly correlate with infection severity

#install.packages("lookupTable")
library(lookupTable)
#install.packages("plyr")
library(plyr)

data=read.csv("feature-table-6.csv")
#Remove controls as log(0)=Inf
datasub<-subset(data, auc_continuous!=0)
datasubpre<-subset(datasub, timepoint == "pre")
datasublate<-subset(datasub, timepoint =="late")
taxakey=read.csv("Level6TaxaMetadata1.csv")


#lookup do.call
#---------------
#Store the column names of all the taxa into a variable (-1-18)
Taxa<-colnames(datasub[-c(1:26)])
#things tested variable (list)
#thingsTested<-c("sex","line","sex:line","batch","behav","TaxaN")
#I don't know why thingsTested was not working in the below command
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
#Have to organize the row names by the order the ANOVA outputs, otherwise values will not match their name
#The NAs are in the sex:line row because the ANOVA pops the value for that up where the TaxaN was
colnames(x)<-Taxa
#Filtering out the NA taxa should still result in removing the singularity-taxa.. since the NA result due to that?
removedx<-na.omit(t(x))
#removedx<-t(removedx)
#write.table(removedx, "pvalR.txt", sep = "\t")

#BonFeroni correction - most conservative pg 457
#divide all p values by # of tests
#lose a lot of power
#Look at this leaflet for FDR as well, and look up that paper

#Calculate Bonferroni p value by # of taxa left in removedx
a<-0.05/833

#How many taxa have a p value less than or equal to a?
tm<-removedx[,"TaxaN"]<=a
table(tm)["TRUE"]
print(tm)
#Which taxa are those?
sigdata<-which(removedx[,"TaxaN"] <=a, arr.ind=TRUE)
#list of row numbers where this is true
print(sigdata)
#returns row numbers as well as Taxa numbers
#In this case, Taxa952 from row 771

#Subset removedx by including only rows in sigdata
tempsubsig<-removedx[sigdata, , drop=FALSE]
#print(tempsubsig)
subsig<-t(tempsubsig[c(),])
#print(subsig)

#Rename Taxa by taxa name
subsig2<-data.frame(subsig)
namedsig<-setnames(subsig2, as.character(taxakey$ID), as.character(taxakey$Name), skip_absent = TRUE)

write.table(namedsig, "namedsig2.txt", sep="\t")
#-----------
#PRE ONLY
#Store the column names of all the taxa into a variable (-1-18)
Taxa<-colnames(datasubpre[-c(1:26)])
#things tested variable (list)
#thingsTested<-c("sex","line","sex:line","batch","behav","TaxaN")
#I don't know why thingsTested was not working in the below command
xp<-data.frame(row.names= c("sex","line","batch","behav","TaxaN","seqbatch","sex:line"))

for(i in Taxa){
  tempcolindex<-as.numeric(substr(i,5,nchar(i)))
  tempcolindex<- tempcolindex+26
  tempdatasub<-subset(datasubpre, select = c(1:26,tempcolindex))
  names(tempdatasub)[27]<-"TaxaN"
  templm<-lm(log(auc_continuous)~sex+line+sex:line+batch+behav+TaxaN+seqbatch, data=tempdatasub)
  tempanova<-anova(templm)
  xp<-cbind(xp,tempanova$`Pr(>F)`[1:7])
}
#Have to organize the row names by the order the ANOVA outputs, otherwise values will not match their name
#The NAs are in the sex:line row because the ANOVA pops the value for that up where the TaxaN was
colnames(xp)<-Taxa
#Filtering out the NA taxa should still result in removing the singularity-taxa.. since the NA result due to that?
removedxp<-na.omit(t(xp))
#removedx<-t(removedx)
#write.table(removedx, "pvalR.txt", sep = "\t")

#BonFeroni correction - most conservative pg 457
#divide all p values by # of tests
#lose a lot of power
#Look at this leaflet for FDR as well, and look up that paper

#Calculate Bonferroni p value by # of taxa left in removedx
ap<-0.05/587

#How many taxa have a p value less than or equal to a?
tmp<-removedxp[,"TaxaN"]<=0.05
table(tmp)["TRUE"]
print(tmp)
#Which taxa are those?
sigdatap<-which(removedxp[,"TaxaN"] <=0.05, arr.ind=TRUE)
#list of row numbers where this is true
print(sigdatap)
#returns row numbers as well as Taxa numbers

#Subset removedx by including only rows in sigdata
tempsubsigp<-removedxp[sigdatap, , drop=FALSE]
#print(tempsubsig)
subsigp<-t(tempsubsigp)
#print(subsig)

psubsetsig<-subsigp[c(5), , drop=FALSE]
psubsetsigcol<-colnames(psubsetsig)
#fdrcorrectionp<-p.adjust(psubsetsig, method = "fdr", n=length(psubsetsig))
library(magrittr)
fdrcorrectionp<-(matrix(p.adjust(as.vector(as.matrix(psubsetsig)), method='fdr'),ncol=62))
colnames(fdrcorrectionp)<-psubsetsigcol

class(fdrcorrectionp)
class(psubsetsigcol)
class(psubsetsig)

dfr<-as.data.frame(t(fdrcorrectionp))
dtaxa<-as.list(t(psubsetsigcol))

dfr$taxa<-dtaxa

View(dfr)
??transpose

help(colnames)
#Rename Taxa by taxa name
subsig2p<-data.frame(subsigp)
namedsigp<-setnames(subsig2p, as.character(taxakey$ID), as.character(taxakey$Name), skip_absent = TRUE)

write.table(namedsigp, "namedsig2.txt", sep="\t")
#----------------
#make a temp column index (as.numeric) substr i, 5, nchar(i)
#then add 18 to the temp column index (this will tell it to LOOK at that column for that taxa data)
#build a temp data frame (datasub) including c 1:18, plus the column I care about tempcolumnindex
#Rename column 19 to be TaxaN (names(sub)[])

#then your lm - TaxaN referencing datasub
#temp anova (lm)

#into matrix (x) cbind (x, Pr, 1:6 or however many rows)
#END LOOP

#rename colnames(x) to Taxa

#to remove the NAs, transpose the dataframe, and then use ,complete.cases(t(x))



#--------
#Notes for Elizabeth
#SMI columns have NA values, need to fix, not included in model for now
#need to add in seqbatch and make sure Im using both in this dataset
#changed behav to y/n and not 0/1
#changed taxa names to IDs, made a metadata key sheet for this
#--------
#Generalized Linear Model testing the following variables for a singular Taxa (Taxa1)
lm1<-lm(log(auc)~sex+line+sex:line+batch+behav+Taxa1, data=datasub)
summary(lm1)
anova(lm1)
plot(lm1)
hist(resid(lm1))

#Purpose: Looking at any relationship between the log transformed area under the curve of infection (log(AUC)), which is our "infection severity"
#and the relative abundance of Taxa1
#Taking into account effects of sex, line, sex:line, batch, and behavior
#---------         
#How to automate? Ideally I want it to loop through the taxa so that I am running a single command and it does the rest.
#PROBLEM 1
#Want to make each loop automatically add one to the Taxa ID (ex Taxa1 -> Taxa2 -> Taxa3)
#OR tell each loop to move one column to the right (to the next Taxa column)
#PROBLEM 2
#Want to tell each loop to fill in the next column of the results matrix
#Want to make a vector of 586 columns and 6 rows and record anova p-values for each taxa? (Each taxa = 1 column)
#I do not want it overwriting the results of the prior loop

#Can I subset to remove singularities? (where only one sample has any abundance for that taxa), That would be far fewer comparisons.