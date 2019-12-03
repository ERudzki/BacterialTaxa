data=read.csv("TaxaTest.csv")
#Remove controls as log(0)=Inf
datasub<-subset(data, auc!=0)
taxakey=read.csv("Level6TaxaMetadata.csv")

#lookup do.call

#Store the column names of all the taxa into a variable (-1-18)
Taxa<-colnames(datasub[-c(1:18)])
#things tested variable (list)
#thingsTested<-c("sex","line","sex:line","batch","behav","TaxaN")
#I don't know why thingsTested was not working in the below command
x<-data.frame(row.names= c("sex","line","batch","behav","TaxaN","sex:line"))

for(i in Taxa){
  tempcolindex<-as.numeric(substr(i,5,nchar(i)))
  tempcolindex<- tempcolindex+18
  tempdatasub<-subset(datasub, select = c(1:18,tempcolindex))
  names(tempdatasub)[19]<-"TaxaN"
  templm<-lm(log(auc)~sex+line+sex:line+batch+behav+TaxaN, data=tempdatasub)
  tempanova<-anova(templm)
  x<-cbind(x,tempanova$`Pr(>F)`[1:6])
}
#Have to organize the row names by the order the ANOVA outputs, otherwise values will not match their name
#The NAs are in the sex:line row because the ANOVA pops the value for that up where the TaxaN was
colnames(x)<-Taxa
#Filtering out the NA taxa should still result in removing the singularity-taxa.. since the NA result due to that?
removedx<-na.omit(t(x))
#removedx<-t(removedx)

#How many taxa have a p value less than or equal to 0.02?
tm<-removedx[,"TaxaN"]<=0.05
table(tm)["TRUE"]

#Which taxa are those?
sigdata<-which(removedx[,"TaxaN"] <=0.05, arr.ind=TRUE)
#list of row numbers where this is true
print(sigdata)
#returns row numbers as well as Taxa numbers

#Subset removedx by including only rows in sigdata
tempsubfinal<-removedx[sigdata,]
subfinal<-t(tempsubfinal)

help(which)
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