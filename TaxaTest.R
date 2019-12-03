data=read.csv("TaxaTest.csv")
#Remove controls as log(0)=Inf
datasub<-subset(data, auc!=0)

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