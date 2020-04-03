#Principial Component Analysis (PCA) and Factor Analysis in R

#Libraries
library(Hmisc) #Describe Function
library(psych) #Multiple Functions for Statistics and Multivariate Analysis
library(GGally) #ggpairs Function
library(ggplot2) #ggplot2 Functions
library(vioplot) #Violin Plot Function
library(corrplot) #Plot Correlations
library(REdaS) #Bartlett's Test of Sphericity
library(psych) #PCA/FA functions
library(factoextra) #PCA Visualizations
library("FactoMineR") #PCA functions
library(ade4) #PCA Visualizations


#set working directory
setwd("C:\\Users\\eniet\\Documents\\School\\Current\\DSC424\\Project\\")

PCA_Plot = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = data.frame(pcaData$rotation, .names = row.names(pcaData$rotation))
  p + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
}

PCA_Plot_Secondary = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = data.frame(pcaData$rotation, .names = row.names(pcaData$rotation))
  p + geom_text(data=loadings, mapping=aes(x = PC3, y = PC4, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC3", y = "PC4")
}

PCA_Plot_Psyc = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = as.data.frame(unclass(pcaData$loadings))
  s = rep(0, ncol(loadings))
  for (i in 1:ncol(loadings))
  {
    s[i] = 0
    for (j in 1:nrow(loadings))
      s[i] = s[i] + loadings[j, i]^2
    s[i] = sqrt(s[i])
  }
  
  for (i in 1:ncol(loadings))
    loadings[, i] = loadings[, i] / s[i]
  
  loadings$.names = row.names(loadings)
  
  p + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
}

PCA_Plot_Psyc_Secondary = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = as.data.frame(unclass(pcaData$loadings))
  s = rep(0, ncol(loadings))
  for (i in 1:ncol(loadings))
  {
    s[i] = 0
    for (j in 1:nrow(loadings))
      s[i] = s[i] + loadings[j, i]^2
    s[i] = sqrt(s[i])
  }
  
  for (i in 1:ncol(loadings))
    loadings[, i] = loadings[, i] / s[i]
  
  loadings$.names = row.names(loadings)
  
  print(loadings)
  p + geom_text(data=loadings, mapping=aes(x = PC3, y = PC4, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC3", y = "PC4")
}

############################################
###########   Import Data & EDA   ##########


data <- read.csv('data.csv', header = TRUE, sep = ',')

#Check Sample Size and Number of Variables
dim(data)
#569-Sample Size and 33 variables

#Show Structure of Dataset; shows all variables in dataset
str(data, list.len=ncol(data))

#column names
names(data)

#Preview Data 
head(data)

#remove id and assignment and preview
data2 <-  data[,3:32]
head(data2)

#Show descriptive statistics
library(Hmisc)
describe(data2)

#Check for missing values
sum(is.na(data))

#Exploratory Analysis Graphing: Histograms, Correlations
#histogram
library(ggplot2)
ggplot2

ggplot(data = data2,aes(data2$radius_mean)) + geom_histogram()

# To save the ggplot as png
ggsave("p1.png")


#Check Correlation Plot
library(corrplot)
c = cor(data2)
corrplot(c, method = 'ellipse', order = 'AOE')

# Run a correlation test to see how correlated the variables are.  Which correlations are significant
options("scipen"=100, "digits"=5)
round(cor(data2), 2)
MCorrTest = corr.test(data2, adjust="none")
MCorrTest

M = MCorrTest$p
M

# Now, for each element, see if it is < .01 (or whatever significance) and set the entry to 
# true = significant or false
MTest = ifelse(M < .01, T, F)
MTest

# Now lets see how many significant correlations there are for each variable.  We can do
# this by summing the columns of the matrix
colSums(MTest) - 1  # We have to subtract 1 for the diagonal elements (self-correlation)


#all variables are signigicatly correlated to 17-29 other variables
####################################################################
#PCA

#test for factorability

#KMO sampling adequacy
library(psych)
KMO(data2)
#Overall MSA =  0.83 this is goood

#Test Bartlett's Test of Sphericity
library(REdaS)
bart_spher(data2)
#p-value < 2.22e-16 (Very Small Number) meaning, 

#Test for Reliability Analysis using Cronbach's Alpha # check last 10 min of lecture. thatis where she changed this code
library(psych)
alpha(data2,check.keys=TRUE)
#raw_alpha = 0.59

####################################################################

#create PCA and normalize

dataPCA = prcomp(data2, scale = T)
summary(dataPCA)

#scree plot
plot(dataPCA, main= 'Scree Plot - Kaiser Meyer Method')
abline(1, 0)

########################################################

#Check PCA visualizations
plot(dataPCA) #Scree Plot
PCA_Plot(dataPCA) #PCA_plot1
PCA_Plot_Secondary(dataPCA) #PCA_Plot2
biplot(p2) #Biplot

#########################################################


#Best Way to Conduct PCA Analysis

p2 = psych::principal(data2, rotate="varimax", nfactors=3, scores=TRUE)
p2
print(p2$loadings, cutoff=.653, sort=T)

#initial cutoff value of .568 and n = 6 comps with 2 vars


#PCAs Other Available Information

ls(p2)

p2$values
p2$communality
p2$rot.mat

#scores of PCA 
#these are the new values representing the original data
#but instead of having a avlue for each feature it has a value for each set of 
#correlated features aka each component
scores <- p2$scores
dim(scores) #569 3
head(scores)

#####   merge eigenvalues with class    ######
#get diagnosis cols to merge w/ scores
diagnosis <- data[,2]
head(diagnosis)
#M = 2 and B = 1

#bind scores with diagnosis and review
all<- cbind(scores, diagnosis)
dim(all)
head(all)

###### Visualize #########

library(scatterplot3d)
with(all, {
  scatterplot3d(RC1,   # x axis
                RC2,     # y axis
                RC3,    # z axis
                main="3-D Scatterplot Example 1")
})



#scores_1 <- scores[,1]

#min_score <- min(scores_1)
#min_score

#max_score <- max(scores_1)
#max_score

#scores_2 <- scores[,2]
#scores_3 <- scores[,3]
#scores_4 <- scores[,4]
#scores_5 <- scores[,5]
#scores_6 <- scores[,6]


#dont run this yet
round(dataPCA$rotation[,1:3], 2)