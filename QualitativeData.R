#####################################################
# Numerical Ecology, Legendre & Legendre, UBC group
# Qualitative data - Practice code for: 
# Contingency tables, entropy and hierarchical models
# 20 November 2015
#####################################################

# Inspiration from:
# http://www.statmethods.net/stats/frequencies.html
# http://alumni.media.mit.edu/~tpminka/courses/36-350.2001/lectures/day12/
# http://stats.stackexchange.com/questions/14158/how-to-generate-random-categorical-data

####################################################
##########    CONTINGENCY TABLES         ###########

### 2-way table
Money <- factor(c("Broken","Broken","Ok","Broken","Ok","Ok"))
Edu <- factor(c("PhD","PhD","PhD","MSc","MSc","MSc"))
ct <- table(Money,Edu)
ct 
t (ct) #Transpose
mosaicplot(ct) # Mosaicplot

### Silly data 
x <- c(rep("Ok",0.3*500),rep("Broken",0.7*500))
y <- c(rep("PhD",0.4*500),rep("MSc",0.6*500))
Money <- sample(x) 
Edu <- sample(y) 
mytable <- table(Money,Edu)
mytable

margin.table(mytable, 1) # x frequencies (summed over y) 
margin.table(mytable, 2) # y frequencies (summed over x)

prop.table(mytable, 1) # row percentages (Probability of Broken or ok)
prop.table(mytable, 2) # column percentages (Probability of PhD or MSc)
prop.table(mytable) # cell percentages (Probability of broken/ok being PhD or MSc))

### 3-Way Frequency Table 
z <- c(rep("happy",0.6*500),rep("unhappy",0.4*500))
Happi <- sample(z) 

mytable2 <- ftable(Happi,Money,Edu) #3 way table 
mytable2

mosaicplot(mytable2,shade=T) #Standarized resioduals

install.packages('vcd')
library (vcd)
cotabplot(mytable2, panel = cotab_coindep, shade = TRUE, legend = FALSE, type = "assoc")


### Chi-square, fisher, Mantelhaen, correspondence
chisq.test(mytable2)  # Ho: Descriptors are independent
fisher.test(mytable2) #
mantelhaen.test(Happi,Money,Edu) #CM Chi-square of interaction 

####################################################
##########       ENTROPY                 ###########

# From http://stackoverflow.com/questions/27254550/calculating-entropy
info <- function(CLASS.FREQ){
  freq.class <- CLASS.FREQ
  info <- 0
  for(i in 1:length(freq.class)){
    if(freq.class[[i]] != 0){ # zero check in class
      entropy <- -sum(freq.class[[i]] * log2(freq.class[[i]]))  #I calculate the entropy for each class  here
    }else{ 
      entropy <- 0
    } 
    info <- info + entropy # sum up entropy from all classes
  }
  return(info)
}

freqs <- table(Edu)/length(Edu)
freqs2 <- table(Money)/length(Money) 
freqs3 <- table(Happi)/length(Happi)

#Calculate entropy: 
info(freqs)  # (Bits)
info(freqs2) 
-sum(freqs * log2(freqs)) 
entropy.empirical(freqs, unit="log2") 

#With package entropy
install.packages("entropy") # http: //strimmerlab.org/software/entropy/
library (entropy)
help (entropy)

freqs.empirical(freqs) #Edu
entropy(freqs, method="ML") # Also "MM", "Jeffreys", "Laplace", "SG", "minimax", "CS", "NSB", "shrink"
entropy.empirical(freqs, unit=c("log"))  #Nats
entropy.empirical(freqs, unit=c("log2"))  #Bits
entropy.empirical(freqs, unit=c("log10")) #Logits

# Kullback-Leiber (KL) divergence 
KL.plugin (freqs, freqs2)  #from Happi to Money. 
KL.plugin (freqs2, freqs3)  #from Money to Edu
KL.plugin (freqs, freqs3)  #from Happi to Edu

####################################################
##########    LOGLINEAR - HIERARCHICAL   ###########

library(MASS)
mytable <- xtabs(~ Happi + Money + Edu)  #3-way contingency table 

# Mutual Independence: 
loglm(~ Happi + Money + Edu, mytable) #Ho: Pairwise independent.

# Partial Independence:
loglm(~ Happi + Money + Edu + Money * Edu, mytable)  
# Ho: Happiness is partially independent of Money and Education
# (i.e., Happi is independent of the composite variable MoneyEdu).

#Conditional Independence: 
loglm(~Happi+Money+Edu+ Happi*Edu + Money*Edu, mytable) 
#Ho: Happiness is independent of Money, given Edu.

# Ho: No Three-Way Interaction
loglm(~Happi+Money+Edu+Happi*Money+Happi*Edu+Money*Edu, mytable)





##### diversity
install.packages('vegan')
library(vegan)

data(dune)
data("dune.env")
?dune

barplot(sort(colSums(dune)))


s <- specnumber(dune)

rar <- rarefy(dune, min(rowSums(dune)))
hist(rar, breaks = 15)
barplot(sort(rar))


spac <- specaccum(dune)
plot(spac, ci.type = 'polygon', ci.col = 'green')
spac.rar <- specaccum(dune, method = 'rarefaction')
plot(spac.rar, ci.type = 'polygon', ci.col = 'yellow',add =TRUE)


