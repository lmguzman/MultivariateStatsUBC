library(vegan)
library(ade4)
library(gclus)
library(FD)
library(dplyr)

source("coldiss.R")

##########################################################
######## COEFFICIENTS OF ECOLOGICAL RESEMBLANCE ##########

#SYMMETRICAL COEFFICIENTS:
#In symmetrical coefficients, the state zero for two objects is treated in exactly the
#same way as any other pair of values, when computing a similarity. 
#These coefficients can be used in cases where the state "zero" is a valid basis for 
#comparing two objects and represents the same kind of information as any other value.

#ASYMMETRICAL COEFFICIENTS:
#NOT drawing any ecological conclusion from the absence of a species at two sites. 
#In numerical terms, this means to skip double zeros altogether when computing similarity 
#or distance coefficients using species presence/absence or abundance data. 
#Coefficients of this type are called asymmetrical because they treat zeros in a different
#way than other values.

#aviurba data set
data(aviurba)
?aviurba

aviurba$fau #community data
aviurba$traits #trait data
View(aviurba$mil) #site environmental data

# Create binary environmental dataset
aviurba.env <- aviurba$mil %>%
  select(-veg.cover)
aviurba.env <- gsub("no", 0, as.matrix(aviurba.env))
aviurba.env <- gsub("yes", 1, as.matrix(aviurba.env))
aviurba.env <- as.data.frame(apply(aviurba.env, 2, as.numeric))

# Make a NMDS plotting function
nmds_plot <- function(dist.matrix, colorby){
  m <- metaMDS(dist.matrix) #NMDS object
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- rainbow(nlevels(colorby)) #vector of colors needed
  ordiplot(m, display = c("sites"), type = "n")
  points(m, col = cols[colorby], pch = 16)
  text(m, labels = row.names(aviurba.env), cex = 0.8) 
  legend("topright", title ="factor of choice", legend=levels(colorby), col=cols, pch = 16)
}
#########################################
### 1. BINARY DATA (presence/absence) ###
?dist.binary
?vegdist #with "binary = TRUE"

# 1.1 Symmetrical binary coefficients (for environmental data)
# Choose method = 2, 4, 6, 8 or 9

#Simple matching (method = 2)
dist <- dist.binary(aviurba.env, method = 2)
coldiss(dist)
nmds_plot(dist, aviurba.env$fields)


# 1.2 Asymmetrical binary coefficients (for species data)
# Choose method = 1, 3, 5, 7 or 10

#Jaccard
dist <- dist.binary(aviurba$fau, method = 1)
dist <- vegdist(aviurba$fau, method = "jaccard", binary = TRUE) #Same as using dist.binary()
coldiss(dist)
nmds_plot(dist, aviurba.env$fields)


############################
### 2. QUANTITATIVE DATA ###
?vegdist

# 2.1 Symmetrical quantitative coefficients (for environmental data)

#Gower
dist <- vegdist(aviurba.env, method = "gower")
coldiss(dist)
nmds_plot(dist, aviurba.env$fields)

#Gower for mixed types (mix of qual, semi-quant and quant descriptors)
#function gowdis() is better than vegdist() for mixed types
dist <- gowdis(aviurba$mil) #now includes 
coldiss(dist)
nmds_plot(dist, aviurba.env$fields)


# 2.2 Asymmetrical quantitative coefficients (for species data)

#Bray-Curtis, Jaccard, Kulczynski
dist <- vegdist(aviurba$fau, method = "bray") #try other methods
coldiss(dist)
nmds_plot(dist, colorby = aviurba.env$fields)
    #NOTE: Jaccard index is metric, and probably should be preferred instead of 
    #the default Bray-Curtis which is semimetric.



