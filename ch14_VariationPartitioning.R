library(vegan)
#library(ade4)

setwd("C:/Users/User/Dropbox/UBC/Multivariate Study Group/NEwR data")

#### VARIATION PARTITIONING ####
mite <- read.table("mite.txt")
mite.env <- read.table("mite_env.txt")
mite.xy <- read.table("mite_xy.txt")

mite.hell<-decostand(mite, "hellinger")

#Centered xy coordinates
mite.xy.c <- scale(mite.xy, center = TRUE, scale = FALSE)

#Figure of variation partitioning (empty)
showvarparts(2)

#Spatial component using polynomials

mite.poly <- poly(as.matrix(mite.xy.c), degree = 3, raw=TRUE) #raw=FALSE is orthogonal
colnames(mite.poly) <- c("X", "X2", "X3", "Y", "XY", "X2Y", "Y2", "XY2", "Y3")

# RDA of all 9 polynomials (unrealistic)
rda.space <- rda(mite.hell ~ ., data = as.data.frame(mite.poly))
rda.space

#Forward selection: Which terms of polynomial to keep?
# evaluates variables according to two criteria: 
# 1. if their inclusion into the model leads to significant increase of explained variance
# 2. if the AIC of the new model is lower than AIC of the more simple model. 

fwd <- ordistep(rda(mite.hell ~ 1, data = as.data.frame(mite.poly)), 
         scope = formula(rda.space), direction = 'forward')

#Retained variables: Y + Y3 + Y2 + X2Y + X + XY

# Partitioning environmental and spatial variables
mm1 <- model.matrix(~ ., mite.env)[,-1]
mm2 <- model.matrix(~ Y + Y3 + Y2 + X2Y + X + XY, mite.poly)[,-1]
variation.part <- varpart(mite.hell, mm1, mm2)
variation.part

indiv.frac <- as.character(round(variation.part$part$indfract[,3], 3))

showvarparts(2, labels = indiv.frac, cex = 1, bg = c("green","hotpink"))

# test fraction [a] (environmental)
rda.result <- rda(mite.hell ~ mm1 + Condition(mm2))
anova(rda.result, step=200, perm.max=200)

# test fraction [c] (spatial)
rda.result <- rda(mite.hell ~ mm2 + Condition(mm1))
anova(rda.result, step=200, perm.max=200)




