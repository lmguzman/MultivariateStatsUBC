library(MASS)
library(vegan)
library(cluster)

##########################################
######## CH.11 CANONICAL ANALYSIS ########
data(varechem)
data(varespec)
?varechem 

View(varechem) #site environmental data
View(varespec) #community data


#### 1. RDA ####

# 1.1 Normal RDA
?rda

#First, let's do a simple PCA (unconstrained ordination, ch.9)
my.pca <- rda(varespec, data = varechem)
my.pca

plot(my.pca, type="t") #To see the species axes

#Now, let's do the RDA (constrained ordination)
my.rda <- rda(varespec ~ Ca + Fe + N + S + Baresoil, data = varechem)
my.rda

plot(my.rda, type="t") #To see the species axes


# Test significance of axes:
  # Kaiser-Guttman criterion (Axes above the mean)
my.rda$CCA$eig[my.rda$CCA$eig > mean(my.rda$CCA$eig)]

  # OR do an ANOVA of each ordination axis
?anova.cca
anova.cca(my.rda, by = "axis")

# Proportion of variation explained by each ordination axis 
summ.rda <- summary(my.rda)
summ.rda$cont$importance

# Test significance of explanatory variables (sequentially!)
anova.cca(my.rda, by = "terms")
anova.cca(my.rda, by = "margin") #Not sequentially


# 1.2 Partial RDA

# Looking at the effect of one variable (maintaining others constant)
my.part.rda <- rda(varespec ~ N + Condition(Fe) + Condition(Baresoil) + 
                Condition(S) + Condition(Ca), data = varechem)
my.part.rda



#### 2. CCA ####

# 2.1 Normal CCA

#First, let's do a simple CA (unconstrained ordination, ch.9)
my.ca <- cca(varespec, data = varechem)
my.ca

plot(my.ca, type="t") #To see the species axes

#Now, let's do the CCA (constrained ordination)
my.cca <- cca(varespec ~ Ca + Fe + N + S + Baresoil, data = varechem)
my.cca

plot(my.cca, type="t") #To see the species axes


# Test significance of axes:
# Kaiser-Guttman criterion (Axes above the mean)
my.cca$CCA$eig[my.cca$CCA$eig > mean(my.cca$CCA$eig)]

# OR do an ANOVA of each ordination axis
anova.cca(my.cca, by = "axis")

# Proportion of variation explained by each ordination axis 
summ.cca <- summary(my.cca)
summ.cca$cont$importance

# Test significance of explanatory variables (sequentially)
anova.cca(my.cca, by = "terms")


# 2.2 Partial CCA

# Looking at the effect of one variable (maintaining others constant)
my.part.cca <- cca(varespec ~ Baresoil + Condition(Fe) + Condition(N) + 
                     Condition(S) + Condition(Ca), data = varechem)
my.part.cca






#### 3. Linear Discriminant Analysis ####
# Note: assumes same variation for each group (dispersal on ordination)

# Create an a priori group, but could be known group (hypothesis testing)
#Partition sites according to their communities (Partitioning Around Medoids, ch.8)
varespec.pam <- pam(vegdist(varespec, method = "bray"), 4) 
group <- varespec.pam$clustering
group

#calculate the dispersion of each group
bdist <- betadisper(vegdist(varespec, method = "bray"), group = group)
TukeyHSD(bdist)


# plot sites on a basic CA ordination
varespec.ca <- cca(varechem, scale = TRUE)
plot(varespec.ca, type = "n")
plot(Ca ~ pH, type = "n", data = varechem)
points(Ca ~ pH, col = group, data = varechem)


# Perform LDA determined by TWO variables in varechem
varechem.lda <- lda(group ~ Ca + pH, data = varechem)

# Perform LDA determined by ALL variables in varechem
varechem.lda <- lda(group ~ ., data = varechem)

# A posteriori prediction: How well do the explanatory variables of your lda 
# predict which groups sites belong to?
varechem.lda.p <- predict(varechem.lda)$class

  # Table of matching predictions. Numbers in diagonal are perfect match
(varechem.ldatable <- table(varechem.lda.p, group)) 

  # Proportion of good matches for each group
diag(prop.table(varechem.ldatable,1))




#### 4. Canonical Correlation Analysis (CCorA) ####
?CCorA

Y <- data.frame(varechem$Mo, varechem$Fe, varechem$Al, varechem$pH)
X <- data.frame(varechem$P, varechem$Ca, varechem$Mg, varechem$K, varechem$Humdepth)

CCorA(Y, X)

# ... enjoy!







