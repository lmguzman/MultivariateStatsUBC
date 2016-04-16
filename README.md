# Multivariate Stats UBC

![](image1.JPG)

## Key Concepts

- **Similarity/ Dissimilarity matrix:** A square matrix of association among objects or descriptors
- **Correlation:** Aims at establishing the interdependence between two random variables
- **Regression:** Fits a straight line through the set of n points in such a way that makes the sum of squared residuals of the model as small as possible
- **Multiple correlation:** Measures the intensity of the relationship between a response variable and a linear combination of several explanatory variables. 
- **Multiple regression:** An extension of simple linear regression. It is used when we want to predict the value of a variable based on the value of two or more other variables. 
- **Partial correlation:** Measures the intensity of the relationship between two variables, while taking into account their relationships to other variables
- **Partial regression:** A way of estimating how much of the variation of the response variable can be attributed exclusively to one set of factors, once the effect of the other set has been taken into account. 
- **Qualitative:** Categorical data. Eg: male vs female
- **Quantitative:** Continuous axis of real numbers. Eg: length
- **Semiquantitative:** Small number of ordered classes. Eg: small, medium, large
- **Clustering:** Partitions a collection of objects. 
	- *Single Linkage:* Type of hierarchical agglomerative clustering. 
	- *Partitioning by K-means:* Determines a partition into groups or clusters determined by the user. 
	- *Hierarchical agglomerative:* Clusters are formed hierarchically, starting from the most similar objects and then relaxing the similarity criterion. 
- **Ordination:** Characterizes the main trends of variation of the objects with respect to all descriptors. Start from the scaling of the objects in full-dimensional space and attempt to represent them in a few dimensions while preserving the distance relationship among the objects. 
	- *PCA:* *Linear* combinations of the original descriptors
	- *PCoA:* This method obtains an euclidean representation from all types of variables, even sets of variables of mixed levels of presicion. Principal coordinates are also functions of the original variables. 
	- *MDS:*The priority is to represent the objects in a small and specified number of dimensions, preserving the order of relationships among objects. It is not limited to euclidean distance matrices. Not an eigenvector method => axes are arbitrary. 
	- *CA:* Uses the $X^2$ distance which is a coefficient that excludes double zeroes. Frequently used on species abundance data. Preserves the $X^2$ distance between the rows or columns of the contingency table. 
- **Q mode analysis:** Association matrix among objects
- **R mode analysis:** Association matriz among descriptors
- **Objects:** Sites, sampling units, subjects defined a priori. 
- **Descriptors:** Also called variable, attributes or characters used to describe or compare the objects. ie: species
- **Descriptive analysis:** Method available for explaining the structure of one or several ecological descriptors. The purpose here is data exploration, not hypothesis testing.
- **Forecasting analysis:** Extend into the future or to different situations, structuralr elationships among descripts that have been quantifed for a given data set. Also called correlative models.
- **Predictive analysis:** When the relationships are assumed to be causal and to describe a process. 
- **Path analysis:** An extension of multiple linear regression which allows the decomposition and interpretation of linear relationships among a small number of descriptors. Used to test alternative hypothesis. 
- **Matrix comparison:** To perform direct comparison when the response descriptors form a multivariable dataset. Uses the similarity matrices.
	- *Mantel test:* Calculates the "correlation coefficient" between two matrices. If you have a similarity matrix for species composition, and another for environmental data for the same sites, or geographic distances among sites. (Matrices are the same size)
	- *ANOSIM test:* (could also use PERMANOVA - Prob more robust?) Test if there are significant differences in community composition between two groups. Analogous to an ANOVA for community composotion. 
	- *Procrustes analysis:* Ordination technique to find a compromise ordination for two data matrices concerning the same objects. Minimizes the SS distances between corresponding sites in the matrix. 
- **4th corner problem:** How do the biological and behavioral characteristics of species determine their relative locations in an ecosystem? - You have 3 matrices: species x site, environment x site and species x behaviour/functional traits. You want to estimate a environment x behaviour/functional traits matrix. How do we associate behavioural/functional traits to habitat characteristics? -> estimates the parameters and test for significance. 
- **Canonical analysis:** Allows to perform a direct comparison of two data matrices. Brings out all the variance in matrix Y that is related to matrix X. <- Ordination of Y under constraint of X. Combines the concepts of ordination and regression. 
  - *RDA:* Direct extension of multiple regression to the modelling of multivariate response data. Can also be seen as an extension of principal component analysis. Preserves the euclidean distance among objects.
  - *CCA:* Similar to RDA but preserved the $X^2$ distance (as in correspondence analysis) among objects.
  - *CCorA:* Direct extension of linear correlation.
  - *Canonical* discriminant analysis: Objects are divited into k groups described by a qualitative descriptor. 
  - *Partial RDA/CCA:* The extension of partial linear regression. Control the effect of matrix W on Y when analysing the effect of Y on a set of variables X. 
- **Time series**
	- Detecting cycles
	- Trend extractions
  - *Correlograms:* The coefficients of autocorrelation (or autocovariance) are plotted as a function of the lag (k)
  - *Peridogram:* Harmonic analysis
  - *Periodic variability:* Spectral analysis- correlograms and peridograms lead to spectral analysis. 
- **Spatial analysis** 
  - *Correlograms:* Coefficients of correlation via Moran's I or Mantel correlograms
  - *Moran's I:* Measure of spatial autocorrelation
  - *Trend surface analysis:* Produces smoothed maps. Express observed values as a polynomial function of the greographic coordinates X and Y of the sapling sites. 
- **Variation partitioning/ Partial canonical analysis:** 3 data sets, Y with response variables (species composition), X with explanatory environmental variables and W with explanatory spatial variables. 
- **dbMEM:** Captures the variation in geographic coordinates via eigenfunctions. Can replace trend surface analysis



#R

CRAN Task View: Analysis of Ecological and Environmental Data:
https://cran.r-project.org/web/views/Environmetrics.html

CRAN Task View: Multivariate Statistics:
https://cran.r-project.org/web/views/Multivariate.html

CRAN Task View: Cluster Analysis & Finite Mixture Models
https://cran.r-project.org/web/views/Cluster.html




