PARALLEL PRINCIPAL COMPONENT ANALYSIS<br>
<br>
Script version: 2.1<br>
Author: Davide Massidda<br>
e-mail: davide.massidda@humandata.it<br>
Date: August 14, 2014<br>
URL: http://www.insular.it, http://www.humandata.it<br>
License: GPLv3<br>
Description: this script provides the function parallelPCA, which performs a parallel
principal component analysis on the correlation matrix from continuous and ordinal data.<br>
Dependencies: psych<br>
<br>
USAGE<br>
parallelPCA(x, iter = 1000, ordinal = FALSE, method = "perm",
            alpha = 0.05, standard = FALSE, plot = TRUE, FUN = eigen, ...)<br>
<br>
ARGUMENTS<br>
x: row data matrix.<br>
iter: the number of iterations for parallel analysis.<br>
ordinal: specifies if the data matrix contains ordinal data. If ordinal data are
         present, tetrachoric or polycoric correlations are used instead of the
         Pearson's indices. The number of categories for ordinal data is calculated
         as max(x)-(min(x)-1).<br>
method: specifies if the parallel analysis must be performed using random data
        (method="random") or permutations of observed data (method="perm").<br>
alpha: alpha level to calculate quantiles and confidence intervals.<br>
standard: specifies if the data must be standardized.<br>
plot: plots the scree test.<br>
FUN: function to calculate eigen values. From package psych, options are: principal
     (for Principal Component Analysis) and fa (for Exploratory Factor Analysis).<br>
...: optional arguments to be passed to FUN. When FUN=eigen, specifying symmetric=TRUE
     and only.values=TRUE can be helpful to reduce the system time.<br>
<br>
VALUE<br>
The function returns and object of class 'parpca', containing the slots:<br>
@correlation: the type of correlation index used.<br>
@method: the type of data generation (random or permutations).<br>
@synthetic.eigen: the matrix of eigen values from each iteration of parallel analysis.<br>
@pca.eigen: estimated eigenvalues from observed data.<br>
@parallel.CI = averages and confidence intevals of eigenvalues estimated by parallel analysis.<br>
@parallel.quantiles = quantiles of eigenvalues estimated by parallel analysis.<br>
The plot method is implemented on the class 'parpca'. The function plot() applied
to class 'parpca' returns the scree plot overlapping parallel quantiles (gray
points and bars).<br>
<br>
REFERENCES<br>
[1] Buja A., Eyuboglu N. (1992). Remarks on parallel analysis. Multivariate Behavioral
      Research, 27(4), 509–540.<br>
[2] Crawford A.V., Green S.B., Levy R., Lo W.J. Scott L., Svetina D., Thompson M.S. (2010)
      Evaluation of parallel analysis methods for determining the number of factors.
      Educational and Psychological Measurement, 70(6), 885–901.<br>
[3] Horn J.L. (1965). A rationale and test for the number of factors in factor analysis.
      Psychometrika, 30, 179–185.<br>
[4] Franklin S.B., Gibson D.J., Robertson P.A., Pohlmann J.T., Fralish, J.S. (1995).
      Parallel analysis: a method for determining significant principal components.
      Journal of Vegetation Science, 6(1), 99-106.<br>
[5] Peres-Neto P.R., Jackson D.A., Somers K.M. (2005). How many principal components?
      Stopping rules for determining the number of non-trivial axes revisited.
      Computational Statistics and Data Analysis, 49, 974-997.<br>
[6] Weng L.J., Cheng C.P. (2005). Parallel analysis with unidimensional binary data.
      Educational and Psychological Measurement, 65(5), 697-716.<br>

