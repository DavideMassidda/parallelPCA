PARALLEL PRINCIPAL COMPONENT ANALYSIS

Script version: 1.0
Author: Davide Massidda
e-mail: davide.massidda@humandata.it
Date: June 11, 2013
URL: http://www.insular.it, http://www.humandata.it
License: GPLv3
Description: this script provides the function parallelPCA, which performs a parallel
principal component analysis on the correlation matrix from continuous and ordinal data.
Dependencies: psych

USAGE
parallelPCA(x, iter = 1000, ordinal = FALSE, method = "random",
            alpha = 0.05, standard = FALSE, plot = TRUE)

ARGUMENTS
x: row data matrix.
iter: the number of iterations for parallel analysis.
ordinal: specifies if the data matrix contains ordinal data. If ordinal data are
         present, tetrachoric or polycoric correlations are used instead of the
         Pearson's indices. The number of categories for ordinal data is calculated
         as max(x)-(min(x)-1).
method: specifies if the parallel analysis must be performed using random data
        (method="random") or permutations of observed data (method="perm").
alpha: alpha level to calculate quantiles and confidence intervals.
standard: specifies if the data must be standardized.
plot: plots the scree test.

VALUE
The function returns and object of class 'parpca', containing the slots:
@correlation: the type of correlation index used.
@method: the type of data generation (random or permutations).
@synthetic.eigen: the matrix of eigen values from each iteration of parallel analysis.
@pca.eigen: estimated eigenvalues from observed data.
@parallel.CI = averages and confidence intevals of eigenvalues estimated by parallel analysis.
@parallel.quantiles = quantiles of eigenvalues estimated by parallel analysis.
The plot method is implemented on the class 'parpca'. The function plot() applied
to class 'parpca' returns the scree plot overlapping parallel quantiles (gray
points and bars).

REFERENCES
[1] Buja A., Eyuboglu N. (1992). Remarks on parallel analysis. Multivariate Behavioral
      Research, 27(4), 509–540.
[2] Crawford A.V., Green S.B., Levy R., Lo W.J. Scott L., Svetina D., Thompson M.S. (2010)
      Evaluation of parallel analysis methods for determining the number of factors.
      Educational and Psychological Measurement, 70(6), 885–901.
[3] Horn J.L. (1965). A rationale and test for the number of factors in factor analysis.
      Psychometrika, 30, 179–185.
[4] Franklin S.B., Gibson D.J., Robertson P.A., Pohlmann J.T., Fralish, J.S. (1995).
      Parallel analysis: a method for determining significant principal components.
      Journal of Vegetation Science, 6(1), 99-106
[5] Peres-Neto P.R., Jackson D.A., Somers K.M. (2005). How many principal components?
      Stopping rules for determining the number of non-trivial axes revisited.
      Computational Statistics and Data Analysis, 49, 974-997
[6] Weng L.J., Cheng C.P. (2005). Parallel analysis with unidimensional binary data.
      Educational and Psychological Measurement, 65(5), 697-716.

