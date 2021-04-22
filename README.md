# Parallel Analysis

R function aimed to identify the number of components/factors starting by a raw data matrix by using the principles of Horn's Parallel Analysis.

## Usage

```
parallel(x, iter = 1000, ordinal = FALSE, method = c("perm", "random"),
         alpha = 0.05, standard = FALSE, plot = TRUE, fn = eigen, ...)
```

+ `x` Data frame or matrix of raw data.

+ `iter` Number of iterations.

+ `ordinal` If TRUE, the function uses the polychoric/tetrachoric correlation instead of the Pearson's index (very slow).}

+ `method` Method for random data generation. When method is `perm`, random permutations of observed data are used (permutation is performed within each column independently). When method is `random`, random data normally distributed are generated.

+ `alpha` Alpha threshold to select the number of components/factors.

+ `standard` If TRUE, the analysis is performed on standardized data.

+ `plot` Shows the scree plot overlapping the observed eigen values.

+ `fn` Function to calculate eigenvalues. The default `eigen` uses the principal component analysis with eigen decomposition, while `psych::fa` calculates eigenvalues according to factor analysis.

+ `...` Further arguments for the function `fn`.

## Output

The function returns an object of S3 class `parallel` listing elements:

+ `correlation` type of correlation index used.

+ `method` method of data generation (random or permutations).

+ `synthetic.eigen` matrix of eigen values from iterations of parallel analysis.

+ `pca.eigen` estimated eigenvalues from observed data.

+ `parallel.CI` averages and confidence intervals of eigen values estimated by parallel analysis.

+ `parallel.quantiles` quantiles of eigenvalues estimated by parallel analysis.

+ `suggest.ncomp` suggested number of components/factors.

The `plot` method returns the scree plot overlapping parallel quantiles (gray points and bars).

## Dependencies

The basic functionalities are stand-alone, but the package [`psych`](https://cran.r-project.org/web/packages/psych/) is required to run analyses for ordinal data.

## References

Buja A., Eyuboglu N. (1992). Remarks on parallel analysis. Multivariate Behavioral Research, 27(4), 509-540.

Crawford A.V., Green S.B., Levy R., Lo W.J. Scott L., Svetina D., Thompson M.S. (2010). Evaluation of parallel analysis methods for determining the number of factors. Educational and Psychological Measurement, 70(6), 885-901.

Horn J.L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30, 179-185.

Franklin S.B., Gibson D.J., Robertson P.A., Pohlmann J.T., Fralish, J.S. (1995). Parallel analysis: a method for determining significant principal components. Journal of Vegetation Science, 6(1), 99-106.

Peres-Neto P.R., Jackson D.A., Somers K.M. (2005). How many principal components? Stopping rules for determining the number of non-trivial axes revisited. Computational Statistics and Data Analysis, 49, 974-997.

Weng L.J., Cheng C.P. (2005). Parallel analysis with unidimensional binary data. Educational and Psychological Measurement, 65(5), 697-716.
