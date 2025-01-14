## fasano-franceschini-test

A small Python package implementing the $N$-dimensional generalization of the two-sample Kolmogorov-Smirnov test described by [Fasano and Franceschini (1987)](https://ui.adsabs.harvard.edu/abs/1987MNRAS.225..155F/abstract). The test evaluates the null hypothesis $H_0$ that two i.i.d. random samples $S_1$ and $S_2$ are drawn from the same underlying probability distribution $F_1 = F_2$, against the alternative hypothesis $H_1$ that the samples are drawn from different probability distributions $F_1 \neq F_2$. The implementation is based on the method for two and three dimensions first outlined by [Peacock (1983)](https://ui.adsabs.harvard.edu/abs/1983MNRAS.202..615P/abstract) and expanded upon by [Fasano and Franceschini (1987)](https://ui.adsabs.harvard.edu/abs/1987MNRAS.225..155F/abstract). We have naturally extended their method of computing the two-sample test statistic $D_n$ to arbitrary dimension here. In their original paper, [Fasano and Franceschini (1987)](https://ui.adsabs.harvard.edu/abs/1987MNRAS.225..155F/abstract) do not attempt any significance testing, instead using Monte Carlo sampling to estimate critical values of the test statistic for two- and three-dimensional distributions. Here we adopt the permutation approach for significance testing of [Puritz et al. (2023)](https://ui.adsabs.harvard.edu/abs/2021arXiv210610539P/abstract) in their `R` implementation of the multivariate K-S test.

If you use this package in your research, please cite [Chow & Brown (2025)](https://doi.org/10.1016/j.icarus.2024.116444) and link to this repository on GitHub or Zenodo.

## Installation:
Clone the GitHub repository and install via setup tools:

```
git clone https://github.com/wmpg/fasano-franceschini-test.git
cd fasano-franceschini-test
python setup.py install
```

or simply download `fftest/fftest.py` and import it into your Python file.

## Requirements:
[`numpy`](https://numpy.org/)

## Example
To perform the K-S test for two samples `S1` and `S2`:

```
from fftest import fftest_2samp
Dn, pval = fftest_2samp(s1, s2, n_perms, threads, seed)
```

where `Dn` is the two-sample test statistic and `pval` is the p-value computed using the permutation test approach of [Puritz et al. (2023)](https://ui.adsabs.harvard.edu/abs/2021arXiv210610539P/abstract). The input parameters are defined as follows:

```
Parameters
------------
s1 : array-like
    n1 x d array of samples with length n1 and dimension d
s2: array-like
    n2 x d array of samples with length n2 and dimension d
n_perms: int
    The number of permutations to use for the permutation test for significance testing. If set to 0, only the test statistic Dn is returned (default: 100)
threads: int
    The number of threads to use for the permutation testing. If set to "auto", the number of threads is determined by multiprocessing.cpu_count() - 1 (default: 1)
seed: NumPy seed
    Value to seed the RNG for the permutation testing to reproducibly compute p-values (default: None)
```

See `docs` for further examples.

## References
1. Peacock, J. A. (1983). Two-dimensional goodness-of-fit testing in astronomy. Monthly Notices of the Royal Astronomical Society, 202(3), 615-627. http://dx.doi.org/10.1093/mnras/202.3.615.
2. Fasano, G., & Franceschini, A. (1987). A multidimensional version of the Kolmogorov–Smirnov test. Monthly Notices of the Royal Astronomical Society, 225(1), 155-170. http://dx.doi.org/10.1093/mnras/225.1.155.
3. Puritz, C., Ness-Cohn, E. & Braun, R. (2023). fasano.franceschini.test: An Implementation of a Multivariate KS Test in R. The R Journal, 15(3), 159-171. http://dx.doi.org/10.32614/RJ-2023-067.
