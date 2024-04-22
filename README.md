### Multivariate Kolmogorov-Smirnov Test (Fasano-Franceschini)
The `ff_test_2sample()` function in `fasano_franceschini_test.py` implements an $N$-dimensional generalization of the two-sample Kolmogorov-Smirnov test. The test evaluates the null hypothesis $H_0$ that two i.i.d. random
samples $S_1$ and $S_2$ are drawn from the same underlying probability distribution $F_1 = F_2$, against the alternative hypothesis $H_1$ that the samples are drawn from different probability distributions $F_1 \neq F_2$. 
The implementation is based on the method for two and three-dimensions first outlined by Peacock (1983) and expanded upon by Fasano and Franceschini (1987). 
We extend their method of computing the two-sample test statistic $D_n$ to arbitrary dimension here. In their original paper, Fasano and Franceschini (1987) do not attempt any significance testing, 
instead using Monte Carlo sampling to estimate critical values of the test statistic for two- and three-dimensional distributions. We adopt the permutation approach for significance testing of Puritz et al. (2023) in their 
`R` implementation of the same multivariate K-S test outlined by Fasano and Franceschini (1987).

### References
1. Peacock, J. A. (1983). Two-dimensional goodness-of-fit testing in astronomy. Monthly Notices of the Royal Astronomical Society, 202(3), 615-627
2. Fasano, G., & Franceschini, A. (1987). A multidimensional version of the Kolmogorov–Smirnov test. Monthly Notices of the Royal Astronomical Society, 225(1), 155-170.
3. Puritz, C., Ness-Cohn, E. & Braun, R. (2023). fasano.franceschini.test: An Implementation of a Multivariate KS Test in R. The R Journal, 15(3), 159-171.
