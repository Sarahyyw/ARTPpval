
# ARTPpval

ARTPpval provides fast and accurate p-value calculation methods
for rank truncated product (RTP) and adaptive rank truncated product
(ARTP) test, including integral approximation, cross-entropy importance sampling and Ultra-fast Interpolation approaches.

## Installation

You can install the development version of ARTPpval like so:

``` r
install.packages("remotes") 
remotes::install_github("Sarahyyw/ARTPpval")
```

## Example

```r
library(ARTPpval)
p <- c(0.01, 0.02, 0.05, 0.5, 0.53, 0.7, 0.9, 0.92)
p <- c(0.0001, 0.002, 0.005, rep(0.5, 100))
```

### RTP p-value calculation

```r
# RTP integral approximation
rtp_integral(p.values = p, i = 3)

# RTP importance sampling
rtp_isce(p.val = p, J = 3)
```

### ARTP p-value calculation

```r
# ARTP importance sampling
artp_isce(p.values = p)

# ARTP ultra-fast interpolation (UFI)
artp_ufi(p.values = p)
```

The returned list contains:

* `p.value`: estimated ARTP p-value
* `cutpoint`: RTP truncation point selected by ARTP
* `observed_stat`: observed ARTP test statistic (minimum RTP p-value)
* `rtp_pvalues`: RTP p-values across all truncation points


## References

* Yu, K., Li, Q., Bergen, A. W., Pfeiffer, R. M., Rosenberg, P. S.,
  Caporaso, N., Kraft, P., & Chatterjee, N. (2009).
  *Pathway analysis by adaptive combination of P-values.*
  Genetic Epidemiology, 33(8), 700–709.

* Vsevolozhskaya, O. A., Hu, F., & Zaykin, D. V. (2019).
  *Detecting weak signals by combining small P-values in genetic association studies.*
  Frontiers in Genetics, 10, 1051.

* Fang, Y., Chang, C., & Tseng, G. (2022).
  *On p-value combination of independent and frequent signals:
  asymptotic efficiency and Fisher ensemble.*
  arXiv preprint arXiv:2203.11748.










