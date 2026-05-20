
# ARTPpval

ARTPpval provides fast and accurate p-value calculation methods
for rank truncated product (RTP) and adaptive rank truncated product
(ARTP) statistics, including integral approximation and
cross-entropy importance sampling approaches.

## Installation

You can install the development version of ARTPpval like so:

``` r
install.packages("remotes") 
remotes::install_github("Sarahyyw/ARTPpval")
```

## Example

```r
library(ARTPpval)
p <- c(0.01, 0.02, 0.05, 0.5, , 0.53, 0.7, 0.9, 0.92)
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
```


