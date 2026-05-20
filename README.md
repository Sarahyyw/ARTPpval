
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

This is a basic example:

``` r
library(ARTPpval) 
p <- c(0.01, 0.02, 0.5, 0.7, 0.9, 0.9) 

# RTP integral approximation 
rtp_integral(p.values = p, i = 3) 

# RTP importance sampling 
rtp_isce(p.val = p, J = 3)
```

