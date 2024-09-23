
# NMLmulti

<!-- badges: start -->
<!-- badges: end -->

Calculate NMP Penalities as described in Kellen and Klauer (2020).

## Installation

You can install the development version of NMLmulti from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("davidkellen/NMLmulti")
```

## Basic Example

The basic example uses the SCR model discussed in the paper. Note that in contrast to the code examples given in the paper, the function name needs to be passed without the package name in front (i.e., `SCR` instead of `examplemodels::SCR` or `NMLmulti::SCR`).

``` r
library(NMLmulti)
## basic example code
nml_scr <- run_nml(fun= SCR, parl=10, ks=rep(3,5), Ns=c(rep(100,4),200), 
                   fits = 2, batchsize=5000, burn=10000, precision=0.1)
```

(To speed-up this process for testing purposes, you can increase the precision value.)

## Example Using Custom C++ Likelihood Function

Save the C++ code of your likelihood function in a `.cpp` file. We follow the instruction of the paper and save the following code in file `examplemodels.cpp`. Note that in order to make sure we use the new function, we name it `SCR2` (instead of `SCR`).

``` C++
#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
double SCR2(NumericVector Q0, NumericVector chd, NumericVector NN){

  NumericVector Q = pnorm(Q0);
  
 double Vt = Q[0];
 double Gt = Q[1];
 double Vr = Q[2];
 double Gr = Q[3];
 double a  = Q[4];
 double b  = Q[5];
 double Vtx = Q[6];
 double Gtx = Q[7];
 double Vrx = Q[8];
 double Grx = Q[9];

  Rcpp::NumericVector ee(15);

	ee[0] = Vt + (1-Vt)*Gt*a + (1-Vt)*(1-Gt)*b*a;
	ee[1] = (1-Vt)*Gt*(1-a) + (1-Vt)*(1-Gt)*b*(1-a);
	ee[2] = (1-Vt)*(1-Gt)*(1-b);

	ee[3] = (1-Vr)*Gr*a + (1-Vr)*(1-Gr)*b*a;
	ee[4] = Vr + (1-Vr)*Gr*(1-a) + (1-Vr)*(1-Gr)*b*(1-a);
	ee[5] = (1-Vr)*(1-Gr)*(1-b);

	ee[6] = Vtx + (1-Vtx)*Gtx*a + (1-Vtx)*(1-Gtx)*b*a;
	ee[7] = (1-Vtx)*Gtx*(1-a) + (1-Vtx)*(1-Gtx)*b*(1-a);
	ee[8] = (1-Vtx)*(1-Gtx)*(1-b);

	ee[9] = (1-Vrx)*Grx*a + (1-Vrx)*(1-Grx)*b*a;
	ee[10] = Vrx + (1-Vrx)*Grx*(1-a) + (1-Vrx)*(1-Grx)*b*(1-a);
	ee[11] = (1-Vrx)*(1-Grx)*(1-b);

	ee[12] = b*a;
	ee[13] = b*(1-a);
	ee[14] = (1-b);
	
	ee = ee*NN;
	
  double LL = 0.0;
  
  for(int ii=0; ii<15; ii++){
    if(chd[ii] > 0.0){
      LL = LL + chd[ii]*(log(chd[ii])-log(ee[ii]));
    }
  }
  return 2.0*LL;
}
```

We then create and install a package based on this `.cpp` file:

``` r
Rcpp::Rcpp.package.skeleton("examplemodels", attributes = TRUE, 
                            cpp_files = "examplemodels.cpp")
devtools::install("examplemodels")

```

We can then use the likelihood function from this package **as long as we pass the package name in the `packages_multicore` argument**: `packages_multicore = "examplemodels"`. This ensures that the functions in the package are available at each core.

``` r
library("NMLmulti")
library("examplemodels")
nml_scr <- run_nml(fun= SCR2, parl=10, ks=rep(3,5), Ns=c(rep(100,4),200),
                   packages_multicore = "examplemodels",
                   fits = 2, batchsize=5000, burn=10000, precision=0.1)
```

## Known Limitations

At the moment, `run_nml()` does not support `cores = 1`.
