![OpSel](https://github.com/makoto-yamashita/OpSel.jl/blob/logo-registration/logo/opsel-logo.png)

# OpSel.jl
Efficient optimal selection for tree breeding

Optimal selection problems are to find optimum of selection of genotypes that maximizes genetic gains under a constraint on genetic diversity which involves Wright's numerator relation ship.

Using a conic programming approach, this package provides efficient numerical methods for optimal selection problems arising from tree breeding. This package now implements two methods
1. the compact SOCP formulation for unequally deployment problem
2. the steepest-ascent method for equally deployment problem

## Installation 

## Usage

## Basic Formulation

An unequally-type of optimization problem is of form:

max g'*x subject to x'*A*x <= 2 theta, e'*x = 1, l <=x <=u

Here, the decision variable is "x". The constant vector "g" is estimated breeding values. 
The matrix "A" is Wright's numerator relationship matrix, and the theta is a threashold.
The vector "e" is the vector of all ones. The vectors "l" and "u" are the lower and upper bounds of "x", respectively.

This optimization problem was defined in 
  - T.H.E. Meuwissen, "Maximizing the response of selection with a predefined rate of inbreeding", Journal of Animal Science, Vol. 75, pp. 934-940, 1997.

An equally-type of optimization problem is of form:

max g'*x subject to x'*A*x <= 2 theta, e'*x = 1, x_i in {0,1/N}

Here, N is the given parameter, thus each genotype should contribute nothing (0) or the same amount (1/N).


For more details, please refer to the two papers below at "Papers."

For optimal selection problems, GENCONT by Meuwissen (http://www.genebankdata.cgn.wur.nl/gencont/gencont.html) is often used. The main advantage of this package is its computation speed. The compact SOCP formulation is also available through OPSEL (https://www.skogforsk.se/opsel/)



## Papers
The two methods were proposed in the two papers below.

1. the compact SOCP formulation for unequally deployment problem 
  - Makoto Yamashita, Tim J. Mullin, Sena Safarina, "An efficient second-order cone programming approach for optimal selection in tree breeding," Optimization Letters, Vol. 12 , No. 7, pp. 1683-1697, 2018. https://link.springer.com/article/10.1007/s11590-018-1229-y
```  
@article{Yamashita2018,
author="Yamashita, Makoto and Mullin, Tim J. and Safarina, Sena",
title="An efficient second-order cone programming approach for optimal selection in tree breeding",
journal="Optimization Letters",
year="2018",
month="Oct",
day="01",
volume="12",
number="7",
pages="1683--1697"
}
```
  
2. the steepest-ascent method for equally deployment problem
  - Sena Safarina, Satoko Moriguchi, Tim J. Mullin, and Makoto Yamashita, "Conic relaxation approaches for equal deployment problems," To appear in Discrete Applied Mathematics, 2019. https://www.sciencedirect.com/science/article/pii/S0166218X19304184

```
@article{SAFARINA2019,
title = "Conic relaxation approaches for equal deployment problems",
journal = "Discrete Applied Mathematics",
year = "2019",
issn = "0166-218X",
doi = "https://doi.org/10.1016/j.dam.2019.04.032",
url = "http://www.sciencedirect.com/science/article/pii/S0166218X19304184"
}
```

## Note

1. Ipopt, ECOS
2. Quaas algorithm

## Data

