![OpSel](https://github.com/makoto-yamashita/OpSel.jl/blob/logo-registration/logo/opsel-logo.png)

![travis](https://travis-ci.org/makoto-yamashita/OpSel.jl.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/makoto-yamashita/OpSel.jl/badge.svg?branch=master)](https://coveralls.io/github/makoto-yamashita/OpSel.jl?branch=master)

# OpSel.jl
Efficient optimal selection for tree breeding

Optimal selection problems are to find optimum of selection of genotypes that maximizes genetic gains under a constraint on genetic diversity which involves Wright's numerator relationship matrix.

Using a conic programming approach, this package provides efficient numerical methods for optimal selection problems arising from tree breeding. This package now implements two methods
1. the compact SOCP formulation for unequally deployment problem
2. the steepest-ascent method for equally deployment problem

## Installation 

```import Pkg; Pkg.add("OpSel")```

or in the package mode by `]`,

```pkg> https://github.com/makoto-yamashita/OpSel.jl```

## Usage

The dataset of the package includes the sizes of Z = 200, 1050, 2045, 5050, 5255, 10100, 15100, and 15222; and N = 50 and 100.

1. To execute the compact SOCP formulation with the dataset of the package (for <a href="https://www.codecogs.com/eqnedit.php?latex=Z=2045" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Z=2045" title="Z=2045" /></a>)
```
Using OpSel
OpSel.testUnequal(2045);
```
Afther the execution of Ipopt, the summary will be reported as follows:
```
=== Result Summary ===
JuMP status = LOCALLY_SOLVED
Z = 2045, N_s = 14, N = 2800
gx = 439.116285, xAx = 0.071429, 2theta = 0.071429
time(s): build = 0.049, solver = 1.392, total = 1.441
```
The obtained objective value is <a href="https://www.codecogs.com/eqnedit.php?latex=g^T&space;x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?g^T&space;x" title="g^T x" /></a>=439.116825, and we can see that the quadratic constraint <a href="https://www.codecogs.com/eqnedit.php?latex=x^T&space;A&space;x&space;\le&space;2&space;\theta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x^T&space;A&space;x&space;\le&space;2&space;\theta" title="x^T A x \le 2 \theta" /></a> is satisfied.
The computation time to build the SOCP is 0.049 seconds, the time by Ipopt is 1.392 seconds, and the entire time is 1.441 seconds.


2. To execute the steepest-ascent method with the dataset of the package (for <a href="https://www.codecogs.com/eqnedit.php?latex=Z=200,&space;N=50" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Z=200,&space;N=50" title="Z=200, N=50" /></a>)
```
Using OpSel
OpSel.testEqual(200,50);
```
Afther the execution, the summary will be reported as follows:
```
SOCP = 26.155540, STEEP = 25.090400, gap = 4.072330%
time(s): build = 0.016, solver = 0.088, steep = 0.820, total = 0.925
```
The objective value of the SOCP relaxation is 26.155540, and the objective value after the steepst-ascent method is 25.090440.
The gap is computed from these two objective values.
The exact optimal value must be between these two objective values, thus if the gap is 0, we can obtain the exact optimal value.
In the computation time, steep corresponds to the steepest-ascent part.

3. If an input CSV file is available, the methods can be executed as follows:

```
sp_csv = OpSel.loadFile(filename)
(x_result, info_result) = OpSel.compactSOCP(sp_csv)
```
or
```
sp_csv = OpSel.loadFile(filename)
(x_result, info_result) = OpSel.steepestAscent(sp_csv, N=50)
```

The columns of the input CSV file should be:
```i(id), p(parent1), p(parent2), g(EBV), u(upperbound), h(inbreeding)```
For example, the following line in a CSV
```1040, 782, 751, 3.1313800000000001, 50, 0.0410156250000000```
indicates that the parents of 1040 are 782 and 751. The EBV of 1040 is 3.13138, and the upperbound is 50 (this will be divided by 2800) and the inbreeding value is 0.041015625.

## Basic Formulation

An unequally-type of optimization problem is of form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\max:&space;g^T&space;x&space;\&space;\&space;\text{subject&space;to:}&space;x^TAx&space;\le&space;2\theta,&space;e^T&space;x&space;=&space;1,&space;l\le&space;x&space;\le&space;u" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\max:&space;g^T&space;x&space;\&space;\&space;\text{subject&space;to:}&space;&space;x^TAx&space;\le&space;2\theta,&space;e^T&space;x&space;=&space;1,&space;l\le&space;x&space;\le&space;u" title="\max: g^T x \ \ \text{subject to:} \  x^TAx \le 2\theta, e^T x = 1, l\le x \le u" /></a>

Here, the decision variable is <a href="https://www.codecogs.com/eqnedit.php?latex=x&space;\in&space;\mathbb{R}^Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x&space;\in&space;\mathbb{R}^Z" title="x \in \mathbb{R}^Z" /></a>. The constant vector <a href="https://www.codecogs.com/eqnedit.php?latex=g&space;\in&space;\mathbb{R}^Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?g&space;\in&space;\mathbb{R}^Z" title="g \in \mathbb{R}^Z" /></a> is estimated breeding values. 
The matrix "A" is Wright's numerator relationship matrix, and the theta is a threashold.
The vector <a href="https://www.codecogs.com/eqnedit.php?latex=e&space;\in&space;\mathbb{R}^Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?e&space;\in&space;\mathbb{R}^Z" title="e \in \mathbb{R}^Z" /></a> is the vector of all ones. The vectors <a href="https://www.codecogs.com/eqnedit.php?latex=l&space;\in&space;\mathbb{R}^Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?l&space;\in&space;\mathbb{R}^Z" title="l \in \mathbb{R}^Z" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=u&space;\in&space;\mathbb{R}^Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u&space;\in&space;\mathbb{R}^Z" title="u \in \mathbb{R}^Z" /></a> are the lower and upper bounds of <a href="https://www.codecogs.com/eqnedit.php?latex=x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x" title="x" /></a>, respectively.

This optimization problem was defined in 
  - T.H.E. Meuwissen, "Maximizing the response of selection with a predefined rate of inbreeding", Journal of Animal Science, Vol. 75, pp. 934-940, 1997.

An equally-type of optimization problem is of form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\max:&space;g^T&space;x&space;\&space;\&space;\text{subject&space;to:}&space;\&space;x^TAx&space;\le&space;2\theta,&space;e^T&space;x&space;=&space;1,&space;x_1,&space;\ldots,&space;x_n&space;\in&space;\left\{0,\frac{1}{N}\right\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\max:&space;g^T&space;x&space;\&space;\&space;\text{subject&space;to:}&space;\&space;x^TAx&space;\le&space;2\theta,&space;e^T&space;x&space;=&space;1,&space;x_1,&space;\ldots,&space;x_n&space;\in&space;\left\{0,\frac{1}{N}\right\}" title="\max: g^T x \ \ \text{subject to:} \ x^TAx \le 2\theta, e^T x = 1, x_1, \ldots, x_n \in \left\{0,\frac{1}{N}\right\}" /></a>

Here, <a href="https://www.codecogs.com/eqnedit.php?latex=N" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N" title="N" /></a> is the given parameter, thus each genotype should contribute nothing <a href="https://www.codecogs.com/eqnedit.php?latex=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?0" title="0" /></a> or the same amount <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{1}{N}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{1}{N}" title="\frac{1}{N}" /></a>.


For optimal selection problems, GENCONT by Meuwissen (http://www.genebankdata.cgn.wur.nl/gencont/gencont.html) is often used. The main advantage of this package is its computation speed. The compact SOCP formulation is also available through OPSEL (https://www.skogforsk.se/opsel/).
For more details, please refer to the two papers below at "Papers."

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

## Data

The datasets with sizes 2045 and 15222 were avaiable at the Dryad Digital Repository (http://dx.doi.org/10.5061/dryad.9pn5m). The other
data were produced with POPSIM (https://www.skogforsk.se/popsim/). In each dataset, the candidates are sorted so that all the ids are sequential numbers.

## Notes

1. In the two papers above, ECOS (https://github.com/JuliaOpt/ECOS.jl) was used as the SOCP solver. However, recent versions of ECOS was not numerically stable for optimal selection problems. Instead, Ipopt (https://github.com/JuliaOpt/Ipopt.jl) is used in this package.
2. The inbreeding value (the last column of the input CSV file) was computed by Quaas's algorithm
  - R. L. Quaas, "Computing the diagonal elements and inverse of a large numerator relationship matrix," Biometrics, Vol. 32, pp. 949–953, 1976.

