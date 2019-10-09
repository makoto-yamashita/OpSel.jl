![OpSel](https://github.com/makoto-yamashita/OpSel.jl/blob/logo-registration/logo/opsel-logo.png)

# OpSel.jl
Efficient optimal selection for tree breeding

Optimal selection problems are to find optimum of selection of genotypes that maximizes genetic gains under a constraint on genetic diversity which involves Wright's numerator relation ship.

Using a conic programming approach, this package provides efficient numerical methods for optimal selection problems arising from tree breeding.
This package now implements two methods
1. the compact SOCP formulation for unequally deployment problem <sup>[1](#myfootnote1)</sup>
2. the steepest-ascent method for equally deployment problem <sup>[2](#myfootnote2)</sup>

<a name="myfootnote1">1</a>: Makoto Yamashita, Tim J. Mullin, Sena Safarina, "An efficient second-order cone programming approach for optimal selection in tree breeding," Optimization Letters , Vol. 12 , No. 7. pp 1683-1697, 2018.

<a name="myfootnote2">1</a>: Sena Safarina, Satoko Moriguchi, Tim J. Mullin, and Makoto Yamashita, "Conic relaxation approaches for equal deployment problems," To appear in Discrete Applied Mathematics, 2019.



## Installation 

## Usage

## Basic Formulation

An unequally-type of optimization problem is of form:

max g'*x subject to x'*A*x <= 2 theta, e'*x = 1, l <=x <=u

Here, the decision variable is "x". The constant vector "g" is estimated breeding values. 
The matrix "A" is Wright's numerator relationship matrix, and the theta is a threashold.
The vector "e" is the vector of all ones. The vectors "l" and "u" are the lower and upper bounds of "x", respectively.

For more details, 

For optimal selection problems, GENCONT by Meuwissen (http://www.genebankdata.cgn.wur.nl/gencont/gencont.html) is often used. The numerical methods implemented this package can solve the problems much faster than GENCONT. 



## Citation
BibTeX file is avialbe from
https://citation-needed.springer.com/v2/references/10.1007/s11590-018-1229-y?format=bibtex&flavour=citation

## Note

1. Ipopt, ECOS
2. Quaas algorithm

## Data

