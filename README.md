# metapro


## Introduction
In meta-analysis, the results obtained from multiple experiments are combined to increase the statistical power. The p-value combination has been widely used for meta-analysis when effect sizes are not available. 
<i>metapro </i> is a CRAN R package that provides functions for p-value combinations. There are four functions including

1. MMPO : Minimum marginal p-value in joint order distribution
2. wFisher : weighted Fisher's method
3. wZ : weighted Z-method
4. lancaster : Lancaster's procedure

The original Fisher's method and (weighted) Z-method have been widely used for p-value combination. In particular, the Z-method is effective when the effect sizes of most experiments are not zero. However, this attributes is disadvantageous when detecting partially associated patterns (e.g., detection of African-specific features in trans-ethnic analysis).Therefore, in this study, we designed the MMPO and wFisher method that work effectively in that situation. The detailed formulas are described in our paper 
  
## Installation 
```
install.packages('devtools') # install 'devtools'
library(devtools)
install_github('unistbig/metapro')
library(metapro)
```
or 
```
install.packages('metapro') 
library(metapro)
```

## Usage

__1. MMPO__
-----
MMPO combines p-value based on the minimum marginal p-value in the joint order distribution. In result, it reports the optimal rank. 
#### Example
`MMPO(p=c(0.01, 0.02, 0.8, 0.25), is.onetail=FALSE, eff.sign = c(1,1,1,-1))`
will return

```
$p
[1] 0.004427233

$optimal_rank
[1] 2

$eff.p.idx
[1] 1 2

$MMP
[1] 0.00059203

$overall.eff.direction
[1] "+"
```
#### Input Arguments and Values
|<b>Arguments</b>|   |
|---|---|
| <b> p </b>    |  A numeric vector of p-values |
|  <b> is.onetail	</b> |  Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE. |
|  <b> eff.sign </b>  | A vector of signs of effects. It works when is.onetail = FALSE  |

| <b>Value</b>  	|   |
|---	|---	|
|  <b> p </b> 	|   Combined p-value	|
|  <b> optimal_rank </b>	| The optimal rank where minimum marginal p-value exists.   	|
|   <b> eff.p.idx </b>	|   Indexes of effective p-values 	|
|  <b>MMP</b> 	|   Minimum marginal p-value  	|
|   <b>overall.eff.direction</b>	|  The direction of combined effects. 	|

__2. wFisher__
------
wFisher is designed to assign weights to each experiment based on the sample size. <br>

#### Example
`wFisher(p=c(0.01, 0.02, 0.8, 0.25), weight = c(200, 500, 100, 80), is.onetail=FALSE, eff.sign = c(1,1,1,-1))`
will return

```
$p
[1] 0.003060523

$overall.eff.direction
[1] "+"
```
#### Input Arguments and Values
| Arguments  	|   	|
|---	|---	|
|   <b> p </b>	|   	A numeric vector of p-values	|
|   <b> is.onetail	</b> 	|    Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.	|
|   <b> weight </b>	|  A numeric vector of weight or sample size for each experiment. <b>Note! </b> If no weight option is given, the original Fisher method is performed.  	|
|  <b> eff.sign </b> 	| A vector of signs of effects. It works when is.onetail = FALSE  |


| <b>Value</b>  	|   |
|---	|---	|
|  <b> p </b> 	|   Combined p-value	|
|   <b>overall.eff.direction</b>	|  The direction of combined effects. 	|

__3. wZ__
------
Weighted Z-method. This function has been modified from sumz function in <i>metap</i> package. <br>

#### Example
`wZ(p=c(0.01, 0.02, 0.8, 0.25), weight = c(200, 500, 100, 80), is.onetail=FALSE, eff.sign = c(1,1,1,-1))`
will return

```$p
[1] 0.001798156

$overall.eff.direction
[1] "+"

$sumz
[1] 2.911558
```

#### Input Arguments and Values
| Arguments  	|   	|
|---	|---	|
|   <b> p </b>	|   	A numeric vector of p-values	|
|   <b> is.onetail	</b> 	|    Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.	|
|   <b> weight </b>	|  A numeric vector of weight or sample size for each experiment. <b>Note! </b> If no weight option is given, the Stouffer's method is performed.  	|
|  <b> eff.sign </b> 	| A vector of signs of effects. It works when is.onetail = FALSE  |

| <b>Value</b>  	|   |
|---	|---	|
|  <b> p </b> 	|   Combined p-value	|
|   <b>overall.eff.direction</b>	|  The direction of combined effects. 	|
|   <b>sumz</b>	|  Transformed sum of z-values 	|


__3. lancaster__
------
Lancaster's procedure - the generalized version Fisher's method <br>

#### Example
`lancaster(p=c(0.01, 0.02, 0.8, 0.25), weight = c(200, 500, 100, 80), is.onetail=FALSE, eff.sign = c(1,1,1,-1))`
will return

```$p
$p
[1] 0.005694935

$overall.eff.direction
[1] "+"
```
#### Input Arguments and Values
| <b>Arguments</b> 	|   	|
|---	|---	|
|   <b> p </b>	|   	A numeric vector of p-values	|
|   <b> is.onetail	</b> 	|    Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.	|
|   <b> weight </b>	|  A numeric vector of weight or sample size for each experiment. <b>Required! </b>  	|
|  <b> eff.sign </b> 	| A vector of signs of effects. It works when is.onetail = FALSE  |
  
| <b>Value</b>  	|   |
|---	|---	|
|  <b> p </b> 	|   Combined p-value	|
|   <b>overall.eff.direction</b>	|  The direction of combined effects. 	|
