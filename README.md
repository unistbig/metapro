# metapro


## Introduction
The meta-analysis is conducted to increase the statistical power by combining evidences (e.g., effect sizes and p-values) obtained from multiple experiments. The p-value combination has been widely used for meta-analysis when effect sizes are not available. 
<i>metapro </i> is a CRAN R package that provides functions for p-value combinations. There are four functions including

1. ordmeta : Minimum marginal p-value in joint order distribution
2. wFisher : weighted Fisher's method
3. wZ : weighted Z-method
4. lancaster : Lancaster's procedure

The original Fisher's method and (weighted) Z-method have been commonly used for p-value combination. In particular, the Z-method is effective when a feature is associated with most of the cohorts. However, this attributes is disadvantageous detecting partially associated patterns (e.g., detection of African-specific features in trans-ethnic analysis).Therefore, in this study, we designed the ordmeta and wFisher method that work effectively in those situations. The detailed formulas are described in our paper. 
  
## Installation 
<b> NOTE! </b><br>
If JAVA is not installed on your PC, please <b>install JAVA first</b> because this package requires <i>'rJava'</i> package. <br> 
Download: <a href="https://www.java.com/en/">https://www.java.com/en/</a><br>

After the JAVA is installed, open R and install the <i>metapro</i> package by typing 

```
install.packages('devtools') # install 'devtools'
library(devtools)
install_github('unistbig/metapro', INSTALL_opts=c("--no-multiarch"))
library(metapro) # Load metapro package
sympyStart()
```
or 
```
install.packages('metapro', INSTALL_opts=c("--no-multiarch")) 
library(metapro)
sympyStart()
```
<b> Tip: Dealing with Errors </b><br>
1. Encoding Error
If you are using Windows, the following errors may occur:
```
Error in nchar(object, type = "chars") : 
  invalid multibyte string, element 1
```
In this case, type following in R console, and install the package again.
```
Sys.setlocale('LC_ALL','C')
```
2. Error from JAVA setting <br>
Most error will occur with rJAVA. 
If the error message contains 'JAVA', please try followings. <br>
  * Check whether both JAVA and R are 64-bit.
  * Set proper environment variable for JAVA_HOME, CLASSPATH, and RPATH. For example,
  ```
  Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jdk1.8.0_211")
  Sys.setenv(CLASSPATH="C:\\Program Files\\Java\\jdk1.8.0_211\\jre\\lib\\ext")
  Sys.setenv(RPATH="C:\\Program Files\\R\\R-3.6.0\\bin\\x64")
  ```
  Next, check whether rJAVA is installed correctly, and try installing metapro again.
  ```
 install.packages('rJava')
 library(rJava)
  ```

## Usage

__1. ordmeta__
-----
ordmeta combines p-value based on the minimum marginal p-value in the joint order distribution. In the result, it reports the combined p-value, optimal rank, effective p-value indices, marginal p-value at the optimal rank (MMP), and the direction of combined effect size (for two-tailed test). <br>
#### Example
`ordmeta(p=c(0.01, 0.02, 0.8, 0.25), is.onetail=FALSE, eff.sign = c(1,1,1,-1))`
<br><span font-size:1em>  will return </span>

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
<br><span font-size:1em>  will return </span>

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
<br><span font-size:1em>  will return </span>

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
<br><span font-size:1em>  will return </span>

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

## Contact
Sora Yoon: <yoonsora1@unist.ac.kr> <br>
Department of Biological Sciences, UNIST
