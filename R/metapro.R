ordmeta_mono <- function(p) {
  .Call('_metapro_ordmeta_mono', PACKAGE = 'metapro', p)
}

ordmeta_signed <- function(p, effect_size, sort_by_decreasing = TRUE) {
  .Call('_metapro_ordmeta_signed', PACKAGE = 'metapro', p, effect_size, sort_by_decreasing)
}

#' @title MMPO
#' @description Minimum Marginal P-value in joint order distribution
#' @param p A vector of p-values
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @import Rcpp
#' @return A list of Combined p-value based on MMPO method, the indices of effecitve p-values and NA values.
#' @examples MMPO(p=c(0.01, 0.02, 0.2, 0.5, 0.8), is.onetail=FALSE, eff.sign=c(1,1,1,1,1))
#' @export
#' @useDynLib metapro

MMPO = function(p, is.onetail=TRUE, eff.sign)
{
  if(is.onetail){res = ordmeta_mono(p)}
  if(!is.onetail){
    res1 = ordmeta_signed(p, eff.sign, sort_by_decreasing = TRUE)
    res2 = ordmeta_signed(p, eff.sign, sort_by_decreasing = FALSE)
    if(res1$p <= res2$p){
      res = res1
      res$p = res$p * 2
      overall.eff.direction = "+"
    } else {
      res = res2
      res$p = res$p * 2
      overall.eff.direction = "-"
    }
    res$overall.eff.direction = overall.eff.direction
  }
  idx = which(names(res) == "index.NA")
  res = res[-idx]
  return (res)
}

#' @title wFisher
#' @description sample size-weighted Fisher's method
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weight or sample size for each experiment
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effects, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @return A Combined p-value
#' @examples wFisher(p=c(0.01,0.2,0.8), weight = c(50,60,100),is.onetail=FALSE, eff.sign=c(1,1,1))
#' @importFrom "stats" "pgamma" "qgamma"
#' @export

wFisher = function(p, weight = NULL, is.onetail = TRUE, eff.sign)
{
  if(is.null(weight)){weight = rep(1, length(p))}
  idx.na = which(is.na(p))
  if(length(idx.na)>0){
    p = p[-idx.na];
    weight = weight[-idx.na];
    if(!is.onetail)
    {
      eff.sign = eff.sign[-idx.na]
    }
  }
  NP = length(p)
  NS = length(weight)
  if(NP!=NS){stop("The length of p and weight vector must be identical.")}
  N = NS
  Ntotal = sum(weight)
  ratio = weight/Ntotal
  Ns = N*ratio
  G = c()

  if(is.onetail)
  {
    for(i in 1:length(p))
    {
      G = append(G, qgamma(p = p[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign > 0)
    idx_neg = which(eff.sign < 0)
    # positive direction
    G = c()
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1-p[idx_neg]/2
    for(i in 1:length(p1))
    {
      G = append(G, qgamma(p = p1[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP1 = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
    # negative direction
    G = c()
    p2[idx_pos] = 1-p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2
    for(i in 1:length(p2))
    {
      G = append(G, qgamma(p = p2[i], shape = Ns[i], scale=2, lower.tail=F))
    }
    Gsum = sum(G)
    resultP2 = pgamma(q=Gsum, shape=N, scale=2, lower.tail=F)
    resultP = 2* min(resultP1, resultP2)
    overall.eff.direction = if(resultP1<=resultP2){"+"}else{"-"}
  }
  RES = if(is.onetail){list(p=min(1,resultP))}else{list(p=min(1,resultP), overall.eff.direction=overall.eff.direction)}
  return(RES)
}


#' @title wZ
#' @description P-value combination based on weighted Z-method
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weights (e.g., sample sizes)
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @return Combined p-value
#' @examples wZ(p=c(0.01,0.2,0.8), weight = c(20,10,40), is.onetail=FALSE, eff.sign=c(1,-1,1))
#' @import "metap"
#' @export
wZ = function(p, weight = NULL, is.onetail = TRUE, eff.sign)
{
  if(length(p) == 0){stop("No input p-values exist.")}
  if(is.null(weight)){weight = rep(1, length(p))}
  idx.na = which(is.na(p))
  if(length(idx.na)>0){
    p = p[-idx.na];
    weight = weight[-idx.na]
    if(!is.onetail)
    {
      eff.sign = eff.sign[-idx.na]
    }
  }
  if(is.onetail){res = sumz(p = p, weights = weight)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign > 0)
    idx_neg = which(eff.sign < 0)
    # positive direction
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1 - p[idx_neg]/2
    res1 = sumz(p = p1, weights = weight)
    # negative direction
    p2[idx_pos] = 1 - p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2
    res2 = sumz(p = p2, weights = weight)

    if(res1$p <= res2$p){res = res1;overall.eff.direction = "+"}else{res = res2; overall.eff.direction = "-"}
  }
  if(is.onetail)
  {
    RES = list(p = res$p[,1], sumz = res$z[,1])
  }else{
    RES = list(p = res$p[,1], overall.eff.direction = overall.eff.direction, sumz = res$z[,1])
  }

  return(RES)
}

#' @title Lancaster
#' @description P-value combination based on Lancaster's procedure
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weights (e.g., samples sizes)
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes (1 or -1). It works when is.onetail = FALSE
#' @return Combined p-value
#' @examples lancaster(p=c(0.01,0.2,0.8), weight=c(20,50,10), is.onetail=FALSE, eff.sign=c(1,1,1))
#' @importFrom "stats" "pchisq" "qchisq"
#' @export

lancaster = function(p, weight, is.onetail=TRUE, eff.sign)
{
  if(is.null(p)){stop("Input p-values are required.")}
  if(is.null(weight)){stop("Weights are required.")}
  idx.na = which(is.na(p))
  if(length(idx.na)>0){
    p = p[-idx.na];
    weight = weight[-idx.na];
    if(!is.onetail)
    {
      eff.sign = eff.sign[-idx.na]
    }
  }

  NP = length(p)
  NS = length(weight)
  if(NP!=NS){stop("The length of pvalue and samplesize vector must be identical")}
  N = NP
  Ntotal = sum(weight)
  SS = weight
  #Ns = N*ratio
  if(is.onetail)
  {
    X = c()
    for(i in 1:N)
    {
      X = append(X, qchisq(p = p[i],df = SS[i], lower.tail=F))
    }
    Xsum = sum(X)
    resultP = pchisq(q=Xsum, df = Ntotal, lower.tail=F)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign > 0)
    idx_neg = which(eff.sign < 0)
    # positive direction
    X = c()
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1-p[idx_neg]/2
    for(i in 1:N)
    {
      X = append(X, qchisq(p = p1[i],df = SS[i], lower.tail=F))
    }
    Xsum = sum(X)
    resultP1 = pchisq(q=Xsum, df = Ntotal, lower.tail=F)
    # negative direction
    X = c()
    p2[idx_pos] = 1-p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2
    for(i in 1:N)
    {
      X = append(X, qchisq(p = p2[i],df = SS[i], lower.tail=F))
    }
    Xsum = sum(X)
    resultP2 = pchisq(q=Xsum, df = Ntotal, lower.tail=F)
    resultP = 2* min(resultP1, resultP2)
    overall.eff.direction = if(resultP1<=resultP2){"+"}else{"-"}
  }
  RES = if(is.onetail){list(p=min(resultP,1))}else{list(p=min(resultP,1), overall.eff.direction=overall.eff.direction)}
  return(RES)
}
