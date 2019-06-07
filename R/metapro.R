#' @title Beta probability
#' @param p p-value
#' @param i rank
#' @param n The number of inputs
#' @importFrom "stats" "pbeta" "qbeta"
F_i = function(p, i, n)
{
  a = i
  b = n-i+1
  res = pbeta(q = p, shape1 = a, shape2 = b, lower.tail = T)
  return(res)
}


#' @title MMPO
#' @description Minimum Marginal P-value in joint order distribution
#' @param p A vector of p-values
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @import rJava rSymPy
#' @return p : Combined p-value
#' @return optimal_rank : Optimal rank where minimum marginal p-value exists.
#' @return eff.p.idx : Index of effective p-values
#' @return MMP : Minimum marginal p-value
#' @return overall.eff.direction : The direction of combined effects.
#' @examples MMPO(p=c(0.01, 0.02, 0.2, 0.5, 0.8), is.onetail=FALSE, eff.sign=c(1,1,1,1,1))
#' @export

MMPO = function(p, is.onetail = TRUE, eff.sign=NULL)
{
  direc = eff.sign
  if(is.null(p)){stop("Input p-values are required.")}
  if(!is.onetail & is.null(eff.sign)){stop("Input the direction of effects.")}
  idx_na = which(is.na(p))
  if(length(idx_na)>0){p = p[-idx_na]; eff.sign = eff.sign[-idx_na]}
  ordmeta = function(p2)
  {
    ord = order(p2, decreasing = F)
    pord = sort(p2, decreasing = F)

    # get alpha = MIN(F_(i)(x)) {i={1..n}}
    N = length(p2)
    alpha = 1.01 # an arbitrary number larger than 1
    for(i in 1:N)
    {
      alpha_temp = F_i(pord[i], i, N)
      if(alpha_temp < alpha){idx_minimum = i; alpha = alpha_temp}
    }
    # symbolic integral
    for(i in 1:N)
    {
      x = Var("x")
      y = Var("y")
      if(i==1)
      {
        templete = paste(i,"*integrate(1, (x, lob, y))")
        lob = qbeta(p = alpha,shape1 = i, shape2 = N-i+1, lower.tail = T)
        templete = gsub("lob", lob, templete)
      }else if(i>1 & i<N){
        integ = gsub(pattern = "y", replacement = "x", x = integ)
        templete = paste(i, "*integrate(",integ,", (x, lob, y))")
        lob = qbeta(p = alpha,shape1 = i, shape2 = N-i+1, lower.tail = T)
        templete = gsub("lob", lob, templete)
      }else if(i==N)
      {
        integ = gsub(pattern = "y", replacement = "x", x=integ)
        templete = paste(i, "*integrate(",integ,", (x, lob, 1))")
        lob = qbeta(p = alpha,shape1 = i, shape2 = N-i+1, lower.tail = T)
        templete = gsub("lob", lob, templete)
      }
      #print(templete)
      integ = sympy(templete)
    }
    res = 1-as.numeric(integ)
    return(list(p=res, optimal_rank = idx_minimum, eff.p.idx = ord[1:idx_minimum], MMP = alpha))
  }
  if(is.onetail)
  {
    RES = ordmeta(p2 = p)
    return(RES)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign > 0)
    idx_neg = which(eff.sign < 0)
    p1[idx_pos] = p/2
    p1[idx_neg] = 1-p/2
    p2[idx_pos] = 1-p/2
    p2[idx_neg] = p/2

    RES1 = ordmeta(p2 = p1)
    RES2 = ordmeta(p2 = p2)
    if(RES1$p<=RES2$p){
      RES = RES1; RES$overall.eff.direction = "+"
    }else{
        RES = RES2; RES$overall.eff.direction = "-"
    }
    RES$p = RES$p * 2
    return(RES)
  }
}


#' @title wFisher
#' @description sample size-weighted Fisher's method
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weight or sample size for each experiment
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effects, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @return p : Combined p-value
#' @return overall.eff.direction : The direction of combined effects.
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
#' @return p : Combined p-value
#' @return overall.eff.direction : The direction of combined effects.
#' @return sumz : Sum of transformed z-score
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
#' @return p : Combined p-value
#' @return overall.eff.direction : The direction of combined effects.
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
