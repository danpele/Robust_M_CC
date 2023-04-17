
pturnoverDN  = function(weights, rets, freq){
  results = Return.portfolio(R = rets,   
                             weights = weights, 
                             rebalance_on = freq, verbose = T)
  bop = results$BOP.Weight #beginning of period weights
  bop
  eop = results$EOP.Weight #end of period weights
  eop
  f = abs(bop - eop)
  out = sum(f)*(1/(length(ep) - 1)) #  
  return(out)
}

mad_ew = function(x){
  a = abs(x - 1/ncol(ret))
  return(a)
}
trans_cost = function(weights, rets, freq, c){
  results = Return.portfolio(R = rets,   
                             weights = weights, 
                             rebalance_on = freq, verbose = T)
  
  bop = results$BOP.Weight #beginning of period weights
  bop
  eop = results$EOP.Weight #end of period weights
  eop
  out = c*row_sums(abs(eop - bop))
  return(out)
} 

calculatePerformanceMeasures = function(start,end){
  collNumbers = vector()
  collNumbers_tc = vector()
  collres = xts()
  weightsNumbers = vector()
  for (stratloop in 1:length(strats)){
    
    Rb = as.xts(rescoll[,  which(strats %in% c("EqualWeight","EW"))], 
                order.by = as.Date(rownames(rescoll)))
    portfolioret_net = na.omit(rescoll[,stratloop])
    strat_weights = weightscoll[[stratloop]] 
    strat_weights[is.nan(strat_weights)] = 0.0
    portfolioret_net_xts = as.xts(as.matrix(na.omit(rescoll[,stratloop])), 
                                  order.by = as.Date(na.omit(rownames(rescoll))))
    portfolioEquity_net = 1 + cumsum(portfolioret_net)
    cumWealth = tail(portfolioEquity_net, 1)
    firstsignal = start
    rettt = portfolioret_net[firstsignal:end]
    rettt_xts = portfolioret_net_xts[firstsignal:end]
    ret_data = rets_log
    stock_rets = ret_data[firstsignal:end]
    tc = trans_cost(strat_weights[firstsignal:end,], stock_rets, freq, c = transi)
    portfolioret_net_tc = portfolioret_net[firstsignal:(end - 1)] - tc
    portfolioEquity_net_tc = 1 + cumsum(portfolioret_net_tc)
    cumWealth_tc = tail(portfolioEquity_net_tc, 1)
    T = (commonDate[end] - commonDate[firstsignal])/365
    Return.ann = (portfolioEquity_net[end]/portfolioEquity_net[firstsignal - 1])^(1/T) - 1
    Return.ann_tc = (tail(portfolioEquity_net_tc, 1)/portfolioEquity_net_tc[firstsignal - 1])^(1/T) - 1
    Vola.ann = sd(rettt)*sqrt(252);
    Vola.ann_tc = sd(portfolioret_net_tc)*sqrt(252);
    Sharpe.ann = Return.ann/Vola.ann
    Sharpe.ann_tc = Return.ann_tc/Vola.ann_tc
    target_turnover = vector();
    for (i in 2:dim(strat_weights)[1]) {
      target_turnover[i] = sum(abs(matrix(strat_weights[i, ]) - matrix(strat_weights[i - 1,])))/dim(strat_weights)[2] 
    }
    Turnover = mean(na.omit(target_turnover))
    value = portfolioEquity_net 
    ma = unlist(lapply(c(2:length(value)),function(x) max(value[1:x])))
    dddisc = value[-1]/ma - 1
    datums = commonDateR[firstsignal:end]
    num = as.numeric(tail(datums,1)-datums[1])
    PR = as.vector(PainRatio(rettt_xts))
    TurnoverDM  = pturnoverDN(strat_weights[firstsignal:end,], stock_rets, freq)
    #Return_annual = as.vector(Return.annualized(rettt_xts, geometric = F))
    AverageDrawdown = as.numeric(AverageDrawdown(rettt_xts))
    Sharpe =  as.numeric(SharpeRatio(rettt_xts))[1]
    StdDev.annualized = as.numeric(StdDev.annualized(rettt_xts))
    
    ###from Performace analytic
    
    Return_anual = as.vector(Return.annualized(rettt_xts))
    AverageDrawdown = as.numeric(AverageDrawdown(rettt_xts))
    AverageLength = as.numeric(AverageLength(rettt_xts))
    AverageRecovery = as.numeric(AverageRecovery(rettt_xts))
    BernardoLedoitRatio = as.numeric(BernardoLedoitRatio(rettt_xts))
    BurkeRatio = as.numeric(BurkeRatio(rettt_xts))
    CalmarRatio = as.numeric(CalmarRatio(rettt_xts))
    CDD = as.numeric(CDD(rettt_xts))
    CVaR = as.numeric(CVaR(rettt_xts))
    DownsideDeviation = as.numeric(DownsideDeviation(rettt_xts))
    DownsideFrequency = as.numeric(DownsideFrequency(rettt_xts))
    DownsidePotential = as.numeric(DownsidePotential(rettt_xts))
    DRatio = as.numeric(DRatio(rettt_xts))
    HurstIndex = as.numeric(HurstIndex(rettt_xts))
    #Kappa = as.numeric(Kappa(rettt_xts,MAR = 0.00,l = 2))
    KellyRatio = as.numeric(KellyRatio(rettt_xts))
    #Kurtosis = apply(rettt_xts,2,kurtosis) as.numeric(kurtosis(rettt_xts))
    MartinRatio = as.numeric(MartinRatio(rettt_xts))
    #Omega = as.numeric(Omega(rettt_xts))
    OmegaSharpeRatio = as.numeric(OmegaSharpeRatio(rettt_xts))
    #PainRatio = as.numeric(PainRatio(rettt_xts))
    ProspectRatio = as.numeric(ProspectRatio(rettt_xts,MAR = 0.00))
    SemiDeviation = as.numeric(SemiDeviation(rettt_xts))
    # Skewness = apply(rettt_xts,2,skewness) as.numeric(skewness(rettt_xts))
    SkewnessKurtosisRatio = as.numeric(SkewnessKurtosisRatio(rettt_xts))
    SortinoRatio = as.numeric(SortinoRatio(rettt_xts))
    StdDev.annualized = as.numeric(StdDev.annualized(rettt_xts))
    UpsidePotentialRatio = as.numeric(UpsidePotentialRatio(rettt_xts))
    UpsideRisk = as.numeric(UpsideRisk(rettt_xts))
    #pfInformationRatio=as.numeric(InformationRatio(rettt_xts,))
    VolatilitySkewness = as.numeric(VolatilitySkewness(rettt_xts))
    ###Check dimensions of collNumbers and latex table!!!
    collNumbers = cbind(collNumbers, as.vector(c(cumWealth, 100*Sharpe.ann,  
                                                 Turnover, TurnoverDM, AverageDrawdown , 
                                                 StdDev.annualized, CVaR, AverageLength , AverageRecovery , 
                                                 BernardoLedoitRatio , BurkeRatio ,
                                                 CDD , DownsideDeviation , DownsideFrequency , 
                                                 DownsidePotential , DRatio , HurstIndex , KellyRatio,
                                                 MartinRatio , 
                                                 OmegaSharpeRatio ,  ProspectRatio ,
                                                 SemiDeviation ,  SkewnessKurtosisRatio , 
                                                 SortinoRatio  , UpsidePotentialRatio , 
                                                 UpsideRisk , VolatilitySkewness)))#, pfInformationRatio)))#, 
    collNumbers_tc = cbind(collNumbers_tc, as.vector(c(cumWealth_tc,  
                                                       100*Sharpe.ann_tc, Turnover, TurnoverDM)))
    weightcoll = as.data.frame(strat_weights[first_signal:end])
    weightsNumbers = cbind(weightsNumbers, as.vector(c((mean(apply(weightcoll, 1, min))), mean(apply(weightcoll, 1, max)), 
                                                       mean(apply(weightcoll, 1, sd)), mean(apply(weightcoll, 1, mad_ew)), mean(diff(apply(weightcoll, 1, range))))))
    collNumbers = round(collNumbers,4)
    weightsNumbers = round(weightsNumbers,4)
    res = as.xts(portfolioret_net[first_signal:end],order.by=index(ret)[first_signal:end])
    
    res_tc = as.xts(portfolioret_net_tc[first_signal:end],order.by=index(ret)[first_signal:end])
    collres = cbind(collres,res)
    first(index(res))
    last(index(res))
    
  }
  return(list(collNumbers,collres,weightsNumbers, collNumbers_tc, res_tc))
}

#### HRP portfolio construction ####
getHRP = function(cov, corr, max = 1, min = 0, return_raw = NULL, robust_cov = FALSE) {
  #  Sigma,CR,max=maxweight, min=minweight, robust_cov = FALSE, return_raw = rettt
  #cov = Sigma#abs(cov)
  #corr = CR#cov2cor(cov)
  vols = StdDev(return_raw)
  covs = t(vols) %*% vols * corr
  # Construct a hierarchical portfolio
  if (robust_cov == T) {
    cov = cov_shrink(return_raw)
    corr = cov2cor(cov)
  }
  
  # Set the constraint matrix
  if (is.null(max)) max = rep(1,ncol(cov))
  if (is.null(min)) min = rep(0,ncol(cov))
  if (length(max)==1) max = rep(max,ncol(cov)) else if (length(max)<ncol(cov)) stop("Provide correct weights")
  if (length(min)==1) min = rep(min,ncol(cov)) else if (length(min)<ncol(cov)) stop("Provide correct weights")
  const = rbind(max, min)
  
  # check constraints
  if (sum(const[1,]) < 1 | sum(const[2,]) > 1) stop("Incompatible weights")
  diag(corr) = 1
  distmat = ((1 - corr) * 2)^0.5    # ((1 - corr) / 2)^0.5 #1-abs(corr) #
  #clustOrder = fastcluster::hclust(dist(distmat), method = dendprocedure)$order
  #clustOrder = fastcluster::hclust(dist(distmat, method ="euclidean"),dendprocedure)$order
  # hrp weights
  clustOrder = hclust(dist(corr), method = 'single')$order
  out = getRecBipart(covs, clustOrder)
  #out = getRecBipart(cov, clustOrder, const, type)
  return(out)
}

##getRecBipart function from https://www.r-bloggers.com/2017/05/the-marcos-lopez-de-prado-hierarchical-risk-parity-algorithm/

getRecBipart <- function(covMat, sortIx) {
  # keeping track of w in the global environment
  assign("w", value = rep(1, ncol(covMat)), envir = .GlobalEnv)
  recurFun(covMat, sortIx)
  return(w)
}

recurFun <- function(covMat, sortIx) {
  subIdx <- 1:trunc(length(sortIx)/2)
  cItems0 <- sortIx[subIdx]
  cItems1 <- sortIx[-subIdx]
  cVar0 <- getClusterVar(covMat, cItems0)
  cVar1 <- getClusterVar(covMat, cItems1)
  alpha <- 1 - cVar0/(cVar0 + cVar1)
  
  # scoping mechanics using w as a free parameter
  w[cItems0] <<- w[cItems0] * alpha
  w[cItems1] <<- w[cItems1] * (1-alpha)
  
  if(length(cItems0) > 1) {
    recurFun(covMat, cItems0)
  }
  if(length(cItems1) > 1) {
    recurFun(covMat, cItems1)
  }
}
getIVP <- function(covMat) {
  invDiag <- 1/diag(as.matrix(covMat))
  weights <- invDiag/sum(invDiag)
  return(weights)
}

getClusterVar <- function(covMat, cItems) {
  covMatSlice <- covMat[cItems, cItems]
  weights <- getIVP(covMatSlice)
  cVar <- t(weights) %*% as.matrix(covMatSlice) %*% weights
  return(cVar)
}


