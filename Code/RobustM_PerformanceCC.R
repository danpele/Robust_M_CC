###### SETTINGS #####
# Clean the environment 
graphics.off()
rm(list = ls(all = TRUE))
setwd("D:\\PROIECTE\\Robust Markowitz CC\\Robust_M_CC_19_01_2023")
# Load Functions and other Files
source('./RobustM_Packages.R')
source('./RobustM_Functions.R')
source('factor_an_static.R')
#install.packages("fitHeavyTail")

# Working directory
print('Please use the code from the project in order to maintain a correct functionality of the code')

#Choose dataset to analyse
df = "Assets"
load("Assets.RData")
load("Assets_classes.RData")

#Choose the starting point of analysis  
start_point = "2016-01-04" 
myData = switch(df,
                "Assets" = Assets,
                "SP100" = SP100)
insample = F
random   = T 

#### Python environment ####

path_to_python = "C:\\Users\\dell\\anaconda3\\envs\\r-reticulate\\python.exe"

#define the path to Python
Sys.which("python")
Sys.setenv(RETICULATE_PYTHON = path_to_python)
use_python(path_to_python,  required = TRUE)
reticulate::py_config()
use_condaenv("r-reticulate")
# 
# 
# #Uncomment if set up the Python-R connection for the first time
# #create a new environment if you run the code for the first time
# conda_create("r-reticulate",python_version="3.10")
# #install Packages to the environment
# conda_install("r-reticulate", "scipy")
# conda_install("r-reticulate", "pandas")
# conda_install("r-reticulate", "numpy")
# conda_install("r-reticulate", "matplotlib")
# conda_install("r-reticulate", "scikit-learn")
# conda_install("r-reticulate", "statsmodels")
# conda_install("r-reticulate", "seaborn")
# conda_install("r-reticulate", "cvxopt")#from cvxopt import solvers, matrix, spmatrix
# # import ("pprint")

np     = import("numpy")
pd     = import("pandas")
mpl    = import("matplotlib")
plt    = import("matplotlib.pyplot")
cvxopt = import("cvxopt")
mrkw   = import("markowitz")

##Create Vectors of prices, returns 

  rets        = myData$Rets[, ]
  rets_log    = myData$LogRets[, ]
  commonDateR = index(myData$Data)
  prices      = as.matrix(myData$Data)
  commonDate  = as.numeric(commonDateR)

####Parameters for Portfolio Strategies####

ret                          = rets_log
freq                         = "months"
transi                       = 0.005 # transactional costs ratio
end                          = dim(ret)[1]
lambda                       = 2
hel                          = ret;
hel[]                        = 0
minweight                    = 0
maxweight                    = 1 
dendprocedure                = "complete"
#gammas = 1.0*np$exp(-np$linspace(0, 4, as.integer(2000)))

####Portfolios constraction####

  w_length = 1

if (insample) {
  {first_signal = endpoints(commonDateR,on = "years")[1]}
} else {first_signal = endpoints(commonDateR,on = "years")[w_length + 1]}
ep = endpoints(commonDateR,on = freq) # rebalancing dates
ep = ep[(ep >= first_signal) & (ep < (end-1))]
ret_data   = list()
ret_random = list()

    for (t in commonDateR[ep]) {#building data for estimation of parameters
      t = as.POSIXlt(as.Date(t)) 
      start_date = t
      start_date$year = t$year - w_length # date of begining of the estimation window
      ret_data[[as.character(as.Date(t))]] = ret[which(index(ret) > as.Date(start_date) & index(ret) <= as.Date(t)),]
    }


  strats = c("EW",   "GMV_long", "GMV",
             "GMV_lin", "GMV_nlin", "invvol", "ERC", "MD","GMV_robust","HRP")

#Settings for parallel computing
rescoll     = list()
weightscoll = list()
UseCores  = detectCores() - 1
cl = parallel::makeCluster(UseCores, setup_strategy = "sequential")
registerDoParallel(cl)
ptm = proc.time()


results = foreach::foreach(stratloop = strats, .packages=libraries, .combine = cbind) %dopar% {
 #for (stratloop in strats) {
  weightcoll = matrix()
  retttscoll = list()
  ind =  which(stratloop %in% strats)
  end = dim(ret)[1]
  weight = hel;
  period = NA
  period <- c(first_signal:end)#(end - 1)
  if (insample) period = c(1:(end))
  for (t in period) {
    if (!is.na(match(t,ep))) { 
      retT = ret_data[[as.character(as.Date(commonDateR[t]))]]

      n = dim(retT)[2]
      rettt = retT
      w1 = NA
      mu = meanEstimation(na.omit(rettt))


      # Covariance estimation
      Sigma = covEstimation(na.omit(rettt))
      CR    = cor(na.omit(rettt))
      # Semi-deviation estimation
      semiDev = semidevEstimation(na.omit(rettt))
      #for robust Markowitz
      block_n = as.integer(5)
      n = length(mu)
      if (qr(Sigma)$rank < n) Sigma = covEstimation(na.omit(rettt), control = list(type = 'cor'))
      
      if (stratloop == "EW")  w1 = rep(1/dim(rettt)[2], dim(rettt)[2])
      if (stratloop == "GMV")  w1 = optimalPortfolio(Sigma = Sigma,  
                                                     control = list(type = 'minvol', constraint = 'none'))
      if (stratloop == "GMV_long")  w1 = optimalPortfolio(Sigma = Sigma,  
                                                          control = list(type = 'minvol', constraint = 'lo'))
      if (stratloop == "GMV_robust") {
        mrkw = import("markowitz") 
        gammas = mrkw$recommended_gamma(Sigma) 
        w1 = mrkw$robust_qp(rettt, block_n, gamma = gammas, lmbd = 0)}
      if (stratloop == "GMV_lin") {          
        Sigma_lin = linshrink_cov(na.omit(rettt), k = 0)
        w1 = optimalPortfolio(Sigma = Sigma_lin, 
                              control = list(type = 'minvol'))}
      if (stratloop == "GMV_nlin") {
        Sigma_nlin = nlshrink_cov(na.omit(rettt), k = 0)
        w1 = optimalPortfolio(Sigma = Sigma_nlin, 
                              control = list(type = 'minvol'))
      } 
      if (stratloop == "GMV_mcd") {Sigma_mcd = assetsMeanCov(na.omit(rettt), method = "mcd") 
      w1 = optimalPortfolio(Sigma = Sigma_mcd, mu = mu, control = list(type = 'minvol'))}
      if (stratloop == "GMV_mve")  {Sigma_mve = assetsMeanCov(na.omit(rettt), method = "mve")
      w1 = optimalPortfolio(Sigma = Sigma_mve, mu = mu, control = list(type = 'minvol'))}
      if (stratloop == "GMV_ogk")  {Sigma_ogk = assetsMeanCov(na.omit(rettt), method = "ogk")
      w1 = optimalPortfolio(Sigma = Sigma_ogk, mu = mu, control = list(type = 'minvol'))}
      if (stratloop == "GMV_nnve")  {Sigma_nnve = assetsMeanCov(na.omit(rettt), method = "nnve")
      w1 = optimalPortfolio(Sigma = Sigma_nnve, mu = mu, control = list(type = 'minvol'))}
      if (stratloop=="ERC"){w1 = Weights(PERC(Sigma = Sigma, par = NULL, 
                                              percentage = FALSE))}
     
      if (stratloop=="MD"){w1 = as.vector(PMD(na.omit(rettt))@weights)/100}
      # Inverse volatility portfolio
      if (stratloop == "invvol")  {w1 = round(optimalPortfolio(Sigma = Sigma,
                                       control = list(type = 'invvol')), 4)}
      
      # Equal-risk-contribution portfolio with the long-only constraint
      if (stratloop == "ERC ")    {w1 =  round(optimalPortfolio(Sigma = Sigma,
                                    control = list(type = 'erc', constraint = 'lo')), 4)}
      
      # Maximum diversification portoflio with the long-only constraint
     # maxdiv  = round(optimalPortfolio(Sigma = Sigma,
      #                                 control = list(type = 'maxdiv', constraint = 'lo')), 4)
      
      # Risk-efficient portfolio 
      if (stratloop == "riskeff")  {w1 = round(optimalPortfolio(Sigma = Sigma, semiDev = semiDev,
                                       control = list(type = 'riskeff')), 4)}

      
      # Maximum decorrelation portoflio with the long-only constraint
      maxdec  = round(optimalPortfolio(Sigma = Sigma,
                                       control = list(type = 'maxdec', constraint = 'lo')), 4)

# HRP DePrado
if (stratloop=="HRP"){
  w1 = getHRP(Sigma,CR,max=maxweight, min=minweight, robust_cov = FALSE, return_raw = rettt)
}
weight[t,] = w1 
    } 
    else{    
      weight[t,] = weight[(t - 1),]
    }
  }
  res                                   = as.vector(row_sums(weight*ret)) 
  rescoll[[stratloop]]                  = res;
  weightscoll[[stratloop]]              = weight; 
  print(list(returns_str = res, weights = weight)) 
}
calculation_time = proc.time() - ptm
weightscoll                       = results[seq(2,length(results), 2)]
names(weightscoll)                = strats
rescoll                           = as.data.frame(results[seq(1,length(results), 2)])
rownames(rescoll)                 = index(weightscoll[[1]])
colnames(rescoll)                 = strats 
#### Compute performance####
#long period
collNumbers                            = NA
collNumbers_tc                         = NA
cp                                     = NA
weightsNumbers                         = NA
start                                  = first_signal
end                                    = dim(ret)[1] 
cp                                     = calculatePerformanceMeasures(start,end);
collNumbers                            = cp[[1]];
collres                                = cp[[2]]; 
weightsNumbers                         = cp[[3]];
collNumbers_tc                         = cp[[4]]; 
collres_tc                             = cp[[5]];
colnames(collNumbers)                  = strats; 
rownames(collNumbers_tc)                  = c("CumWealth", "SR", "TTO", "TO")

collNumbers_tc                            = collNumbers
colnames(collNumbers_tc)               = strats; 

rownames(collNumbers_tc)              = c("CumWealth-TC", 
                                         "Sharpe.Ann-TC",  "TTO", "TO",
                                         "AverageDrawdown ",
                                         "StdDev.annualized",
                                         "CVaR",
                                         "AverageLength",
                                         "AverageRecovery ",
                                         "BernardoLedoitRatio ",
                                         "BurkeRatio ",
                                         "CDD ",
                                         "DownsideDeviation ",
                                         "DownsideFrequency ",
                                         "DownsidePotential ",
                                         "DRatio ",
                                         "HurstIndex ",
                                         "KellyRatio",
                                         "MartinRatio ",
                                         "OmegaSharpeRatio ",
                                         "ProspectRatio ",
                                         "SemiDeviation ",
                                         "SkewnessKurtosisRatio ",
                                         "SortinoRatio ",
                                         
                                         "UpsidePotentialRatio ",
                                         "UpsideRisk ",
                                         "VolatilitySkewness")

colnames(collres)                      = strats 
n                                      = length(strats)
cols                                   = rainbow(length(strats), s = 1, v = 1,
                                          start = 0, end = max(1, n - 1)/n, alpha = 1)#c("red","blue","pink","darkgreen","lightgreen","lightblue","black")
weightsNumbers = t(weightsNumbers)
colnames(weightsNumbers) = c("min", "max", "sd", "mad-ew", "max-min");


rownames(weightsNumbers) = strats#[b][-4];
#### Visualization of performace####

  ind = match(as.vector(c("EW",  "GMV_long", 
                          "GMV_lin", "GMV_nlin", "GMV",   
                          "ERC",      "MD" ,"GMV_robust" ,"HRP" , "invvol")), strats) #GMV_robust  "HRP"


write.csv(collNumbers_tc, "indicatori.csv")
number <- nlevels(strats)
my_colors = RColorBrewer::brewer.pal(length(strats),"Paired")
# count the needed levels of a factor


# repeat the given colors enough times
palette <- rep(c("color1", "color2"), length.out = number)
# EW - "#E41A1C" RColorBrewer::brewer.pal(7,"Set1")[1]
# GMV - "#377EB8" RColorBrewer::brewer.pal(7,"Set1")[2]
# GMV_lin - "#4DAF4A" RColorBrewer::brewer.pal(7,"Set1")[3]
# GMV_long - "#984EA3" RColorBrewer::brewer.pal(7,"Set1")[4]
# GMV_nlin - "#FF7F00" RColorBrewer::brewer.pal(7,"Set1")[5]
# GMV_robust - "#A65628" RColorBrewer::brewer.pal(7,"Set1")[7]

#IV "#A65628"
#ERC "#F781BF"
#MD "#999999"

write.csv(collNumbers_tc,"results.csv")
cumWealth = 1 + cumsum(collres[,ind])
write.csv(cumWealth, "cumWealth.csv")
cumWealth_graph=tidy(cumWealth)
cumWealth_graph$`Portfolio strategies`=cumWealth_graph$series
plot13 = cumWealth_graph %>% ggplot(aes(x = index, y = value, color = `Portfolio strategies`)) + 
  geom_line() + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", 
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent")) + # get rid of legend panel bg#

  #, legend.position = "none"
  ylab("Cumulative Wealth") + xlab("Time") + 
  scale_color_manual(values = my_colors) 

pdfname =  paste0("CumWealth", df, ".pdf")
pngname =  paste0("CumWealth", df, ".png")
pdf(file = pdfname)
plot13
dev.off()
png(file = pngname,  width=200,height=100,units='mm', res=300)
plot13
dev.off()


plot_13 = tidy(cumWealth)%>% ggplot(aes(x = index, y = value, color = series)) + #
  theme_bw() + scale_color_manual(values = my_colors) + geom_line() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", 
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent")) + # get rid of legend panel bg#
  ylab("Cumulative Wealth") + xlab("Time")
plot_13

# asset_classes = vector(mode = "character", 102)
# asset_classes =c(rep("cryptos",35), rep("XR",13),
#                  rep("bonds",7),  rep("stocks",27),rep("com",18),
#                  rep("real estate",2))
#
asset_classes = Assets_classes$Tickers
#which(asset_classes=="cryptos")1:28
#which(asset_classes=="stocks")43:61
#which(asset_classes=="bonds")62
#which(asset_classes=="real estate")
#which(asset_classes=="com")78:88
#which(asset_classes=="XR")65:77
#### Weights vizualization ####
GMV_lin = as.data.frame(weightscoll$GMV_lin[first_signal:end])
GMV_lin_assets = GMV_lin
colnames(GMV_lin_assets) = c(asset_classes)
GMV_lin_assets = as.data.frame(t(rowsum(t(GMV_lin_assets), 
                                        group = colnames(GMV_lin_assets), na.rm = T)))
GMV_lin$Date = as.Date(index(weightscoll$GMV_lin[first_signal:end]))
GMV_lin_assets$Date = as.Date(index(weightscoll$GMV_lin[first_signal:end]))
weights_GMV_lin = gather(GMV_lin, Ticker, P.weight, -Date)
weights_GMV_lin_assets = gather(GMV_lin_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot6 = ggplot(weights_GMV_lin, aes(x=Date, y=P.weight, fill=Ticker)) +
  geom_area() + theme_bw() + 
  ggtitle("GMV_lin") + theme(legend.position = "none", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)
plot6
plot6_A = ggplot(weights_GMV_lin_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("GMV_lin") + theme(legend.position = "right", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)
plot6_A
GMV_nlin = as.data.frame(weightscoll$GMV_nlin[first_signal:end])
GMV_nlin_assets = GMV_nlin
colnames(GMV_nlin_assets) = c(asset_classes)
GMV_nlin_assets = as.data.frame(t(rowsum(t(GMV_nlin_assets), 
                                         group = colnames(GMV_nlin_assets), na.rm = T)))
GMV_nlin$Date = as.Date(index(weightscoll$GMV_nlin[first_signal:end]))
GMV_nlin_assets$Date = as.Date(index(weightscoll$GMV_nlin[first_signal:end]))
weights_GMV_nlin = gather(GMV_nlin, Ticker, P.weight, -Date)
weights_GMV_nlin_assets = gather(GMV_nlin_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot7 = ggplot(weights_GMV_nlin, aes(x=Date, y=P.weight, fill=Ticker)) +
  geom_area() +
  theme_bw() + 
  ggtitle("GMV_nlin") + theme(legend.position = "none", 
                              plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[5])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)

plot7_A = ggplot(weights_GMV_nlin_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("GMV_nlin") + theme(legend.position = "right", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)


GMV_long = as.data.frame(weightscoll$GMV_long[first_signal:end])
GMV_long_assets = GMV_long
colnames(GMV_long_assets) = c(asset_classes)
GMV_long_assets = as.data.frame(t(rowsum(t(GMV_long_assets), 
                                         group = colnames(GMV_long_assets), na.rm = T)))
GMV_long$Date = as.Date(index(weightscoll$GMV_long[first_signal:end]))
GMV_long_assets$Date = as.Date(index(weightscoll$GMV_long[first_signal:end]))
weights_GMV_long = gather(GMV_long, Ticker, P.weight, -Date)
weights_GMV_long_assets = gather(GMV_long_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot9 = ggplot(weights_GMV_long, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +
  theme_bw() + ggtitle("GMV_long") +
  theme(legend.position = "none", plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[4])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)
plot9_A = ggplot(weights_GMV_long_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("GMV_long") + theme(legend.position = "right", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)


GMV_robust = as.data.frame(weightscoll$GMV_robust[first_signal:end])
GMV_robust_assets = GMV_robust
colnames(GMV_robust_assets) = c(asset_classes)
GMV_robust_assets = as.data.frame(t(rowsum(t(GMV_robust_assets), 
                                           group = colnames(GMV_robust_assets), na.rm = T)))
GMV_robust$Date = as.Date(index(weightscoll$GMV_robust[first_signal:end]))
GMV_robust_assets$Date = as.Date(index(weightscoll$GMV_robust[first_signal:end]))
weights_GMV_robust = gather(GMV_robust, Ticker, P.weight, -Date)
weights_GMV_robust_assets = gather(GMV_robust_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot10 = ggplot(weights_GMV_robust, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +
  theme_bw() + ggtitle("GMV_robust") + 
  theme(legend.position = "none", plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[7])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)
plot10_A = ggplot(weights_GMV_robust_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("GMV_robust") + theme(legend.position = "right", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)


EW = as.data.frame(weightscoll$EW[first_signal:end])
EW_assets = EW
colnames(EW_assets) = c(asset_classes)
EW_assets = as.data.frame(t(rowsum(t(EW_assets), 
                                   group = colnames(EW_assets), na.rm = T)))
EW$Date = as.Date(index(weightscoll$EW[first_signal:end]))
EW_assets$Date = as.Date(index(weightscoll$EW[first_signal:end]))
weights_EW = gather(EW, Ticker, P.weight, -Date)
weights_EW_assets = gather(EW_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot11 = ggplot(weights_EW, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +   ggtitle("EW") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[1])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)
plot11_A = ggplot(weights_EW_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("EW") + theme(legend.position = "right", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)


GMV = as.data.frame(weightscoll$GMV[first_signal:end])
GMV_assets = GMV
colnames(GMV_assets) = c(asset_classes)
GMV_assets = as.data.frame(t(rowsum(t(GMV_assets), 
                                    group = colnames(GMV_assets), na.rm = T)))
GMV$Date = as.Date(index(weightscoll$GMV[first_signal:end]))
GMV_assets$Date = as.Date(index(weightscoll$GMV[first_signal:end]))
weights_GMV = gather(GMV, Ticker, P.weight, -Date)
weights_GMV_assets = gather(GMV_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot12 = ggplot(weights_GMV, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +  ggtitle("GMV") + theme_bw() +
  theme(legend.position = "none", plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[2])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)

plot12_A = ggplot(weights_GMV_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("GMV") + theme(legend.position = "right", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)
ERC = as.data.frame(weightscoll$ERC[first_signal:end])
ERC_assets = ERC
colnames(ERC_assets) = c(asset_classes)
ERC_assets = as.data.frame(t(rowsum(t(ERC_assets), 
                                    group = colnames(ERC_assets), na.rm = T)))
ERC$Date = as.Date(index(weightscoll$ERC[first_signal:end]))
ERC_assets$Date = as.Date(index(weightscoll$ERC[first_signal:end]))
weights_ERC = gather(ERC, Ticker, P.weight, -Date)
weights_ERC_assets = gather(ERC_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot14 = ggplot(weights_ERC, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +  ggtitle("ERC") + theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)

plot14_A = ggplot(weights_ERC_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("ERC") + theme(legend.position = "right") +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)

MD = as.data.frame(weightscoll$MD[first_signal:end])
MD_assets = MD
colnames(MD_assets) = c(asset_classes)
MD_assets = as.data.frame(t(rowsum(t(MD_assets), 
                                   group = colnames(MD_assets), na.rm = T)))
MD$Date = as.Date(index(weightscoll$MD[first_signal:end]))
MD_assets$Date = as.Date(index(weightscoll$MD[first_signal:end]))
weights_MD = gather(MD, Ticker, P.weight, -Date)
weights_MD_assets = gather(MD_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot15 = ggplot(weights_MD, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +  ggtitle("MD") + theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)

plot15_A = ggplot(weights_MD_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("MD") + theme(legend.position = "right") +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)

IV = as.data.frame(weightscoll$invvol[first_signal:end])
IV_assets = IV
colnames(IV_assets) = c(asset_classes)
IV_assets = as.data.frame(t(rowsum(t(IV_assets), 
                                   group = colnames(IV_assets), na.rm = T)))
IV$Date = as.Date(index(weightscoll$invvol[first_signal:end]))
IV_assets$Date = as.Date(index(weightscoll$invvol[first_signal:end]))
weights_IV = gather(IV, Ticker, P.weight, -Date)
weights_IV_assets = gather(IV_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot16 = ggplot(weights_IV, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +  ggtitle("IV") + theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)

plot16_A = ggplot(weights_IV_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("MD") + theme(legend.position = "right") +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)

#HRP


HRP = as.data.frame(weightscoll$HRP[first_signal:end])
HRP_assets = HRP
colnames(HRP_assets) = c(asset_classes)
HRP_assets = as.data.frame(t(rowsum(t(HRP_assets), 
                                        group = colnames(HRP_assets), na.rm = T)))
HRP$Date = as.Date(index(weightscoll$HRP[first_signal:end]))
HRP_assets$Date = as.Date(index(weightscoll$HRP[first_signal:end]))
weights_HRP = gather(HRP, Ticker, P.weight, -Date)
weights_HRP_assets = gather(HRP_assets, Asset, P.weight, -Date)## , 1:dim(rets)[2])
plot17 = ggplot(weights_HRP, aes(x=Date, y=P.weight, fill=Ticker)) +
  geom_area() + theme_bw() + 
  ggtitle("HRP") + theme(legend.position = "none", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)
plot17
plot17_A = ggplot(weights_HRP_assets, aes(x=Date, y=P.weight, fill=Asset)) +
  geom_area() + theme_bw() + 
  ggtitle("HRP") + theme(legend.position = "right", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") #+ ylim(-4, 4)



png(file="weights.png", width = 465,
    height = 225, units='mm', res = 300)
    g=wrap_plots(plot16, plot17,
                 plot6,plot7,  plot9, plot10, plot11,
                 plot12, plot14, plot15, nrow = 5)
    g
    dev.off()


    
    png(file="weights_all.png",width = 465,
        height = 225, units='mm', res = 300)
    g_A=wrap_plots(plot16_A, plot17_A,
                 plot6_A,plot7_A,  plot9_A, plot10_A, plot11_A,
                 plot12_A, plot14_A, plot15_A, nrow = 5)
    g_A
    dev.off()
    
    


#### Create latex table ####
performance.table.latex = xtable(collNumbers[ind,], digits = 4)
print(x = performance.table.latex, file = paste0("Performance_",  df, ".tex") ,
      include.rownames = T, booktabs = T, floating = F)

performance_tc.table.latex = xtable(collNumbers_tc, digits = 4)
print(x = performance_tc.table.latex, file = paste0("Performance_tc_", df, ".tex") ,
      include.rownames = T, booktabs = T, floating = F)

weights.table.latex = xtable(weightsNumbers[ind,], digits = 4)
print(x = weights.table.latex, file = paste0("Weights_", df, ".tex") ,
      include.rownames = T, booktabs = T, floating = F)

####Save  Results ####

save.image(file = paste0(Sys.Date(),"Results_RobustM",  df,  ".RData"))  
#Tracking of packages used versions#
writeLines(capture.output(sessionInfo()), paste("sessionInfo.txt"))

####PCA
z <- t(collNumbers_tc)
str(z)



methods=c("EW",   "GMV_long", "GMV",
          "GMV_lin", "GMV_nlin", "invvol", "ERC", "MD","GMV_robust","HRP")
means <- apply(z,2,mean)
sds <- apply(z,2,sd)
nor <- scale(z,center=means,scale=sds)
distance = dist(nor)


#Correlation heatmap
cormat <- round(cor(nor ),2)
head(cormat)

melted_cormat <- melt(cormat)
head(melted_cormat)

# Heatmap

Var1=colnames(z)
Var2=colnames(z)


png("Correlation matrix.png", width=200,height=200,units='mm', res=300)


corr<-ggplot(data = melted_cormat, aes(X1, X2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()+
  theme(axis.title.x=element_blank(),  axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", 
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        text = element_text(size=20))
corr

dev.off()

png("Silhouette.png", width=200,height=100,units='mm', res=300)
fviz_nbclust(nor, kmeans, method = "silhouette", k.max = 8) + theme_minimal() 
dev.off()
###Cluster Dendrogram
mydata.hclust = hclust(distance)
plot(mydata.hclust)
plot(mydata.hclust,labels=methods, main='Default from hclust')

mydata.hclust<-hclust(distance,method="ward") 
plot(mydata.hclust,hang=-1) 


png("dendrogram.png", width=225,height=225,units='mm',res=300)


dev.off()



#Scree plot

png("Scree_plot.png")
###Scree Plot
wss <- (nrow(nor)-1)*sum(apply(nor,2,var))
for (i in 2:8) wss[i] <- sum(kmeans(nor, centers=i)$withinss)
g<-plot(1:8, wss, type="b", xlab="Number of Clusters",
        ylab="Within groups sum of squares")
g
dev.off()

####PCA
library("FactoMineR")
library(factoextra)
res.pca <- PCA(nor, graph = FALSE)

eig.val <- get_eigenvalue(res.pca)
eig.val

var <- get_pca_var(res.pca)

# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
var$contrib

library("corrplot")
corrplot(var$cos2, is.corr=FALSE)


res.desc <- dimdesc(res.pca, axes = c(1,2,3), proba = 0.05)
# Description of dimension 1
dim1<-res.desc$Dim.1
dim1=as.data.frame((dim1))
# Description of dimension 2
dim2<-res.desc$Dim.2
dim2=as.data.frame((dim2))
# Description of dimension 3
dim3<-res.desc$Dim.3
dim3=as.data.frame((dim3))
#### Create latex table ####

print(xtable(dim1, type = "latex"), file = "dim1.tex")

print(xtable(dim2, type = "latex"), file = "dim2.tex")

print(xtable(dim3, type = "latex"), file = "dim3.tex")
################################

png("Scree plot.png")
# plot(
#   x = seq(1:10), y = eigval[1:10],
#   type = "o",
#   xlab = "Princiale Component", ylab = "Eigenvalue", main="Scree plot")
y = res.pca$eig[1:6]
names(y)=c("1","2","3","4","5","6")

pareto.chart(y, xlab="Principal Component",
             ylab = "Eigenvalue",ylab2 = "Cumulative Variance",main="Pareto Chart")
dev.off()


fviz_pca_ind(res.pca)

fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)


k2 <- kmeans(nor, centers = 5, nstart = 1)

k2

png("clusters.png", width=225,height=150,units='mm', res=300)
fviz_cluster(k2, data = nor,ellipse = TRUE,ellipse.type = "convex",
             ellipse.level = 0.99,
             ellipse.alpha = 0.2,
             shape = NULL,
             pointsize = 1.5,
             labelsize = 12,
             xlab = "Wealth factor",
             ylab = "Risk factor",
             main=NULL, ggtheme = theme_minimal())
dev.off()

#Loadings
loadings=res.pca$var$cor

df <- data.frame(loadings[,1], loadings[,2], loadings[,3])

names(df) <- c("1", "2", "3")
rownames(df) <- colnames(z)
df
df2=melt(t(df), id.vars = rownames(df))

names(df2)[1] <- "Factor"
names(df2)[2] <- "Variable"


# Reverse the order for ggplot
df2$Variable <- factor(df2$Variable, levels = rev(levels(df2$Variable)))

png("Loadings.png",  width=225,height=150,units='mm', res=300)
ggplot(df2,aes(x =Factor ,y=Variable,fill=value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high ="red", mid="white",
                       midpoint = 0, limit = c(-1,1)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 0, vjust = 1, 
                                   size =8, hjust = 1))
dev.off()




####################################
# Factors and Kernel Density Contours
####################################
factors=res.pca$ind$coord

F <- data.frame(factors[,1], factors[,2], factors[,3])
colnames(F)=c("F1","F2","F3")
png("Factors1.png")
p1<-ggplot( F, aes( -F1, F2))+  geom_point(aes(color = as.factor(k2$cluster)),
                                           size = 4,
                                          show.legend = FALSE) +
  scale_color_manual(values = c("red", "green", "blue","black", "brown"))+
  geom_text(aes(label=methods),size=6)+xlim(-10,5)+
  labs(x="Wealth factor", y="Risk factor")
p1

dev.off()

k2$cluster
avg<-distances %>%
  group_by(method) %>%
  summarise(
    Mean=mean(distance, na.rm=TRUE),
    Median=median(distance, na.rm=TRUE),
    SD = sd(distance, na.rm = TRUE),
    Min=min(distance, na.rm=TRUE),
    Max=max(distance, na.rm=TRUE))


