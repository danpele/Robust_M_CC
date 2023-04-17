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
df = "Assets"#"Russell3000""SP100"
load("Assets.RData")
#load("SP100.RData")

#Choose the starting point of analysis  
start_point = "2016-01-04" #for the Russell3000 dataset
myData = switch(df,
                "Assets" = Assets,
                "SP100" = SP100)
insample = F
random   = T #if TRUE then it picks 600 random stocks for Russell3000 dataset

#### Python environment ####
use_condaenv("r-reticulate")
path_to_python = "C:\\Users\\dell\\anaconda3\\envs\\r-reticulate\\python.exe"

#define the path to Python
Sys.which("python")
Sys.setenv(RETICULATE_PYTHON = path_to_python)
use_python(path_to_python,  required = TRUE)
reticulate::py_config()



#Uncomment if set up the Python-R connection for the first time
#create a new environment if you run the code for the first time
# conda_create("r-reticulate", python_version = 3.9)
# #install Packages to the environment
# conda_install("r-reticulate", "scipy")
# conda_install("r-reticulate", "pandas")
# conda_install("r-reticulate", "matplotlib")
# conda_install("r-reticulate", "scikit-learn", pip = T)
# conda_install("r-reticulate", "statsmodels", pip = T)
# conda_install("r-reticulate", "seaborn")
# conda_install("r-reticulate", "cvxopt", pip = T)#from cvxopt import solvers, matrix, spmatrix
# import ("pprint")
# np     = import("numpy")
# pd     = import("pandas")
# mpl    = import("matplotlib")
# plt    = import("matplotlib.pyplot")
# cvxopt = import("cvxopt")
# mrkw   = import("markowitz")

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

  w_length = 1 # number of years of observations to estimate input parameters

if (insample) {
  {first_signal = endpoints(commonDateR,on = "years")[1]}
} else {first_signal = endpoints(commonDateR,on = "years")[w_length + 1]}
ep = endpoints(commonDateR,on = freq) # rebalancing dates
ep = ep[(ep >= first_signal) & (ep < (end))]
ret_data   = list()
ret_random = list()



    for (t in commonDateR[ep]) {#building data for estimation of parameters
      t = as.POSIXlt(as.Date(t)) 
      start_date = t
      start_date$year = t$year - w_length # date of begining of the estimation window
      ret_data[[as.character(as.Date(t))]] = ret[which(index(ret) > as.Date(start_date) & index(ret) <= as.Date(t)),]
    }

  strats = c("EW",   "GMV_long", "GMV",
             "GMV_lin", "GMV_nlin", "invvol", "ERC", "MD","GMV_robust","HRP")#, "HRP")#,"GMV_robust", "GMV_mcd", "GMV_mve", "GMV_ogk", "GMV_nnve")


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
  period = c(first_signal:(end))#(end - 1)
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

gc()
#### Compute performance####
#long period
collNumbers                            = NA
collNumbers_tc                         = NA
cp                                     = NA
weightsNumbers                         = NA
start                                  = first_signal
end                                    = dim(ret)[1] 

cp                                     = calculatePerformanceMeasures(start,end);
gc()
collNumbers                            = cp[[1]];
collres                                = cp[[2]]; 
weightsNumbers                         = cp[[3]];
collNumbers_tc                         = cp[[4]]; 
collres_tc                             = cp[[5]];
colnames(collNumbers)                  = strats; 
rownames(collNumbers_tc)                  = c("CumWealth", "SR", "TTO", "TO")

collNumbers_tc                            = collNumbers

collNumbers_tc                            = collNumbers
colnames(collNumbers_tc)               = strats; 
rownames(collNumbers_tc)               = c("CumWealth-TC", 
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
                                           #"Kappa ",
                                           "KellyRatio",
                                           "MartinRatio ",

                                           "OmegaSharpeRatio ",
                                           "ProspectRatio ",
                                           "SemiDeviation ",
                                           "SkewnessKurtosisRatio ",
                                           "SortinoRatio ",
                                           "UpsidePotentialRatio ",
                                           "UpsideRisk ",
                                           "VolatilitySkewness")#"pfInformationRatio")
collNumbers_tc                         = t(collNumbers_tc)
colnames(collres)                      = strats 
n                                      = length(strats)
cols                                   = rainbow(length(strats), s = 1, v = 1,
                                                 start = 0, end = max(1, n - 1)/n, alpha = 1)#c("red","blue","pink","darkgreen","lightgreen","lightblue","black")
weightsNumbers = t(weightsNumbers)
colnames(weightsNumbers) = c("min", "max", "sd", "mad-ew", "max-min");
rownames(weightsNumbers) = strats#[b][-4];

z <- collNumbers_tc
str(z)

df_total<-as.data.frame(collNumbers_tc)
df_total['time']=as.numeric(1)
df_total['method']=c("EW",  "GMV_long", "GMV","invvol",
                 "GMV_lin", "GMV_nlin",    
                 "ERC",      "MD" ,"GMV_robust" ,"HRP")
methods=c("EW",   "GMV_long", "GMV",
          "GMV_lin", "GMV_nlin", "invvol", "ERC", "MD","GMV_robust","HRP")


means <- apply(z,2,mean)
sds <- apply(z,2,sd)
zz <- scale(z,center=means,scale=sds)

####PCA
library("FactoMineR")
library(factoextra)

p1<- prcomp(z, scale = TRUE)
summary(p1)
f2=p1$rotation ## f2 contains the standardized scoring coefficients;
F=zz%*%f2 #F contains the final scores after the varimax rotation;

scores=as.data.frame(F)
p1<-ggplot( scores, aes( PC1, PC2))+  geom_point() + 
theme_classic()+
  geom_text(aes(label=methods))+
  labs(title = "2016-2020",x="Wealth factor", y="Risk-return factor")
p1


k2 <- kmeans(zz, centers = 5, nstart = 25)
fviz_cluster(k2, data = zz,ellipse = TRUE,ellipse.type = "convex",
             ellipse.level = 0.99,
             ellipse.alpha = 0.2,
             shape = NULL,
             pointsize = 1.5,
             labelsize = 12,
             xlab = "Wealth factor",
             ylab = "Risk-return factor",
             main=NULL, ggtheme = theme_minimal())




n_waves=dim(ret)[1]-first_signal# number of waves


collNumbers_tc                         = collNumbers_tc
colnames(collres)                      = strats 
data<-as.data.frame(collNumbers_tc)
data['time']=1
data['method']=(c("EW",  "GMV_long", "GMV","invvol",
               "GMV_lin", "GMV_nlin",    
               "ERC",      "MD" ,"GMV_robust" ,"HRP"))

f1=start



w=as.integer(end/20)-1
seq<-seq(0, w)
d_f <- data.frame()
et=2*248
while(et<(end-1)) {
for (t in seq)

{
  
  
  print(t)
  
  st   = f1+20*t
  et   =st+248
 

gc()
  cp_t            = calculatePerformanceMeasures(st,et);
  collNumbers_t                           = cp_t[[1]]
  
  colnames(collNumbers_t)                  = strats; 
  
  
  collNumbers_tc_t                         = collNumbers_t
  colnames(collNumbers_tc_t)               = strats; 
  rownames(collNumbers_tc_t)               = c("CumWealth-TC", 
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
                                             #"Kappa ",
                                             "KellyRatio",
                                             "MartinRatio ",
                                              "OmegaSharpeRatio ",
                                             "ProspectRatio ",
                                             "SemiDeviation ",
                                             "SkewnessKurtosisRatio ",
                                             "SortinoRatio ",
                                             "UpsidePotentialRatio ",
                                             "UpsideRisk ",
                                             "VolatilitySkewness")#"pfInformationRatio")
  collNumbers_tc_t                       = t(collNumbers_tc_t)
df_t<-as.data.frame(collNumbers_tc_t)
df_t['time']=as.numeric(t)
df_t['method']=c("EW",  "GMV_long", "GMV","invvol",
               "GMV_lin", "GMV_nlin",    
               "ERC",      "MD" ,"GMV_robust" ,"HRP")
#row.names(df_t)=df_t['method']

d_f<-rbind(d_f,df_t)

}}



methods=c("EW",   "GMV_long", "GMV",
          "GMV_lin", "GMV_nlin", "invvol", "ERC", "MD","GMV_robust","HRP")

#d_f<-d_f[-(1:10), , drop = FALSE]
distances <- data.frame()
seq1=seq(0:59)
for (i in (seq1))
{  z<-d_f[d_f$time==i,]
print(i)
z=z[,1:27]
print(z)
means <- apply(z,2,mean)
sds <- apply(z,2,sd)
zz <- scale(z,center=means,scale=sds)

p1<- prcomp(z, scale = TRUE)
summary(p1)
f2=p1$rotation ## f2 contains the standardized scoring coefficients;
F=zz%*%f2 #F contains the final scores after the varimax rotation;
scores=as.data.frame(F)
k <- kmeans(zz, centers = 5, nstart = 1)

k$cluster
# target point
t_x <- 0
t_y <- 0

scores$distance <- sqrt(((scores$PC1)^2) + ((scores$PC2)^2)) 
scores['time']=as.integer(i)

file_name = paste("figs\\Time_",i,".png",sep="")
png(file_name)
scores$method <- c("EW",  "GMV_long", "GMV","invvol",
                    "GMV_lin", "GMV_nlin",    
                    "ERC",      "MD" ,"GMV_robust" ,"HRP")
scores$color_robust<-ifelse(str_detect(scores$method,"GMV_robust") ,
                            "GMV_robust","other")
distances<-rbind(distances,scores)
ttl=paste("Time: ",i)
fig1<-ggplot( scores, aes( PC1, PC2))+  
  theme_classic()+
  geom_point(aes(color = as.factor(color_robust)),size = 5,
             show.legend = FALSE) +
  scale_color_manual(values = c("red", "green", "blue","black", "brown"))+
  geom_text(aes(label=methods))+
  labs(title = ttl,x="Wealth factor", y="Risk-return factor")  + xlim(-7,7)+
  ylim(-7,7)
print(fig1)
dev.off()
}

#rownames(distances)<-distances$method
my_colors = RColorBrewer::brewer.pal(length(methods),"Paired")

png(file = "distance.png",  width=200,height=100,units='mm', res=300)
plot_dist = distances %>% ggplot(aes(x = time, y = distance, color = color_robust)) + 
  geom_line() + theme_bw() + 


  ylab("Distance to origin") + xlab("Time")

plot_dist




dev.off()

avg<-distances %>%
  group_by(method) %>%
  summarise(
            Mean=mean(distance, na.rm=TRUE),
            Median=median(distance, na.rm=TRUE),
            SD = sd(distance, na.rm = TRUE),
            Min=min(distance, na.rm=TRUE),
            Max=max(distance, na.rm=TRUE))


print(xtable(avg, type = "latex"), file = "distances.tex")
