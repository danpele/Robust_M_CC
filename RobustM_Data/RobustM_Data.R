##### SETTINGS #####

# Clean the environment 
graphics.off()
rm(list = ls(all = TRUE))

# Load Functions and other Files
source('./RobustM_Packages.R')
source('./RobustM_Functions.R')
Assets = list()
Assetsprices = as.matrix(read.csv("Assets_2022.csv", 
      sep = ",", dec = ".", row.names = 1, header = T, check.names = FALSE))
Assetsprices[,(36:101)] = Assetsprices[,(36:101)] %>% na.locf(fromLast=F)
Assetsdates = as.Date(rownames(Assetsprices), "%Y-%m-%d")
Assets$Data = xts::as.xts(Assetsprices,order.by=Assetsdates)
Assets$Data = Assets$Data[ , colSums(is.na(Assets$Data)) == 0]
#Calculate discrete/log returns
Assets$Rets    = CalculateReturns(Assets$Data)[-1]
Assets$LogRets = CalculateReturns(Assets$Data, method = "log")[-1]
Assets$LogRets = Assets$LogRets[ , colSums(is.na(Assets$LogRets )) == 0]#due to negative prices we remove two observations
# Relative returns, i.e. p1/p0
Assets$RelRets = na.omit(Assets$Rets + 1) 
Assets$Tickers = colnames(Assets$LogRets)
save(Assets, file = "Assets.RData")

Assets$Tickers
Assets_classes = list()
Assetsprices = as.matrix(read.csv("Assets_2022.csv", 
                                  sep = ",", dec = ".", row.names = 1, header = T, check.names = FALSE))

Assetsdates = as.Date(rownames(Assetsprices), "%Y-%m-%d")
Assetsprices[,(36:101)] = Assetsprices[,(36:101)] %>% na.locf(fromLast=F)
Assets_classes$Data = xts::as.xts(Assetsprices,order.by=Assetsdates)
colnames(Assets_classes$Data) = c(rep("cryptos",35), rep("bonds",7),
                                  rep("com",18),
                                  rep("XR",13), rep("real estate",2),
                                  rep("stocks",27)
                         )

Assets_classes$Data = Assets_classes$Data[ , colSums(is.na(Assets_classes$Data)) == 0]
#Calculate discrete/log returns
Assets_classes$Rets    = CalculateReturns(Assets_classes$Data)[-1]
Assets_classes$LogRets = CalculateReturns(Assets_classes$Data, method = "log")[-1]
Assets_classes$LogRets = Assets_classes$LogRets[ , colSums(is.na(Assets_classes$LogRets )) == 0]#due to negative prices we remove two observations
# Relative returns, i.e. p1/p0
Assets_classes$RelRets = na.omit(Assets_classes$Rets + 1) 
Assets_classes$Tickers = colnames(Assets_classes$LogRets)
save(Assets_classes, file = "Assets_classes.RData")