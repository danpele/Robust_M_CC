# Install and load packages
#install.packages("IntroCompFinR", repos="http://R-Forge.R-project.org")
libraries = c("matrixcalc", "MASS", "PortfolioAnalytics", "PerformanceAnalytics", "zoo",
              "plotly", "RiskPortfolios", "devtools", "PMwR","Jmisc",
               "FRAPO",  "R.utils", "ade4", "grDevices", "foreach", 
               "doParallel","StepwiseTest", "stringr",  "DescTools","DT",
             "knitr", "tidyr", "plotly", "rmarkdown", "gridExtra", "reticulate",
              "fPortfolio", "xtable", "DEoptim", "ROI", "ROI.plugin.quadprog", "ROI.plugin.glpk",
              "nlshrink", "psych", "quantmod", "RColorBrewer", "broom", 
             "egg", "ggpubr","nloptr", "fAssets", "fastcluster",
             "tidyverse",
             "patchwork", "ggdendro","qcc","factoextra",
             "reshape", "pracma") # "IntroCompFinR", 
              
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)


