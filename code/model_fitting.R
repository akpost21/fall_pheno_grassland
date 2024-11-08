#Run new models

#Install BayesianTools package first or phenor download fails
install.packages("BayesianTools")

#Upload akpost21 forked version of phenor
#Only need to do first time code is run
if(!require(devtools)){install.packages(devtools)}
devtools::install_github("akpost21/phenor")


#load packages
library(phenor)
library(ggplot2)
library(ggpubr)
library(beepr)
library(plotly)


#Upload datasets
setwd("G:/My Drive/fall_pheno_project/R code/Published code")

#Upload example dataset
#Downloaded and processed with R package "phenocamr"
All_1day <- readRDS("./AllSites_fall_1day_90_trans50_noOut_v4.rds")

#Only need code below if using original version of phenor
#If using akpost21 fork of phenor (see above), no need to run this

# #upload new fall models
# source("./New_fall_models_final.R")
# 
# #upload edited functions
# source("./pr_flatten_SM_SOS.R")
# source("./pr_fit_SM_SOS.R")

# #upload parameter file
# params_file <- "./fall_parameter_ranges_final.csv"


####Test single model####

test_model = pr_fit_SM_SOS(
  data = All_1day,
  model = "CDD_W3",
  method = "GenSA",
  control = list(max.call = 100),
  #par_ranges = file.path(params_file)
)

beep()

#See observed vs predicted of single model (combined data)
plot(test_model$measured, test_model$predicted, ylab = "Predicted", xlab = "Observed")
abline(lm(test_model$predicted ~ test_model$measured), col = "red") #regression line


####Run multiple models together####

#include names of all models you want to test in this list

models = c("CDD","CDDs","CDDP","CDD2","CDD3","CDDs3","SOS","W3","SM3","VPD3",
          "DD_W3","DD_SM3","DD_VPD3","CDD_W3","CDD_SM3","CDD_VPD3","CDD_DD_W3",
          "CDD_DD_SM3","CDD_DD_VPD3","W_T_bal3","CDDP3","SMP3","VPDP3","CDDP_DD_W3")


#specify dataset
data_site = All_1day

##Fit models
#max.call = # of iterations (increase to larger number ~100,000)
#specify parameter file

pr_fit_new <- function(model){
  Fit <- pr_fit_SM_SOS(model = model,
                data = data_site,
                method = "GenSA",
                control = list(max.call = 100), 
                #par_ranges = file.path(params_file)
                )
  return(Fit)
}

model_fits <- mapply(pr_fit_new, model= models, SIMPLIFY = FALSE)

beep()


####Extract AIC & RMSE####

#Extract elements (rmse = 5, AIC = 7) from nested lists ([[ = 2nd level list)
RMSE <- lapply(model_fits, '[[', 5)
AIC_list <- lapply(model_fits, '[[', 7)
AIC <- lapply(AIC_list, '[[', 1)

diag <- do.call(rbind, Map(data.frame, AIC = AIC, RMSE = RMSE))

#order models by lowest AIC
diag[order(diag$AIC),]


#####Assign location to each value & graph####

#Make list of all observed transition dates
All_trans_dates <-  do.call(rbind, Map(data.frame, site = lapply(data_site, '[', 1),
                                     observed = lapply(data_site, '[[', 5),
                                      year = lapply(data_site, '[[', 7)))

#Make data frame of all observed and predicted dates
#each = # of site-years of data

results <- do.call(rbind, Map(data.frame, model = rep(unlist(lapply(model_fits, '[', 1)), each = length(All_trans_dates$year)),
                              predicted = unlist(lapply(model_fits, '[[', 4)), 
                              observed = unlist(lapply(model_fits, '[[', 3))))


##Use automatic alphabetical order to add location names to results data frame

#number of sites
s <- data.frame(All_trans_dates$site, All_trans_dates$year)
#number of models
n <- length(diag$AIC)
#Replicate site names based on number of models
Site_names <- do.call("rbind", replicate(n, s, simplify = FALSE))


#Combine site names with results data frame, rename columns
results2 <- cbind(results, Site_names)
names(results2)[names(results2) == 'All_trans_dates.site'] <- 'site'
names(results2)[names(results2) == 'All_trans_dates.year'] <- 'year'

##See all model fits as regressions
All <- ggplot(data = results2, mapping = aes(x = observed, y = predicted ))+
  geom_point(aes(col = site))+
  geom_smooth(method= "lm", se = FALSE, fullrange = F, col = "black")+
  geom_abline(intercept = 0, slope = 1, col = "dark grey")+
  facet_wrap(facets = vars(model), scales = "free") +
  theme(aspect.ratio=1)+
  stat_cor(aes(label =..rr.label..))
All + labs(title="Model Fits", y="Predicted", x = "Observed")

#Make graphs interactive
ggplotly(All)



###Optional to save results: 

#save model fits
setwd("G:/My Drive/fall_pheno_project")
saveRDS(model_fits, file = "./output/model_fits/AllSites_model_output.rds")


#Save results table (RMSE & AIC)
write.csv(results2,"./output/model_fits/AllSites_model_output.csv", row.names = FALSE)



#Random extra
#Put elements in nested list in alphabetical order
All_1day_order <- All_1day[order(names(All_1day))]
