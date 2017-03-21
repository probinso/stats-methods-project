setwd("~/git/kaggle-stats/philip")
source(file.path("..", "setup.R"))
library(caret)
library(glmnet)

target_genes = genes_cov_thresh(1) #%>% genes_cor_thresh(0.7)
length(target_genes) 

x = target_genes %>% train_by_genes %>% hotextend_subtypes
cvfits = lapply(drugs, function(drug) {
  y = success_by_drug(drugs[1]) %>% `==`('Y') %>% ifelse(1, 0) 
  cvfit = cv.glmnet(x=as.matrix(x), y=y, family="binomial", alpha=0.5)
})


plot(cvfits[[drugs[1]]])


see  = part(coef, s = "lambda.min")
tump = lapply(cvfits, function(fitted) see(fitted) %>% as.matrix %>% .[. != 0,])
tump
