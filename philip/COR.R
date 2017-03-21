setwd("~/git/kaggle-stats/philip")
source(file.path("..", "setup.R"))
library(caret)
library(glmnet)

#target_genes = genes_cov_thresh(0.2)
#length(target_genes)



sort_cor_target = function(df, COTR, target)
  apply(df, 2, function(col)cor(col, df[,target])) %>%
  abs %>% sort(decreasing = T) %>% .[-c(1)] %>% .[.>COTR] %>%
  names %>% c(target) %>% df[, .] %>% mutate

df = genes_cov_thresh(0.2) %>% train_by_drug(drugs[1]) %>% hotextend_subtypes %>% 
  mutate_all(as.numeric) %>% sort_cor_target(0.2, "success") %>%
  mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))


####################################################
models =
  lapply(drugs, function(drug) {
    df = genes_cov_thresh(0.2) %>% train_by_drug(drug) %>% hotextend_subtypes %>% 
      mutate_all(as.numeric) %>% sort_cor_target(0.2, "success") %>% 
      mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))
    
    control = trainControl(method="repeatedcv", number=3, repeats=4)
    mtry = sqrt(ncol(df))
    tg = expand.grid(.mtry=mtry)
    rf_default = train(
      success ~ ., data=df,
      method="rf", metric="Accuracy", tuneGrid=tg, trControl=control, ntree=300)
    rf_default
  })


lapply(models, confusionMatrix)

yield = lapply(names(models), function(drug) {
  features = models[[drug]][["coefnames"]]
  features
  
  features[features %ni% colnames(all_test_data)]
  yhat = predict(models[[drug]], all_test_data[,features])
  cbind(
    cellline=test_samples,
    drug=rep(drug, each=length(test_samples)),
    value=(yhat %>% `==`('Y') %>% ifelse(1, 0))
  )
})

yield[[2]]

df = Reduce(rbind, yield) %>% data.frame
df
df$id = apply(
  df,
  function(r) 
    mapping[mapping$drug==r[["drug"]] & mapping$cellline==r[["cellline"]],]$id,
  MARGIN = 1)

df[c("id", "value")] %>% arrange(id) %>% write.table("rf.csv", sep=",",row.names = F, quote=F)

