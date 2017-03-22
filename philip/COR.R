setwd("~/git/kaggle-stats/philip")
source(file.path("..", "setup.R"))
library(caret)
library(glmnet)

####################################################
models =
  lapply(drugs, function(drug) {
    df = genes_cov_thresh(0.3) %>% train_by_drug(drug) %>% hotextend_subtypes %>% 
      mutate_all(as.numeric) %>% sort_cor_target(0.3, "success") %>% 
      mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))
    print(drug)
    print(dim(df))
    
    control = trainControl(method="repeatedcv", number=5, repeats=8)
    mtry = sqrt(ncol(df))
    tg = expand.grid(.mtry=mtry)
    rf_default = train(
      success ~ ., data=df,
      method="rf", metric="Accuracy", tuneGrid=tg, trControl=control, ntree=1001)
    rf_default
  })


lapply(models, confusionMatrix)

"CGC-11047", "Carboplatin", "GSK1070916", "PF-3084014", "PF-4691502"

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

df = Reduce(rbind, yield) %>% data.frame
df$id = apply(
  df,
  function(r) 
    mapping[mapping$drug==r[["drug"]] & mapping$cellline==r[["cellline"]],]$id,
  MARGIN = 1)

df[c("id", "value")] %>% arrange(id) %>% bak_and_save("rf.csv")

file.copy("rf.csv", to = "rf.csv.bak", overwrite = T)
df[c("id", "value")] %>% arrange(id) %>% write.table("rf.csv", sep=",",row.names = F, quote=F)

