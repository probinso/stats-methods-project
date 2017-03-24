setwd("~/git/kaggle-stats/")
source(file.path("setup.R"))
library(caret)
library(doParallel)

registerDoParallel(7)


RFmodels =
  lapply(drugs, function(drug) {
    
    df = genes_cov_thresh(0.2) %>%
      train_by_drug(drug) %>% hotextend_subtypes %>%
      mutate_all(as.numeric) %>% sort_cor_target(0.3, "success") %>%
      mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))

    control = trainControl(method="repeatedcv", number=3, repeats=5)
    mtry = sqrt(ncol(df))
    tg = expand.grid(.mtry=mtry)
    rf_default = train(
      success ~ ., data=df,
      method="rf", metric="Accuracy", tuneGrid=tg, trControl=control, ntree=4000)
    rf_default
  })

varImp(RFmodels[[1]]) 

%>% head(10)

lapply(RFmodels, function(m) m[["results"]][["Accuracy"]]) %>% unlist %>% sort

checkyield = RFmodels %>% get_yield(all_train_data, T)
checkyield

yield = RFmodels %>% get_yield(all_test_data)
df = Reduce(rbind, yield) %>% data.frame
df$id = apply(
  df,
  function(r) 
    mapping[mapping$drug==r[["drug"]] & mapping$cellline==r[["cellline"]],]$id,
  MARGIN = 1)

df[c("id", "value")] %>% arrange(id) %>% bak_and_save("submit.csv")
