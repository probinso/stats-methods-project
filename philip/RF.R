setwd("~/git/kaggle-stats/philip")
source(file.path("..", "setup.R"))
library(caret)

target_genes = genes_cov_thresh(0.4) #%>% genes_cor_thresh(0.9)
length(target_genes) 

models =
  lapply(drugs, function(drug) {
    
  drug = drugs[1]
  df = target_genes %>% train_by_drug(drug)
  #df = dummyVars(" ~ .", data = subtypes) %>%
  #  predict(., newdata = subtypes) %>%
  #  .[rownames(df),] %>% cbind(df)

  control = trainControl(method="repeatedcv", number=4, repeats=3)
  mtry = sqrt(ncol(df))
  tg = expand.grid(.mtry=mtry)
  rf_default = train(
    success ~ ., data=df,
    method="rf", metric="Accuracy", tuneGrid=tg, trControl=control, ntree=100)
  rf_default
  })


models[3]

yield = lapply(names(models), function(drug) {
  yhat = predict(models[[drug]], test_by_genes(target_genes))
  cbind(
    test_samples,
    rep(drug, each=length(test_samples)),
    yhat %>% `==`('Y') %>% ifelse(1, 0)
  )
})
  

df = Reduce(rbind, yield) %>% data.frame
colnames(df) = c("cellline", "drug", "value")
df$id = apply(
  df,
  function(r) 
    mapping[mapping$drug==r[["drug"]] & mapping$cellline==r[["cellline"]],]$id,
  MARGIN = 1)

df[c("id", "value")] %>% arrange(id) %>% write.table("rf.csv", sep=",",row.names = F, quote=F)

