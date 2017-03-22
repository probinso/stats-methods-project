setwd("~/git/kaggle-stats/philip")
source(file.path("..", "setup.R"))
library(caret)
library(memoise)
library(doMC)

registerDoMC(8)

makedf = memoise(function(drug) {
  genes_cov_thresh(0.2) %>%
    train_by_drug(drug) %>% hotextend_subtypes %>%
    mutate_all(as.numeric) %>% sort_cor_target(0.3, "success") %>%
    mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))
})

fmodels = lapply(
  drugs[
    c("CGC-11047", "Carboplatin", 
      "GSK1070916", "PF-3084014", 
      "PF-4691502")],
  function(drug) {
    df = makedf(drug)

    rfctrl  = rfeControl(
      functions=caretFuncs,
      method="repeatedcv", number=5, repeats=8)

    fmodel = rfe(
      x=df %>% dropcols(c("success")), y=df[,"success"],
      rfeControl=rfctrl, method="rf", sizes=seq(3, 30, 5))

    fmodel
  }
)

svmmodels = lapply(
  names(fmodels),
  function(drug) {
    df = makedf(drug)
    features = predictors(fmodels[[drug]])
  
    ctrl = trainControl(method="repeatedcv", number=5, repeats=8, sampling="up")
  
    model = train(success ~ ., df[,c(features, "success")], trControl = ctrl, method = "svmLinear")
    model
  })
names(svmmodels) = c("CGC-11047", "Carboplatin", 
                     "GSK1070916", "PF-3084014", 
                     "PF-4691502")
lapply(svmmodels, confusionMatrix)


confusionMatrix(models[["PF-4691502"]])

replace_models = function(src, repl) {
  dst = src
  for (drug in names(repl)) dst[[drug]] = repl[[drug]]
  dst
}

names()

finalset = replace_models(
  models,
  svmmodels[names(svmmodels) %ni% c("PF-4691502")]
)

yield = lapply(names(finalset), function(drug) {
  features = finalset[[drug]][["coefnames"]]
  features
  
  features[features %ni% colnames(all_test_data)]
  yhat = predict(finalset[[drug]], all_test_data[,features])
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

