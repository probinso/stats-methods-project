source("./../setup.R")

library(caret)
library(memoise)
library(doMC)

registerDoMC(7)

makeTRAIN = memoise(function(drug)
  genes_cov_thresh(0.2) %>%
    train_by_drug(drug) %>% hotextend_subtypes %>%
    mutate_all(as.numeric) %>% sort_cor_target(0.3, "success") %>%
    mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))
)

models = lapply(drugs, function(drug) {
  TRAIN = makeTRAIN(drug)
  rfctrl = rfeControl(functions=caretFuncs, method='cv', number=4)
  fmodel = rfe(
    success ~ ., data=TRAIN,
    rfeControl = rfctrl, methods='rf', sizes=seq(3,30,4))
  features = predictors(fmodel)
  
  tctrl    = trainControl(method='cv', number=2, search='random')
  cgcmod   = train(
    success ~., TRAIN[,c("success",features)],
    trControl=tctrl, method='svmLinearWeights', tuneLength=300)
  
  cgcmod
})

get_yeild = function(models, data, CHECK=F) {
  lapply(names(models), function(drug) {
    features = models[[drug]][["coefnames"]]
    samples = rownames(data)
    yhat = predict(models[[drug]], data[,features])
    tump = cbind(
      cellline=samples,
      drug=rep(drug, each=length(samples)),
      value=(yhat %>% `==`('Y') %>% ifelse(0, 1))
    )
    if (CHECK)
      cbind(tump,Y = success_by_drug(drug))
    else
      tump
  })
}

checkyield = models %>% get_yeild(all_train_data, T)
checkyield

yield = models %>% get_yeild(all_test_data)
df = yield %>% Reduce(rbind, .) %>% data.frame

df$id = apply(
  df,
  function(r)
    mapping[
      mapping$drug==r[["drug"]] &
        mapping$cellline==r[["cellline"]],
      ]$id,
  MARGIN = 1)

df[c("id", "value")] %>% arrange(id) %>% bak_and_save("submit.csv")
