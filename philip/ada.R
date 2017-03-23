
setwd("~/git/kaggle-stats")
library(caret)
library(memoise)
library(doParallel)
source(file.path("setup.R"))

library(ggplot2)
source(file.path("multiplot.R"))


registerDoParallel(8)

makedf = memoise(function(drug) {
  genes_cov_thresh(0.2) %>%
    train_by_drug(drug) %>% hotextend_subtypes %>%
    mutate_all(as.numeric) %>% sort_cor_target(0.3, "success") %>%
    mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))
})

fmodels = lapply(
  drugs,
  function(drug) {
    print("start:feature:" %&% drug)
    df = makedf(drug)

    rfctrl  = rfeControl(
      functions=caretFuncs,
      method="repeatedcv", number=3, repeats=5)

    fmodel = rfe(
      success ~ ., data = df,
      rfeControl=rfctrl, method="rf", sizes=seq(3, 25, 4))
    print("stop:feature:" %&% drug)
    fmodel
  }
)


Adamodels =
  lapply(names(fmodels), function(drug) {
    print("start:training:" %&% drug)

    df = makedf(drug)
    features = predictors(fmodels[[drug]])
    
    control = trainControl(method="repeatedcv", number=3, repeats=5)

    rf_default = train(
      success ~ ., data=df,
      method="ada", trControl=control)
    print("stop:training:" %&% drug)
    rf_default
  })
names(Adamodels) = drugs

