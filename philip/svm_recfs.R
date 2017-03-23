setwd("~/git/kaggle-stats")
source(file.path("setup.R"))
library(caret)
library(memoise)
library(doMC)

library(easyGgplot2)
library(ggplot2)
source(file.path("multiplot.R"))


registerDoMC(7)

makedf = memoise(function(drug) {
  genes_cov_thresh(0.2) %>%
    train_by_drug(drug) %>% hotextend_subtypes %>%
    mutate_all(as.numeric) %>% sort_cor_target(0.3, "success") %>%
    mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))
})

fmodels = lapply(
  drugs,
  function(drug) {
    print("start:" %&% drug)
    df = makedf(drug)

    rfctrl  = rfeControl(
      functions=caretFuncs,
      method="repeatedcv", number=3, repeats=5)

    fmodel = rfe(
      success ~ ., data = df,
      rfeControl=rfctrl, method="rf", sizes=seq(3, 25, 4))
    print("stop:" %&% drug)
    fmodel
  }
)


plots = lapply(
  names(fmodels), 
  function(drug) {
    features = predictors(fmodels[[drug]])
    df = all_train_data[, features] %>%
      cbind(success=success_by_drug(drug)[train_samples]) %>% data.frame %>%
      draw_rownames %>% melt(id=c("rownames", "success")) %>% 
      ggplot2.stripchart(
        data=., xName='variable',yName='value',
        groupName='success', position=position_dodge(0.8),
        backgroundColor="white",
        groupColors=c('#999999','#E69F00'),
        stat="identity", addBoxplot = T
      ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position="none",
            axis.title.x=element_blank())+
      ggtitle(drug)
})

plots[1]
plots %>% multiplot(plotlist = ., cols = 4)


svmmodels = lapply(
  names(fmodels),
  function(drug) {
    df = makedf(drug)
    features = predictors(fmodels[[drug]])
  
    ctrl = trainControl(method="repeatedcv", number=3, repeats=5)
  
    model = train(
      success ~ ., df[,c(features, "success")],
      trControl = ctrl, method = "svmLinear")
    model
  })
names(svmmodels) = drugs
lapply(svmmodels, confusionMatrix)


sum(lapply(svmmodels, function(m) attr(m$finalModel, "nSV")) <=8)
get_SV_count = function(m) attr(m$finalModel, "nSV")

CMPsvmmodels = lapply( # cost random
  names(fmodels),
  function(drug) {
    df = makedf(drug)
    features = predictors(fmodels[[drug]])
    
    ctrl = trainControl(method="repeatedcv", number=3, repeats=5)

    model = train(
      success ~ ., df[,c(features, "success")],
      trControl = ctrl, method = "svmRadial")
    model
  })
names(CMPsvmmodels) = drugs


RFmodels =
  lapply(names(fmodels), function(drug) {
    
    df = makedf(drug)
    features = predictors(fmodels[[drug]])
    
    control = trainControl(method="repeatedcv", number=3, repeats=5)
    mtry = sqrt(ncol(df))
    tg = expand.grid(.mtry=mtry)
    rf_default = train(
      success ~ ., data=df,
      method="rf", metric="Accuracy", tuneGrid=tg, trControl=control, ntree=5001)
    rf_default
  })
names(RFmodels) = drugs


Adamodels =
  lapply(names(fmodels), function(drug) {
    
    df = makedf(drug)
    features = predictors(fmodels[[drug]])
    
    control = trainControl(method="repeatedcv", number=3, repeats=5)

    rf_default = train(
      success ~ ., data=df,
      method="ada", trControl=control)
    rf_default
  })
names(Adamodels) = drugs




lapply(RFmodels,  confusionMatrix)


cbind(svm=lapply(svmmodels, get_SV_count), CMP=lapply(CRsvmmodels, get_SV_count))

confusionMatrix(svmmodels[[1]])
confusionMatrix(CMPsvmmodels[[1]])

checkyield = svmmodels %>% get_yield(all_train_data, T)
checkyield

yield = svmmodels %>% get_yield(all_test_data)
df = Reduce(rbind, yield) %>% data.frame
df$id = apply(
  df,
  function(r) 
    mapping[mapping$drug==r[["drug"]] & mapping$cellline==r[["cellline"]],]$id,
  MARGIN = 1)

df[c("id", "value")] %>% arrange(id) %>% bak_and_save("submit.csv")
