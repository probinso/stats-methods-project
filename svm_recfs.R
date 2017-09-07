#!/usr/bin/env Rscript
############################################################
# `setup.R` provides(") useful data manipulation code.
source(file.path("setup.R"))

############################################################
# library imports
library(caret)
library(memoise)
library(doParallel)

library(devtools)
install_github("kassambara/easyGgplot2", force=T)
library(easyGgplot2)
library(ggplot2)
library(e1071)

library(lutilities)
############################################################
# Yield of all time expensive code blocks are saved out to
#   disk. This allows for easier and iterative interactive
#   workflow after executing each code sections

CORES = 6 #readinteger("Number of cores you can sacrifice:")
registerDoParallel(CORES)

############################################################
# Need for constistant general filtering, in order to limit
#   and standardize learning parameters. This permits
#   an easier basis for model comparison.
# Function takes in a `drug` and produces the dataset for
#   learning with described parameters. Memoised for speedy
#   re-use. Return value later refered to as `loosly filtered`.
#
# drug : name of drug
# CVTR : minimum threshold coeficient of variation
# COTR : minimum threshold correlation to drug success
makedf = memoise(function(drug, CVTR=0.2, COTR=0.3) {
  genes_cov_thresh(CVTR) %>%
    train_by_drug(drug) %>% hotextend_subtypes %>%
    mutate_all(as.numeric) %>% sort_cor_target(COTR, "success") %>%
    mutate(success=as.factor(ifelse(success==1, 'Y', 'N')))
})

############################################################
# Feature Elimination is performed using Random Forests
#   a list of feature models is generated indexed by
#   drugname these features are intended for use in SVM below
fmodels = lapply(
  drugs,
  function(drug) {
    print("start:features:" %&% drug)
    df = makedf(drug)

    tg = expand.grid(.mtry=sqrt(ncol(df)))
    rfctrl  = rfeControl(
      functions=caretFuncs,
      method="repeatedcv", number=3, repeats=8, allowParallel=T)

    fmodel = rfe(
      success ~ ., data = df,
      rfeControl=rfctrl, method="rf", sizes=seq(2, 25, 3), tuneGrid=tg)
    print("features:stop:" %&% drug)
    fmodel
  }
)
save(file="fmodels.bak", fmodels)

############################################################
# Models generated with random forests, using loosley
#   filtered parameters.
RFmodels =
  lapply(drugs, function(drug) {
    print("RF:start:" %&% drug)
    df = makedf(drug)

    control = trainControl(
      method="repeatedcv", number=3, repeats=5, classProbs = T)

    mtry = sqrt(ncol(df))
    tg = expand.grid(.mtry=mtry)
    rf_default = train(
      success ~ ., data=df,
      method="rf", metric="ROC",
      tuneGrid=tg, trControl=control, ntree=4000)
    print("RF:start:" %&% drug)
    rf_default
  })
save(file="RFmodels.bak", RFmodels)


############################################################
# Produces box-plots of model's most influential features.
#   It is expected that all models in one list will have the
#   same form, in order to support a consistant means of
#   extracting features.
#
# fmodels     : models to extract features from
# getfeatures : function that takes in a model and returns
#   its most influential features.
make_plots = function(fmodels, getfeatures) lapply(
  drugs,
  function(drug) {
    features = getfeatures(fmodels[[drug]])
    df = all_train_data[, features] %>%
      cbind(success=success_by_drug(drug)[train_samples]) %>%
      data.frame %>%
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
      ggtitle(drug)+
      ylim(2, 15)
  })


############################################################
# Save plots for each model in lists, in order to generate
#   side by side multiplots.
rfeplots = fmodels  %>% make_plots(function(model) predictors(model))
rfplots  = RFmodels %>% make_plots(function(model) {
  tump = varImp(model)$importance %>% draw_rownames %>%
    arrange(Overall) %>% tail(10)
  tump$rownames %>% rev
})


############################################################
# Multiplots for comparison are saved out to disk.
plots = lapply(drugs, function(drug) {
  png(drug %&% ".png")
  list(rfeplots[[drug]], rfplots[[drug]]) %>%
    multiplot(plotlist = ., cols = 2)
  dev.off()
  }
)

############################################################
# Because the SVM is given preference, we hope to include
#   correlated cell subtypes. We use the same criteria as
#   the `loosly filtered` features
target_features = lapply(
  drugs,
  function(drug) {
    df = makedf(drug)
    colnames(df)[grepl("subtype*", colnames(df))] %>%
    union(predictors(fmodels[[drug]])) %>% unlist
  })

############################################################
# Create svm-models indexed by drugname trained using the
#   features from yield of `recursive feature elimination`.
svmmodels = lapply(
  drugs,
  function(drug) {
    print("SVM:start:" %&% drug)
    df = makedf(drug)
    features = target_features[[drug]]

    ctrl = trainControl(
      method="repeatedcv", number=3, repeats=5,
      classProbs = T)

    model = train(
      success ~ ., df[,c(features, "success")],
      trControl = ctrl, method = "svmLinearWeights", metric="ROC")
    print("SVM:stop:" %&% drug)
    model
  })
save(file="svmmodels.bak", svmmodels)

############################################################
# get_yield produces a dataframe containing predicitons in
#   a consistant manner, given input data and a predictor
#   model.
get_yield = function(models, data, CHECK=F) {
  lapply(names(models), function(drug) {
    features = models[[drug]][["coefnames"]]
    samples  = rownames(data)
    yhat = predict(models[[drug]], data[,features])
    tump = cbind(
      cellline=samples,
      drug  = rep(drug, each=length(samples)),
      value =(yhat %>% `==`('Y') %>% ifelse(0, 1))
    )
    if (CHECK)
      cbind(tump, Expected = success_by_drug(drug))
    else {
      tump
    }
  })
}

############################################################
# Code verifies class assignment looks as expected
checkyield = svmmodels %>% get_yield(all_train_data, T)
checkyield

############################################################
# Function saves model yields out to disk
#
# models   : list of trained predictor models
# filename : filename to print classification yield to
save_yield = function(models, filename) {
  yield = models %>% get_yield(all_test_data)
  df = Reduce(rbind, yield) %>% data.frame
  df$id = apply(
    df,
    function(r)
      mapping[
        mapping$drug==r[["drug"]] & mapping$cellline==r[["cellline"]],
      ]$id,
    MARGIN = 1)

  df[c("id", "value")] %>% arrange(id) %>% bak_and_save(filename)
}

RFmodels  %>% save_yield("rf_submit.csv")
svmmodels %>% save_yield("svm_submit.csv")
