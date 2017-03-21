source(file.path("..", "setup.R"))
library(caret)

# select genes whose coefficient of variation is greater than 0.3
target_genes = genes_cov_thresh(0.2) %>% genes_cor_thresh(0.375)
length(target_genes)

# gene frame
DRUG = drugs[2]
lapply(drugs, function(DRUG) {
gf = target_genes %>% train_by_drug(DRUG)
#sf = subtype_by_drug(DRUG)

control = trainControl(method="repeatedcv", number=4, repeats=3)

ldaFit = train(
  success ~ ., method='qda', data=gf,
  preProcess = c('scale', 'center'), trControl = control
  )

ldaFit$results$Accuracy
})


rbind(lapply(drugs, success_by_drug) %>% unlist)

rbind(targets success_by_drug(DRUG)
      