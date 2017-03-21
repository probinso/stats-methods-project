source(file.path("..", "setup.R"))
library(caret)

# select genes whose coefficient of variation is greater than 0.3
target_genes = genes_cov_thresh(0.2) %>% genes_cor_thresh(0.375)
target_genes = genes_cov_thresh(0.3) %>% genes_cor_thresh(0.6)
length(target_genes)

target_genes %>% data_by_genes %>% matrix_boxplot

# gene frame
DRUG = drugs[2]
lapply(drugs, function(DRUG) {
gf = target_genes %>% train_by_drug(DRUG)

control = trainControl(method="repeatedcv", number=4, repeats=3)

ldaFit = train(
  success ~ ., method='lda', data=gf,
  preProcess = c('scale', 'center'), trControl = control
  )

ldaFit$results$Accuracy
})


