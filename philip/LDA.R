source(file.path("..", "setup.R"))
library(caret)

# select genes whose coefficient of variation is greater than 0.3
target_genes = genes_cov_thresh(0.3) %>% genes_cor_thresh(0.9)

# gene frame
DRUG = drugs[1]
gf = target_genes %>% train_by_drug(DRUG)
sf = subtype_by_drug(DRUG)

