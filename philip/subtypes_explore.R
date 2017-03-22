source(file.path("..", "setup.R"))
library(caret)

melt(targets %>% draw_rownames %>% cbind(success_by_drug(DRUG)))

sf = targets %>% cbind(subtype=subtypes[rownames(targets),])

subtypes

library(dummies)



tl = unique(subtypes$subtype)
tl[1]

sf[sf$subtype==tl[t],][[tl[1]] = 1

#for (st in unique(subtypes$subtype)) {
 
lapply(tl, function(t) sf[sf$subtype==t,] %>% dropcols("subtype") %>% colMeans)
