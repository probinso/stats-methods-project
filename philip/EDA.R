library(tidyverse)
source("./multiplot.R")
library(purrr)
library(ggplot2)
library(reshape2)
library(matrixStats)

library(easyGgplot2)

library(foreach) # use dopar

#########################################################
# Helper Functions
#########################################################

draw_rownames <- function(.data) .data %>%
  do(mutate(.,"rownames"=rownames(.)))

'%&%' = function(x, y)paste0(x,y) # string concatination

part   = function(f, ...) function(X) f(X, ...) # parameter-partial evaluation
'%|%'  = function(f, g) function(X) X %>% f %>% g # function composition
'%T%'  = function(d, funcs) lapply(funcs, part(function(f, d) f(d), d)) # Tee
'%T>%' = function(l, f) lapply(l, f) # lapply after Tee (consistant syntax)

read.tsv = part(read.delim, sep='\t', row.names=1, header=T, check.names = F)

matrix_boxplot = as.data.frame %|% stack %|%
  function(x) ggplot(x) + geom_boxplot() +
  aes(x=ind, y=values) + ylim(c(2,14)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dropcols = function(df, cols) df[,!(colnames(df) %in% cols)]

colCoV = function(x) { # coeficient of variance
  y = x %>% t
  m = rowMeans(y)
  s = rowSds(y)
  s/m
}

pop_header_row = function (df) {
  # malformed data, encoded labels in top row
  df.new = df
  colnames(df.new) = df.new[1,]
  df.new[-1,]
}
#########################################################
# Expression Data
#########################################################
df = read.tsv("./../data/expression.txt") %>% t
TotalSamples = dim(df)[1]; TotalGenes = dim(df)[2]

sample_names = rownames(df)
gene_names   = colnames(df)

gene_by_name = function(gene_name) df[,gene_name]

# genes sorted by coeficient of variance
sort_cov = df %>% colCoV %>% sort(decreasing = T)

# how many have coeficient of variance greater than 0.5
THRESH = 0.3
sum(sort_cov > THRESH)

summary(sort_cov)
which.max(sort_cov) %>% names
which.min(sort_cov) %>% names %>% gene_by_name %>% summary

critical_genes = sort_cov %>% .[.>THRESH] %>% names

sort_cov %>%
  qplot(.) + geom_vline(xintercept = THRESH, col="red")
sort_cov %>% .[critical_genes] %>%
  qplot(.) + geom_vline(xintercept = THRESH, col="red")

COUNT = 15
critical_genes %T%
  list(head=part(head, COUNT), tail=part(tail, COUNT)) %T>%
  gene_by_name %T>% matrix_boxplot
.Last.value %>% multiplot(plotlist = ., cols = 2)

##################################################################
# read subtypes off of disk
subtypes = read.tsv("./../data/subtypes.txt")
# read and re-structure training_set_answers from disk
targets  = read.tsv("./../data/training_set_answers.txt") %>%
  cbind(., subtype=subtypes[rownames(.),])

targets %>% arrange(subtype)

ratio_success = function(df) dmap(df, function(x) sum(x) / length(x))

tump = targets %>%
  mutate(subtype = factor(subtype)) %>%
  slice_rows("subtype") %>% # chunk into groups by "subtype"
  ratio_success %>% # compute success of each drug/subtype
  t %>% pop_header_row %>% data.frame # data munging

ns = rownames(tump)
tump %>%
  mutate_all(as.character) %>%
  mutate_all(as.numeric) %>%
  cbind(ns) %>%
  melt %>%
  ggplot(data=., aes(x=ns, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge())  +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


##########################################
get_drug_data = function(genes, drugname) {
  df = genes %>% gene_by_name
  sample_names = rownames(df)[rownames(df) %in% rownames(targets)]

  df[sample_names,] %>%
  cbind(
    subtype=subtypes[sample_names,],
    success=targets[sample_names, drugname])
}

rownames(targets)

get_drug_model = function(drugname) {
  INFO = get_drug_data(drugname)
  mod = glm(success ~ .,
            data=INFO %>% data.frame,
            family = binomial(link="logit"))
  step(mod)
}

drugs = targets %>% dropcols(c("subtype")) %>% names

badcell = targets %>% dropcols(c("subtype")) %>% rowSums %>% which.min %>% names

get_sample_data = function(samplename)
  df[samplename,] %>% rbind(subtype=subtypes[samplename])

lapply(drugs, function(drugname)
critical_genes %>% head(20) %>%
  get_drug_data(drugname) %>% dropcols(cols=c("subtype")) %>% data.frame %>%
  draw_rownames %>% melt(id=c("rownames", "success")) %>% 
  ggplot2.stripchart(
    data=., xName='variable',yName='value',
    groupName='success',
    position=position_dodge(0.8),
    backgroundColor="white",
    groupColors=c('#999999','#E69F00'),
    stat="identity", addBoxplot = T
    ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(drugname)
) %>% multiplot(plotlist = ., cols = 4)



.f = function () {
"
penalizedSVM uses automatic feature selection
  we should use elastic net, as lasso runs into problems with a poor 
  sample/feature ratio
other option is randomforests

If we first gene-select, then force inclusion of the subtype, then we 
  can force the 

we will use bagging for a voting system

Color code dotplot against drug success for visualizing feature usefulness  


What makes HCC1428 different than other cells?

If there is a good indicator of difference with HCC1428, then we 
  should calculate probability that any medication will work and multiply it
  by our result.
"
}

