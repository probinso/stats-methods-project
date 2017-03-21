#!/usr/bin/env R

library(caret)
library(matrixStats)
library(reshape2)
library(tidyverse)

################################################################

part   = function(f, ...) function(X) f(X, ...) # parameter-partial evaluation

'%ni%' = Negate("%in%")
'%&%'  = function(x, y)paste0(x,y) # string concatination

'%|%'  = function(f, g) function(X) X %>% f %>% g # function composition
'%T%'  = function(d, funcs) lapply(funcs, part(function(f, d) f(d), d)) # Tee
'%T>%' = function(l, f) lapply(l, f) # lapply after Tee (consistant syntax)

draw_rownames <- function(.data) .data %>%
  do(mutate(.,"rownames"=rownames(.)))

read.tsv = part(read.delim, sep='\t', row.names=1, header=T, check.names = F)

pop_header_row = function (df) {
  # malformed data, encoded labels in top row
  df.new = df
  colnames(df.new) = df.new[1,]
  df.new[-1,]
}

dropcols = function(df, cols) df[,colnames(df) %ni% cols]

matrix_boxplot = as.data.frame %|% stack %|%
  function(x) ggplot(x) + geom_boxplot() +
  aes(x=ind, y=values) + ylim(c(2,14)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


################################################################
BASE = normalizePath(file.path("..", "data"))

gene_data = read.tsv(file.path(BASE, "expression.txt")) %>%
  as.matrix %>% t

training = read.tsv(file.path(BASE, "training_set_answers.txt")) %>%
  `==`(1) %>% ifelse(T, F)

subtypes = read.tsv(file.path(BASE, "subtypes.txt"))
# read and re-structure training_set_answers from disk
targets  = read.tsv(file.path(BASE, "training_set_answers.txt"))

mapping = read.csv(file.path(BASE, "scoring_and_test_set_id_mappings.csv")) %>%
  data.frame


drugs = colnames(training)
names(drugs) = drugs

rowCov = function(mat) rowSds(mat) / rowMeans(mat)
colCoV = t %|% rowCov

gene_cov = gene_data %>% colCoV


################################################################
genes_cov_thresh =
  function(CVTR) gene_cov %>% subset(gene_cov > CVTR) %>% names

genes_cor_thresh =
  function(names, COTR) {
    names %>% train_by_genes %>% cor %>%
      findCorrelation(cutoff=COTR) %>%
      sort %>% `*`(-1) %>%
      names[.]
  }

#gene_cov %>% qplot(.) + geom_vline(xintercept = CVTR, col="red")

# target_genes  = gene_cov %>% genes_cov_thresh(CVTR)
train_samples = row.names(training)
test_samples  = rownames(gene_data)[rownames(gene_data) %ni% train_samples]

data_by_genes  = function(genes) gene_data[, genes] %>% data.frame
train_by_genes = function(genes) gene_data[train_samples, genes] %>% data.frame
test_by_genes  = function(genes)
  gene_data[test_samples, genes] %>% data.frame

all_test_data = gene_data[test_samples, ] %>% data.frame %>% hotextend_subtypes

#target_genes %>% train_by_genes
#target_genes %>% test_by_genes

success_by_drug = function(drug)
  training[,drug] %>% `==`(T) %>% ifelse('Y', 'N')

# subtype frame
subtype_by_drug = function(drug)
  subtypes %>%
  cbind(success=success_by_drug(drug)[rownames(subtypes)]) %>%
  remove_missing

hotextend_subtypes = function(df) {
  dummyVars(" ~ .", data = subtypes) %>%
    predict(., newdata = subtypes) %>%
    .[rownames(df),] %>% cbind(df)
}

train_by_drug = function(genes, drug)
  genes %>% train_by_genes %>%
  cbind(success=success_by_drug(drug))

