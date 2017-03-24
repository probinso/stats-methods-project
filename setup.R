#!/usr/bin/env R

library(caret)
library(matrixStats)
library(reshape2)
library(tidyverse)

############################################################
# Helper functions providing infix operators for cleaner code

# parameter-partial evaluation
part   = function(f, ...) function(X) f(X, ...)

# 'x' not in 'y'
'%ni%' = Negate("%in%")

# string concatination
'%&%'  = function(x, y)paste0(x,y)

# function composition
'%|%'  = function(f, g) function(X) X %>% f %>% g # function composition

# T pipe, takes data to multiple functions
'%T%'  = function(d, funcs)
  lapply(funcs, part(function(f, d) f(d), d))

# T map, takes multple datasets to one function
#   equivelent to lapply (using T consistant syntax)
'%T>%' = function(l, f) lapply(l, f)

# to fix tibble removal of rownames, this saves the rownames
#   as a column
draw_rownames <- function(.data) .data %>%
  do(mutate(.,"rownames"=rownames(.)))

# read.tsv provides helper for reading files
#   all data files had the same layout
read.tsv = part(
  read.delim, sep='\t',
  row.names=1, header=T, check.names = F)

# remove specific columns from a dataframe
dropcols = function(df, cols) df[,colnames(df) %ni% cols]

# takes a matrix and generates a simple boxplot
matrix_boxplot = as.data.frame %|% stack %|%
  function(x) ggplot(x) + geom_boxplot() +
  aes(x=ind, y=values) + ylim(c(2,14)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# read integer from terminal
readinteger = function(msg) as.integer(readline(prompt=msg %&% " "))

############################################################
# This is code to read original files in and convert the
#   data to expected types
# gene_data : expression measures
# training  : drug scores
# subtypes  : sample -> subtype lookup table
# drugs     : the drugs by name
# mapping   : id lookup table for kaggle submission
BASE = normalizePath(file.path("."))

gene_data = read.tsv(file.path(BASE, "expression.txt")) %>%
  as.matrix %>% t

training =
  read.tsv(file.path(BASE, "training_set_answers.txt")) %>%
  `==`(1) %>% ifelse(T, F)

subtypes = read.tsv(file.path(BASE, "subtypes.txt"))

mapping =
  read.csv(
    file.path(BASE, "scoring_and_test_set_id_mappings.csv")) %>%
  data.frame

drugs = colnames(training)
names(drugs) = drugs

############################################################
# data transformation and filtering functions, as described
#   by function name.
rowCov = function(mat) rowSds(mat) / rowMeans(mat)
colCoV = t %|% rowCov

gene_cov = gene_data %>% colCoV

genes_cov_thresh =
  function(CVTR) gene_cov %>% subset(gene_cov > CVTR) %>% names

genes_cor_thresh =
  function(names, COTR) {
    names %>% train_by_genes %>% cor %>%
      findCorrelation(cutoff=COTR) %>%
      sort %>% `*`(-1) %>%
      names[.]
  }

sort_cor_target = function(df, COTR, target)
  apply(df, 2, function(col)cor(col, df[,target])) %>%
  abs %>% sort(decreasing = T) %>% .[-c(1)] %>% .[.>COTR] %>%
  names %>% c(target) %>% df[, .]

train_samples = row.names(training)
test_samples  = rownames(gene_data)[rownames(gene_data) %ni% train_samples]

data_by_genes  = function(genes) gene_data[, genes] %>% data.frame

train_by_genes =
  function(genes) gene_data[train_samples, genes] %>% data.frame

test_by_genes  = function(genes)
  gene_data[test_samples, genes] %>% data.frame

success_by_drug = function(drug)
  training[,drug] %>% `==`(T) %>% ifelse('Y', 'N')

subtype_by_drug = function(drug)
  subtypes %>%
  cbind(success=success_by_drug(drug)[rownames(subtypes)]) %>%
  remove_missing

# Takes in a dataframe whose rows are celllines, and adds hot-encoded
#   subtype data to each row, while respecting row order of input frame
hotextend_subtypes = function(df) {
  dummyVars(" ~ .", data = subtypes) %>%
    predict(., newdata = subtypes) %>%
    .[rownames(df),] %>% cbind(df)
}

all_test_data =
  gene_data[test_samples, ] %>% 
  data.frame %>% 
  hotextend_subtypes

all_train_data =
  gene_data[train_samples, ] %>%
  data.frame %>%
  hotextend_subtypes

train_by_drug = function(genes, drug)
  genes %>% train_by_genes %>%
  cbind(success=success_by_drug(drug))

############################################################
# Functions to facilitate generating predictions for kaggle
#
# get_yield produces a dataframe containing predicitons in
#   a consistant manner, given input data and a predictor
#   model.
get_yield = function(models, data, CHECK=F) {
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
      cbind(tump,Expected = success_by_drug(drug))
    else {
      tump
    }
  })
}

# bak_and_save saves the contents of a table to disk at filename
#   without clobbering the most recent saved file of that same
#   name, to allow for verifying saved changes
bak_and_save = function(contents, filename) {
  file.copy(filename, to = filename %&% ".bak", overwrite = T)
  contents %>% write.table(filename, sep=",",row.names = F, quote=F)
}

