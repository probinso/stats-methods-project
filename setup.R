#!/usr/bin/env R

library(caret)
library(matrixStats)
library(reshape2)
library(tidyverse)

library(devtools)
devtools::install_github("probinso/lutilities", force=T)
library(lutilities)

# read integer from terminal
readinteger = function(msg) as.integer(readline(prompt=msg %&% " "))

# read.tsv provides helper for reading files
#   all data files had the same layout
read.tsv = partial(
  read.delim, sep='\t',
  row.names=1, header=T, check.names = F)

# takes a matrix and generates a simple boxplot
matrix_boxplot = as.data.frame %|% stack %|%
  function(x) ggplot(x) + geom_boxplot() +
  aes(x=ind, y=values) + ylim(c(2,14)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


############################################################
# This is code to read original files in and convert the
#   data to expected types
# gene_data : expression measures
# training  : drug scores
# subtypes  : sample -> subtype lookup table
# drugs     : the drugs by name
# mapping   : id lookup table for kaggle submission
BASE <- normalizePath(file.path("~/git/kaggle-stats/"))

gene_data <- read.tsv(file.path(BASE, "data", "expression.txt")) %>%
  as.matrix %>% t

training <-
  read.tsv(file.path(BASE, "data", "training_set_answers.txt")) %>%
  `==`(1) %>% ifelse(T, F)

subtypes <- read.tsv(file.path(BASE, "data", "subtypes.txt"))

mapping <-
  read.csv(
    file.path(BASE, "data" ,"scoring_and_test_set_id_mappings.csv")) %>%
  data.frame

# allowing this list to be indexed by drug name will make looping 
#   easier later
drugs <- colnames(training)
names(drugs) <- drugs

############################################################
# data transformation and filtering functions, as described
#   by function name.
rowCov = function(mat) rowSds(mat) / rowMeans(mat)
colCoV = t %|% rowCov

gene_cov <- gene_data %>% colCoV

# grab gene-names with covariance minimum of threshold
genes_cov_thresh =
  function(CVTR) gene_cov %>% subset(gene_cov > CVTR) %>% names

# grab gene-names with maximum absolute correlation upper bound of threshold
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

train_samples <- row.names(training)
test_samples  <- rownames(gene_data)[rownames(gene_data) %ni% train_samples]

data_by_genes = function(genes) gene_data[, genes] %>% data.frame

train_by_genes =
  function(genes) gene_data[train_samples, genes] %>% data.frame

test_by_genes = function(genes)
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
# bak_and_save saves the contents of a table to disk at filename
#   without clobbering the most recent saved file of that same
#   name, to allow for verifying saved changes
bak_and_save = function(contents, filename) {
  file.copy(filename, to = filename %&% ".bak", overwrite = T)
  contents %>% write.table(filename, sep=",", row.names = F, quote=F)
}

