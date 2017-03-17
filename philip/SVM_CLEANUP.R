#!/usr/bin/env R

library(caret)
library(matrixStats)
library(tidyverse)

part   = function(f, ...) function(X) f(X, ...) # parameter-partial evaluation
'%|%'  = function(f, g) function(X) X %>% f %>% g # function composition
'%T%'  = function(d, funcs) lapply(funcs, part(function(f, d) f(d), d)) # Tee
'%T>%' = function(l, f) lapply(l, f) # lapply after Tee (consistant syntax)

read.tsv = part(read.delim, sep='\t', row.names=1, header=T, check.names = F)

gene_data = read.tsv("./../data/expression.txt") %>% as.matrix %>% t

training  = read.tsv("./../training_set_answers.txt")
training = ifelse(training==1, 'Y', 'N')

rowCov = function(mat) rowSds(mat) / rowMeans(mat)
colCoV = t %|% rowCov

gene_cov = gene_data %>% colCoV
gene_cov %>% qplot

target_genes  = gene_cov %>% subset(gene_cov > 0.3) %>% names
train_samples = row.names(training)

train_by_genes = function(genes) gene_data[train_samples, genes] %>% data.frame
test_by_genes  = function(genes)
  gene_data[!rownames(gene_data) %in% train_samples, genes] %>% data.frame

#target_genes %>% train_by_genes
#target_genes %>% test_by_genes

success_by_drug = function(drug) training[,drug]

train_by_drug  = function(genes, drug)
  genes %>% train_by_genes %>%
  cbind(success=success_by_drug(drug))

##################################################################
SAMPLE_DRUG = "Cisplatin"
df = target_genes %>% train_by_drug(SAMPLE_DRUG)

ctrl = rfeControl(functions=caretFuncs, method='cv',number=5, verbose=T)

svmprofile =
  rfe(x=train_by_genes(target_genes),
    y=success_by_drug(SAMPLE_DRUG),
    method='svmRadial',
    sizes=c(1:5,10,20,40,70,90)
  )

fitcontrol =
  trainControl(method="repeatedcv", number=4, repeats=3, search="grid")

pamprofile = train(success ~ ., data=df, method="pam", trControl=fitcontrol)

mtry<-7
rfprofile = train(
  success ~ .,
  data=df,
  trControl=fitcontrol,
  method='rf',
  tuneGrid=expand.grid(mtry=c(1:70)))

#creating a svm fit to predict
svmgrid<-expand.grid(
  C=c(1,10,100),
  sigma=c(0,1.2,1.5,1.7,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6))

fitcontrol<-trainControl(method="repeatedcv",number=5, repeats=3, search="grid")
svm.fit<-train(
  success ~ C1S + SLC34A2,
  data=df,
  tuneGrid=svmgrid,
  trControl=fitcontrol,
  method='svmRadial',
  verbose=TRUE)
svm.predict<-predict(svm.fit, newdata=target_genes %>% test_by_genes)
print(svm.predict)
