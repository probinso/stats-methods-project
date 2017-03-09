library(tidyverse)
source("./multiplot.R")
library(purrr)
library(ggplot2)
library(reshape2)


'%&%' = function(x, y)paste0(x,y) # string concatination

part   = function(f, ...) function(X) f(X, ...) # parameter-partial evaluation
'%|%'  = function(f, g) function(X) X %>% f %>% g # function composition
'%T%'  = function(d, funcs) lapply(funcs, part(function(f, d) f(d), d)) # Tee
'%T>%' = function(l, f) lapply(l, f) # lapply after Tee (consistant syntax)

g = function(i) {
  as = function(x) {
    filename = "myplot_" %&% i %&% ".pdf"
    ggsave(file = "./" %&% filename, plot = x)
    g(i+1);
  }
  as
}
autosave = g(1)

df = read.delim("./../data/expression.txt", sep="\t", header=T, check.names = F) %>% t
dim(df)

sample_names = rownames(df)
sample_names
gene_names   = colnames(df)

matrix_boxplot = as.data.frame %|% stack %|%
  function(x) ggplot(x) + geom_boxplot() +
  aes(x=ind, y=values) + ylim(c(2,14)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

colCoV = function(x) { # coeficient of variance
  y = x %>% t
  m = rowMeans(y)
  v = rowSums((y - m)^2)/(dim(y)[2] - 1)
  v / m
}

gene_by_name = function(gene_name) df[,gene_name]

sort_cov = df %>% colCoV %>% sort(decreasing = T)

sum(sort_cov > 0.5)

summary(sort_cov)
which.max(sort_cov)
which.min(sort_cov) %>% names %>% gene_by_name %>% summary

target_genes = sort_cov %>% .[.>0.5] %>% names

pl = sort_cov %>% .[.>0.5] %>% qplot(.)
pl



autosave = autosave(pl)

COUNT = 15

bplots =
  sort_cov %>% names %T%
  list(head=part(head, COUNT), tail=part(tail, COUNT)) %T>%
  gene_by_name %T>% matrix_boxplot

autosave = bplots[["head"]] %>% autosave
autosave = bplots[["tail"]] %>% autosave
bplots %>% multiplot(plotlist = ., cols = 2)

draw_rownames <- function(.data) .data %>% do(mutate(.,rownames=rownames(.)))

targets  = read.delim("./../data/training_set_answers.txt", row.names=1, header=T)
subtypes = read.delim("./../data/subtypes.txt", header=T, row.names=1)
targets  = targets %>% cbind(subtype=subtypes[rownames(targets),])
targets  = targets %>% mutate(subtype = factor(subtype)) %>% arrange(subtype)
tump = targets %>% slice_rows("subtype") %>% dmap(function(x) sum(x) / length(x)) %>% t
colnames(tump) = tump[1,]
tump = tump[-1,]
tump = data.frame(tump)
ns = rownames(tump)
tump %>% mutate_all(as.character) %>% mutate_all(as.numeric) %>% cbind(ns) %>%
  melt %>% ggplot(data=., aes(x=ns, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge())  +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))



##########################################3
INFO = df[,target_genes] %>% cbind(subtype=subtypes[rownames(df),])
INFO %>% data.frame

mod = glm(subtype ~ ., data=INFO %>% data.frame, family = binomial(link="logit"))
step(mod)
