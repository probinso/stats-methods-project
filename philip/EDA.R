source("./../setup.R")
source("./../multiplot.R")

library(easyGgplot2)
library(ggplot2)


critical_genes = gene_cov %>% subset(gene_cov > CVTR) %>% names

gene_cov %>%
  qplot(.) + geom_vline(xintercept = CVTR, col="red")
gene_cov %>% .[critical_genes] %>%
  qplot(.) + geom_vline(xintercept = CVTR, col="red")

COUNT = 15
critical_genes %T%
  list(head=part(head, COUNT), tail=part(tail, COUNT)) %T>%
  data_by_genes %T>% matrix_boxplot
.Last.value %>% multiplot(plotlist = ., cols = 2)

##################################################################

ratio_success = function(df) dmap(df, function(x) sum(x) / length(x))

tump = targets %>% cbind(., subtype=subtypes[rownames(.),]) %>%
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
  df = genes %>% data_by_genes
  sample_names = rownames(df)[rownames(df) %in% rownames(targets)]

  df[sample_names,] %>%
  cbind(
    subtype=subtypes[sample_names,],
    success=targets[sample_names, drugname])
}

rownames(targets)

get_drug_model = function(drugname) {
  INFO = critical_genes %>% get_drug_data(drugname)
  mod = glm(success ~ .,
            data=INFO %>% data.frame,
            family = binomial(link="logit"))
  step(mod)
}

drugs[1] %>% get_drug_model

drugs = targets %>% dropcols(c("subtype")) %>% names

badcell = targets %>% dropcols(c("subtype")) %>% rowSums %>% which.min %>% names

get_sample_data = function(samplename)
  df[samplename,] %>% rbind(subtype=subtypes[samplename])


plots = 
  lapply(drugs,
         function(drugname)
           critical_genes %>% head(15) %>%
           get_drug_data(drugname) %>% dropcols(cols=c("subtype")) %>% 
           data.frame %>% draw_rownames %>% melt(id=c("rownames", "success")) %>% 
           ggplot2.stripchart(
             data=., xName='variable',yName='value',
             groupName='success', position=position_dodge(0.8),
             backgroundColor="white",
             groupColors=c('#999999','#E69F00'),
             stat="identity", addBoxplot = T
             ) +
           theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
           theme(legend.position="none",
                 axis.title.x=element_blank())+
           ggtitle(drugname)
         )

plots %>% multiplot(plotlist = ., cols = 4)



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

