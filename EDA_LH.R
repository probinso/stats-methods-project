library(readr)
library(lattice)
library(plyr)
expression <- read_delim("C:/Users/Lawrence Hsu/kaggle-stats/expression.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
subtypes <- read_delim("C:/Users/Lawrence Hsu/kaggle-stats/subtypes.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)



#EDA, I wish to know if there particular discrepancy between the different tissue types
#if there were any we would have to change how Cross validation would have been done
expression<-t(expression)
expression1<-expression
#created a new dataset with the first column removed as those contained the names of the genes
expression1$`184A1`<-NULL
#similar distribution to one another looks normal distribution 
boxplot(expression1)
hist(expression1$`600MPE`)
#looks normal
summary(expression1)
#18632 genes and 38 cell cultures didn't she say there was 39?
#Update: whoever made this expression dataset is an idiot. fills the rows of the cell
#line column with gene names and a headless column if you preview via microsoft excel 
#SAME ERROR in training_set_answers.txt 
dim(expression1)
#how many of different types of cell tissue types,
count(subtypes,'subtype')



