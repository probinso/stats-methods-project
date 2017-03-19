source("./../setup.R")

SAMPLE_DRUG = "Cisplatin"
gene_cor = target_genes %>% train_by_genes %>% cor

gene_cor %>% as.vector %>% qplot

write(gene_cor, "genecor.csv", sep = ',')

COTR = 0.7
gene_cor %>% 
  findCorrelation(cutoff=COTR) %>% 
  sort %>% `*`(-1) %>%
  target_genes[.] %>% train_by_drug(SAMPLE_DRUG) %>% cor %>%
  as.vector %>% qplot(.) + geom_vline(xintercept = COTR, col="red")


ctrl = rfeControl(functions=caretFuncs, method='cv',number=5, verbose=T)

svmprofile =
  rfe(
    x=train_by_genes(target_genes),
    y=success_by_drug(drugs[3]),
    method='svmRadial',
    sizes=c(1:5,10,20,40,70,90)
  )

fitcontrol =
  trainControl(method="repeatedcv", number=4, repeats=3, search="grid")

df = target_genes %>% train_by_drug(SAMPLE_DRUG)

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
