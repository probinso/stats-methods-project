library(caret)
library(matrixStats)
expression<-read.delim("expression.txt", check.names=FALSE)
training<-read.delim("training_set_answers.txt", check.names=FALSE)
expression2<-as.matrix(expression)
expression1<-t(expression)

#pre-feature selection (unsupervised)
mean<-rowMeans(expression)
sd<-rowSds(expression2)
cv<-sd/mean
hist(cv)
cv<-as.data.frame(cv)
cvsub<-subset(cv,cv>0.3)


cvsub<-as.data.frame(cvsub)
cvname<-rownames(cvsub)
expression3<-expression[cvname,]
expression3<-t(expression3)
stuff<-row.names(training)
expression4<-expression3[stuff,]
expression4<-as.data.frame(expression4)
expression3<-as.matrix(expression3)
expression3<-t(expression3)
training$`CGC-11047`<-factor(training$`CGC-11047`)
test1<-training$`CGC-11047`
expression4$test1<-training$`CGC-11047`
expression4$drug1<-as.factor(test1)
is.factor(expression4$drug1)

expression5<-as.data.frame(expression5)
expression5$drug1<-test2
names(testset)<-tolower(names(testset))
#subsetting the cell lines that were used in the test set
#to predict them 
testset<-expression[,-which(names(expression) %in% stuff)]
testset<-t(testset)
testset<-as.data.frame(testset)
#reordering test items to match map
order1<-c('HCC1187','MCF7','MDAMB361','MDAMB231','BT549','600MPE','HCC1954','SKBR3','MCF10A','MCF12A','HCC3153','MDAMB157','LY2','AU565')
test1<-ifelse(test1==0,'N',"Y")
test1<-as.factor(test1)
testsetord<-testset[order1]
order2<-match(order1,row.names(testset))
testsetorder<-testset[order2,]
training[training==0]<-'N'
training[training==1]<-'Y'
coln<-colnames(training)
training[1:12]<-lapply(training[1:12],as.factor)

set.seed(100)
#this first part wi
#ctrl<-trainControl(method="repeatedcv", number=5, repeats=2)
#model<-train(test1~., data=expression4, method="svmRadial", preProcess="scale", trControl=ctrl)
#importance<-varImp(model, scale=FALSE)
#class(importance)
#test3<-rownames(importance$importance)[1:10]
#importance$model
#rownames(expression5)<-NULL

#SVM feature selection
# sets the method of how the algorithm should be run 
ctrl<-rfeControl(functions=caretFuncs, method='cv',number=5, verbose=TRUE)
set.seed(130)
#test2 had to be a worded factor zero or one won't work 
svmprofile<-rfe(expression4[1:287],test1, rfeControl=ctrl,method='svmRadial',sizes=c(1:5,10,20,40,70,90))
print(svmprofile)
svmprofile1<-rfe(expression4[1:287],training$Cisplatin,method='svmRadial',sizes=c(1:5,10,20,40,70,90))

pamprofile<-train(drug1~.
                  ,data=expression4
                  , method='pam'
                  ,trControl=fitcontrol
                  ) 
mtry<-7
rfprofile<-train(drug1~.
                 ,data=expression4
                 , trControl=fitcontrol
                 , method='rf'
                 , tuneGrid=expand.grid(mtry=c(1:70)))

#creating a svm fit to predict 
svmgrid<-expand.grid(C=c(1,10,100), sigma=c(0,1.2,1.5,1.7,1.8,2.0,2.2,2.4,2.6,2.8))
fitcontrol<-trainControl(method="repeatedcv",number=5, repeats=3, search="grid")
svm.fit<-train(drug1~C1S+SLC34A2, data=expression4, tuneGrid=svmgrid, trControl=fitcontrol, method='svmRadial', verbose=TRUE)
svm.predict<-predict(svm.fit, newdata=testset)
print(svm.predict)

#gactrl<-gafsControl(functions = caretGA, method = "cv",number=5, verbose=TRUE)
#obj<-gafs(x=expression5,y=test2, iters=2, gafsControl = gactrl,method="svmRadial")



#names(expression5)<-tolower(names(expression5)[1:287])
#test1<-as.numeric(test1)
#expression4<-as.data.frame(expression4)
#test1<-factor(test1)
#ctrl<-rfeControl(functions=rfFuncs,method='repeatedcv',repeats=5,verbose=FALSE)
#subset<-c(1:10, 15,20,30,40,50,60,80,100,120,160,180,200)
#rfProfile<-rfe(expression4, train, sizes=subset, rfeControl = ctrl)



#colnames(expression4)
#colnames(bbbDescr)
#test2<-as.factor(test1)
#test2<-factor(test2)


