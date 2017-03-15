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
expression5<-expression3[stuff,]
testset<-expression[!stuff]
expression3<-as.matrix(expression3)
expression3<-t(expression3)
training$`CGC-11047`<-factor(training$`CGC-11047`)
test1<-training$`CGC-11047`
expression4$test1<-training$`CGC-11047`
expression4$test1<factor(expression4$test1)
is.factor(expression4$test1)
expression5<-as.data.frame(expression5)
expression5$drug1<-test2
names(testset)<-tolower(names(testset))
#subsetting the cell lines that were used in the test set
#to predict them 
testset<-expression[,-which(names(expression) %in% stuff)]
testset<-t(testset)
testset<-as.data.frame(testset)



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
ctrl<-rfeControl(functions=caretFuncs, method='repeatedcv',number=10,repeats = 1, verbose=TRUE)
set.seed(130)
#test2 had to be a worded factor zero or one won't work 
svmprofile<-rfe(expression5,test2, rfeControl=ctrl,method='svmRadial',sizes=c(1:5,10,20,40,70,90))

#creating a svm fit to predict 
fitcontrol<-trainControl(method="repeatedcv",number=10, repeats=2)
svm.fit<-train(drug1~c1s+slc34a2+flj21986+tmem45a+pi15 , data=expression5, trControl=fitcontrol, method='svmRadial', verbose=TRUE)
svm.predict<-predict(svm.fit, newdata=testset)


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


