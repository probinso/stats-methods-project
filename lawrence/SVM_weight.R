subtypes= c('subtype.Basal','subtype.Claudin.low','subtype.Luminal','subtype.Normal.like')
`CGC-11047`


rfctrl= rfeControl(functions=caretFuncs,method='cv',number =4)
feats=rfe(trsf, training$`CGC-11047`, rfeControl = rfctrl, methods='rf', sample="up",sizes=seq(3,30,4))
features<-predictors(feats)
tctrl=trainControl(method='cv',number=4, search='random')
cgcmod=train(CGC11047~.,final[,c("CGC11047",features)],trControl=tctrl, method='svmLinearWeights', sample='up', tuneLength=300)
cgcmod
confusionMatrix(cgcmod)
cgcmod$finalModel
