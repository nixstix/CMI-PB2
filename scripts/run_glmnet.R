# run regularised regression model
# CV is leave-one-out, across different lamda values
# NOTE: input data is standardised (scaled) as default
# can adjust alpha (alpha 1 = lasso, alpha 0 = ridge, between is elastic net)

run_glmnet = function(p=features, task.name='' ,alpha=1){
    
    predictors.baseline = predictors.baseline[!is.na(predictors.baseline[,task.name]),]
    
    predictors.rmNA = predictors.baseline[,p]
    predictors.rmNA = na.omit(predictors.rmNA)
    dim(predictors.rmNA)
    
    task = predictors.baseline[rownames(predictors.rmNA), task.name]
    
    is.na(task)
    
    
    model<-cv.glmnet(x=as.matrix(predictors.rmNA), 
                     lambda = NULL,
                     task, family='gaussian',
                     alpha=alpha,
                     nfolds=nrow(predictors.rmNA),
                     type.measure="mse", 
                     standardize = TRUE) 
    
    model_coef=coef(model, s = 'lambda.min')[coef(model, s = 'lambda.min')[,1]!= 0]
    names(model_coef)=rownames(coef(model, s = 'lambda.min'))[coef(model, s = 'lambda.min')[,1]!= 0]
    model_coef = model_coef[order(model_coef)]
    model_coef
    
    
    # make predictions on 2021 data
    
    preds<-data.frame(predict(model,newx=as.matrix(test_data.baseline[,p]),s='lambda.min'))
    preds = na.omit(preds)
    preds$rnk = rank(-preds$lambda.min)
    preds$actual = test_data.baseline[rownames(preds), task.name]
    
    
    
    cor(preds$rnk, preds$actual, method = 'spearman')
    c1 = cor.test(preds$rnk, preds$actual, method = 'spearman')
    
    model = list(model=model, model_coef=model_coef, model_cor=c1, preds=preds, nSamples=nrow(predictors.rmNA))
    return(model)
       
    
}
