setwd('~/Documents/Projects/CMI-PB-2ndChallenge/predictions2020+21/')


source('../scripts/libs.R')
source('../scripts/run_glmnet.R')
detach("package:MOFA2", unload=TRUE) # need to detach MOFA in order to use the predict function from glmnet

load('data/predictors.RData')
load('data/test_2020.RData')

# test predictions on 2020 data
# -----------------------------
# -----------------------------



## baseline model
## --------------

# use demographic data, as well as assay data that should be highly predictive of the prediction task

features = c("Day0.IgG_PT", 'Day0.Meta.year_of_birth')
model1 = run_glmnet(p = features, task.name = 'IgG_PT', alpha=1)
model1$model_cor

## integration model
## -----------------

# use assay values from Day 0, as this has already been shown to be the most predictive,
# as well as factors generated from MOFA integration

# build model on 2021 data

f = paste('Day0.Factor', 1:15, sep='')

# test all factors
p.int = c(features, f)
model2 = run_glmnet(p = p.int, task.name = 'IgG_PT', alpha=1)
model2$model_cor
model2$model_coef

# test all factors
p.int = c(f)
model3 = run_glmnet(p = p.int, task.name = 'IgG_PT', alpha=1)
model3$model_cor
model3$model_coef

# test IgG_PT
p.int = c(features, 'Day0.IgG1_PT')
p.int
model4 = run_glmnet(p = p.int, task.name = 'IgG_PT', alpha=1)
model4$model_cor
model4$model_coef

# test IgG_PT
p.int = c(features, 'Day0.IgG1_PT', f)
model5 = run_glmnet(p = p.int, task.name = 'IgG_PT', alpha=1)
model5$model_cor
model5$model_coef



### STILL TO DO
# fold change
# predictors.baseline$fc = predictors.baseline$IgG_PT/predictors.baseline$Day0.IgG_PT
# test_data.baseline$fc = test_data.baseline$I / test_data.baseline$Day0.IgG_PT
# p.int = c(features, 'Day0.IgG1_PT')
# model_fc = run_glmnet(p = p.int, task.name = 'fc', alpha=0)
# model_fc$model_cor


# save best models plus baseline
# ------------------------------
# save(model4, model1, model_fc, file='data/regression_models_IgG_PT.RData')


# predict on 2022 data
# -----------------------------
# -----------------------------

load('data/test_2022.RData')

p = names(model4$model_coef)
p = p[!p %in% '(Intercept)']
p
new_data = na.omit(test_data.baseline[,p])
dim(new_data)

preds<-data.frame(predict(model4$model,newx=as.matrix(new_data), s='lambda.min'))
preds
preds$rnk = rank(-preds$lambda.min)
preds$Subject.ID = as.character(rownames(preds))
preds

# fold change predictions

preds_fc<-data.frame(predict(model_fc$model,newx=as.matrix(new_data), s='lambda.min'))
preds_fc
preds_fc$rnk_fc = rank(-preds_fc$lambda.min)
preds_fc$Subject.ID = as.character(rownames(preds_fc))
preds_fc

# save predictions

samples_to_predict = read.delim2('../data/2ndChallengeSubmissionTemplate.tsv', sep='\t')
samples_to_predict$Subject.ID = as.character(samples_to_predict$Subject.ID)
head(samples_to_predict)

x= left_join(samples_to_predict[,c(1:3)], preds)
x= left_join(x, preds_fc, by = 'Subject.ID')

x
write.table(x, file='data/results_Mn.tsv', sep='\t', quote = F, row.names = T, col.names = NA)
