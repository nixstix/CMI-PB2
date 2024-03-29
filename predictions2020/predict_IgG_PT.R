setwd('~/Documents/Projects/CMI-PB-2ndChallenge/predictions2020/')


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

# select top features from above
p.int = c(features,'Day0.Factor8', 'Day0.Factor13', 'Day0.Factor6')
model3 = run_glmnet(p = p.int, task.name = 'IgG_PT', alpha=1)
model3$model_cor
model3$model_coef


# different alphas
p.int = c(features,'Day0.Factor8', 'Day0.Factor13', 'Day0.Factor6')
model4 = run_glmnet(p = p.int, task.name = 'IgG_PT', alpha=0.5)
model4$model_cor
model4$model_coef

# different alphas
p.int = c(features,'Day0.Factor8', 'Day0.Factor13', 'Day0.Factor6')
model5 = run_glmnet(p = p.int, task.name = 'IgG_PT', alpha=0)
model5$model_cor
model5$model_coef


# fold change
predictors.baseline$fc = predictors.baseline$IgG_PT/predictors.baseline$Day0.IgG_PT
test_data.baseline$fc = test_data.baseline$IgG_PT / test_data.baseline$IgG_PT
#p.fc =c(features)
#p.fc = c(features, f)
#p.fc = c(features,'Day0.IgG1_PT', 'Day0.Factor8')
p.fc = c("Day0.IgG_PT","Day0.Meta.year_of_birth", "Day0.Factor9")

model_fc = run_glmnet(p = p.fc, task.name = 'fc', alpha=0)
model_fc$model_cor
model_fc$preds

# save best models plus baseline
# ------------------------------
save(model4, model1, model_fc, file='data/regression_models_IgG_PT.RData')


# predict on 2022 data
# -----------------------------
# -----------------------------

load('data/test_2022.RData')

#p = names(model4$model_coef)
#p = p[!p %in% '(Intercept)']
#p

p = c("Day0.IgG_PT", 'Day0.Meta.year_of_birth')
new_data = na.omit(test_data.baseline[,p])
dim(new_data)

preds<-data.frame(predict(model1$model,newx=as.matrix(new_data), s='lambda.min'))
preds
preds$rnk = rank(-preds$lambda.min)
preds$Subject.ID = as.character(rownames(preds))
preds

# fold change predictions

new_data = na.omit(test_data.baseline[,p.fc])
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
write.table(x, file='data/results_IgG_PT.tsv', sep='\t', quote = F, row.names = T, col.names = NA)
