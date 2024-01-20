# here we will use the 2021 data as training data, and test models on the 2021 data
# then use this model to predict 2022 tasks

source('../scripts/libs.R')
source('../scripts/generic_prediction_functions.R')

MOFAobject = load_model('../data/MOFA_models/MOFA2_noScale_train2020_21.hdf5')
MOFAobject

# get predictors
# --------------

predictors = format_training_matrix(mofa_object = MOFAobject) # concatenate MOFA factors, clinical data and assay data ; code categorical variables to numeric
predictors[1:5, 1:5]
predictors[1:5, 5110:ncol(predictors)]

# extract day 0 
# --------------

predictors.baseline = predictors[predictors$Meta.timepoint == 0, ]
rownames(predictors.baseline) = predictors.baseline$Meta.subject_id
colnames(predictors.baseline) = paste('Day0.', colnames(predictors.baseline), sep='') 
dim(predictors.baseline)
predictors.baseline[1:5, 1:5]


save(predictors.baseline, file='data/predictors.RData')

