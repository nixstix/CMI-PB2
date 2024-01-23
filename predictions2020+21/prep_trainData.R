

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

# extract tasks
# -------------

# Freq Mn at Day 1
p = predictors[predictors$Meta.timepoint == 1, ]
rownames(p) = p$Meta.subject_id
predictors.baseline$Freq.Monocytes = p[rownames(predictors.baseline), 'Freq.Monocytes']

# CCL3 at Day 3
p = predictors[predictors$Meta.timepoint == 3, ]
rownames(p) = p$Meta.subject_id
predictors.baseline$CCL3 = p[rownames(predictors.baseline), 'Expr.ENSG00000277632.1']


# IgG_PT at Day 14
p = predictors[predictors$Meta.timepoint == 14, ]
rownames(p) = p$Meta.subject_id
predictors.baseline$IgG_PT = p[rownames(predictors.baseline), 'IgG_PT']



# save
# ----

save(predictors.baseline, file='data/predictors.RData')

