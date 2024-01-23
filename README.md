## CMI-PB-2ndChallenge

## Integration model


### Overview

We start with the harmonised and normalised data provided by the CMI-PB team:
`/data/processed_datasets/prediction_dataset/master_processed_prediction_data.RData` and
`/data/processed_datasets/training_dataset/master_processed_training_data.RData`.

For each assay type, this data has been normalised to the median value across all subjects at Day0 (baseline),
and batch-corrected across the 2020-2021 data. More information about this pre-processing may be found here:

https://www.cmi-pb.org/blog/understand-data/#Our%20approach%20to%20data%20standardization ,
with code available here : https://www.cmi-pb.org/blog/prediction-challenge-overview/#Data%20and%20resources

Additional R functions and libraries are loaded from the `scripts` folder

#### Pre-processing

The input data is first processed in such a way as to make it compatible with MOFA -- ensuring each assay (view) has the same number of samples (filled in with NA). We also only select 5000 genes (the most highly variable genes) to try to balance the gene expression data with the other assays that have fewer features (gene expression data is renormalised and where necessary batch-correct). We split the 2020+2021 dataset into training and testing data (65-35). These steps are available here: `pre-processing` for each of : 2020, 2021, 2022, 2020+2021-TRAIN, and 2020+2021-TEST

#### Running MOFA

Scripts to run MOFA may be found in `MOFA`. We integrated assay data from the 2021 dataset, as well as the 2020+2021-TRAIN dataset. The file `first_pass_noScale_train2021.Rmd` and its follow-up, `first_pass_noScale_train2021_explore.Rmd` are annotated to give more information about the pipeline.

We also tried running MOFA, scaling each view, however, performance often appeared worse after the scaling.

#### Model building and predictions

Models were built using two datasets: 2021 (then tested on 2020) ; 2020+2020 TRAIN (then trained on 2020+2021 TEST). Finally, predictions were run on the 2022 data. Scripts prepare the training and test data, build models, and run predictions are in `predictions2020` and `predictions2020+21`. Data associated with predictions are available in `predictions2020/data` and `predictions2020+21/data`

Models were built using LASSO regression. 



