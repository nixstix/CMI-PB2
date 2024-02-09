# REFORMATTING OF TRAINING DATA FOR INPUT INTO MOFA
# -------------------------------------------------

# This data has already been normalised and batch-corrected by the CMI-PB team. 
# 

data_to_process='2020_dataset'

## LOAD DATA

source('../scripts/libs.R')

dat = readRDS('../data/processed_datasets/training_dataset/master_processed_training_data.RDS')
names(dat)


## SELECT DS FROM 2020 AS TRAINING SET
dat$subject_specimen$specimen_id = as.character(dat$subject_specimen$specimen_id)
sampleIdx = dat$subject_specimen[dat$subject_specimen$dataset == data_to_process, 'specimen_id']

dat$subject_specimen = dat$subject_specimen[dat$subject_specimen$specimen_id %in% sampleIdx$specimen_id,]
table(dat$subject_specimen$dataset)
dim(dat$subject_specimen)


# MERGE DATA
## merge all matrices to replace empty values with NAs

meta = dat$subject_specimen
colnames(meta) = paste("Meta.", colnames(meta), sep='')

# ab
x = as.data.frame(t(dat$abtiter_wide$normalized_data))
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# cytokines
x = as.data.frame(t(dat$plasma_cytokine_concentrations$normalized_data))
colnames(x) = paste("Cytokine.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))


# cell_freq
x = as.data.frame(t(dat$pbmc_cell_frequency$normalized_data))
colnames(x) = paste("Freq.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# expression

## NOTE: in a normal workflow, the input for this section should not be batch-corrected or TPM data, but raw count data. If raw count data is not available, the batch-corrected data can be directly fed into MOFA without any kind of filtering or transformation.

counts = dat$pbmc_gene_expression$batchCorrected_data
counts = counts[, which(colnames(counts) %in% meta$Meta.specimen_id)]
counts = apply(counts,2, as.integer)
rownames(counts) = rownames( dat$pbmc_gene_expression$batchCorrected_data)
dim(counts)

# select fewer genes to process, based on most highly variable genes
# renormalise

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = dat$pbmc_gene_expression$metadata[colnames(counts) , ],
                              design = ~ infancy_vac)
dds
dds = estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)
dim(assay(vsd))


# mofa is being dominated by factor 1, this may be a result of poor normalisation, therefore let's limit the geneset even further, based on HVG
topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),5000)
topVarGenes <- rownames(vsd)[topVarGenes]

dds <- dds[topVarGenes,]
dds
dds = estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)
dim(assay(vsd))




x = as.data.frame(t(assay(vsd)))
colnames(x) = paste("Expr.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

#train_data = t(meta)



train_data = list(meta = meta[,grep('Meta', colnames(meta))],
                  geneExp = t(meta[,grep('Expr.', colnames(meta))]) ,
                  ab = t(meta[,grep('^IgG', colnames(meta))]) ,
                  cytokine = t(meta[,grep('Cytokine', colnames(meta))]), 
                  cell_freq = t(meta[,grep('Freq', colnames(meta))])
)

colnames(train_data$geneExp) = train_data$meta$Meta.specimen_id
colnames(train_data$ab) = train_data$meta$Meta.specimen_id
colnames(train_data$cytokine) = train_data$meta$Meta.specimen_id
colnames(train_data$cell_freq) = train_data$meta$Meta.specimen_id

train_data$ab[1:5, 1:13]

lapply(train_data, dim)

save(train_data, file = paste('../data/processed_datasets/training_dataset/train_' , data_to_process, '.RData', sep=''))


