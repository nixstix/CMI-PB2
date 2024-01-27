# REFORMATTING OF TRAINING DATA FOR INPUT INTO MOFA
# -------------------------------------------------

# This data has already been normalised by the CMI-PB team. 
# 
# We will take a quick look to see what the data looks like.

## LOAD DATA

source('../scripts/libs.R')

dat = readRDS('../data/processed_datasets/training_dataset/master_processed_training_data.RDS')
names(dat)


## SELECT DS FROM 2021 AS TRAINING SET
dat$subject_specimen$specimen_id = as.character(dat$subject_specimen$specimen_id)
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

counts = dat$pbmc_gene_expression$batchCorrected_data
counts = counts[, which(colnames(counts) %in% meta$Meta.specimen_id)]
counts = apply(counts,2, as.integer)
rownames(counts) = rownames( dat$pbmc_gene_expression$batchCorrected_data)
dim(counts)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = dat$pbmc_gene_expression$metadata[colnames(counts) , ],
                              design = ~ infancy_vac)
dds
dds = estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)
dim(assay(vsd))

vsd$subject_id = as.factor(vsd$subject_id)

bef_noHVF1 = plotPCA(vsd, "dataset") # PCA before
bef_noHVF2 = plotPCA(vsd, "subject_id") # PCA before

assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=as.vector(dds$dataset))
aft_noHVF1 = plotPCA(vsd, "dataset")
aft_noHVF2 = plotPCA(vsd, "subject_id")

# mofa is being dominated by factor 1, this may be a result of poor normalisation, therefore let's limit the geneset even further, based on HVG
topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),5000)
topVarGenes <- rownames(vsd)[topVarGenes]

dds <- dds[topVarGenes,]
dds
dds = estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)
dim(assay(vsd))

bef_plusHVF=plotPCA(vsd, "dataset") # PCA before

assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=as.vector(dds$dataset))

aft_plusHVF=plotPCA(vsd, "dataset")



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

save(train_data, file = '../data/processed_datasets/training_dataset/train_2020+21.RData')

# split into test and train

idx = unique(train_data$meta$Meta.subject_id)
idx
length(idx)
set.seed(42)
test.idx = sample(train_data$meta$Meta.subject_id, 0.40*length(idx),replace=FALSE)
test.idx = as.vector(train_data$meta[train_data$meta$Meta.subject_id %in% test.idx, 'Meta.specimen_id'])[[1]]
train.idx = as.vector(train_data$meta[!train_data$meta$Meta.specimen_id %in% test.idx, 'Meta.specimen_id'])[[1]]
                                                                                                            
train_data2 = list(
  meta = train_data$meta[train.idx, ],
  geneExp = train_data$geneExp[, train.idx], 
  ab = train_data$ab[, train.idx], 
  cytokine = train_data$cytokine[, train.idx],
  cell_freq = train_data$cell_freq[, train.idx]
)
train_data2$ab[1:5, 1:13]
lapply(train_data2, dim)

test_data = list(
  meta = train_data$meta[test.idx, ],
  geneExp = train_data$geneExp[, test.idx], 
  ab = train_data$ab[, test.idx], 
  cytokine = train_data$cytokine[, test.idx],
  cell_freq = train_data$cell_freq[, test.idx]
)
test_data$ab[1:5, 1:13]
lapply(test_data, dim)

train_data = train_data2
save(train_data, file='../data/processed_datasets/training_dataset/train_2020+21_TRAIN.RData')

train_data = test_data

save(train_data, file='../data/processed_datasets/training_dataset/train_2020+21_TEST.RData')
