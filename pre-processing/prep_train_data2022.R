# REFORMATTING OF TRAINING DATA FOR INPUT INTO MOFA
# -------------------------------------------------

# This data has already been normalised and batch-corrected by the CMI-PB team. 
# 

data_to_process='2022_dataset'

## LOAD DATA

source('../scripts/libs.R')

dat = readRDS('../data/processed_datasets/prediction_dataset/master_processed_prediction_data.RDS')
names(dat)

## SELECT DS FROM 2022 AS TRAINING SET
dat$subject_specimen$specimen_id = as.character(dat$subject_specimen$specimen_id)
dim(dat$subject_specimen)
sampleIdx = dat$subject_specimen[dat$subject_specimen$dataset == '2022_dataset', 'specimen_id']


dat$subject_specimen = dat$subject_specimen[dat$subject_specimen$specimen_id %in% sampleIdx$specimen_id,]
table(dat$subject_specimen$dataset)
dim(dat$subject_specimen)


# MERGE DATA
## merge all matrices to replace empty values with NAs

meta = dat$subject_specimen
colnames(meta) = paste("Meta.", colnames(meta), sep='')

# ab
x = as.data.frame(t(dat$abtiter$processed_similar_to_training))
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# cytokines
x = as.data.frame(t(dat$plasma_cytokine_concentrations$processed_similar_to_training))
colnames(x) = paste("Cytokine.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))


# cell_freq
x = as.data.frame(t(dat$pbmc_cell_frequency$processed_similar_to_training))
colnames(x) = paste("Freq.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# expression

counts_long = dat$pbmc_gene_expression$raw_from_database
head(counts_long)

counts_long = counts_long[,c(1:3)]
head(counts_long)

counts = counts_long %>% 
  pivot_wider(names_from=specimen_id, values_from = raw_count)
counts = as.data.frame(counts)
head(counts)
dim(counts)
rn = counts$versioned_ensembl_gene_id
rownames(counts) = rn
counts[1:5, 1:5]
dim(counts)

# select only counts that have an ID in the metadata table, and genes that are present in the processed data
counts = counts[which(rownames(counts) %in% rownames(dat$pbmc_gene_expression$processed_similar_to_training)) , which(colnames(counts) %in% meta$Meta.specimen_id)]
dim(counts)
counts[1:5, 1:5]
rownames(counts)
colnames(counts)

rownames(meta) = meta$Meta.specimen_id
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta[colnames(counts), ],
                              design = ~ Meta.infancy_vac)
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


