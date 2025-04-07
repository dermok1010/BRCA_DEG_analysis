

# Load libraries
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)

# Load data
metadata_df_filtered <- read.csv("/home/dermot.kelly/Dermot_analysis/BRCA_DEG/data/filtered_metadata_BRCA.csv")
mirna_se_filtered <- readRDS("/home/dermot.kelly/Dermot_analysis/BRCA_DEG/data/filtered_mirna_se.rds")
load("/home/dermot.kelly/Dermot_analysis/BRCA_DEG//data/DE_results_BRCA_miRNA.RData")



# Keep only data that is not NA for pathologic stage
keep <- !is.na(metadata_df_filtered[["pathologic_stage"]])

# Subset data to reflect this
metadata_df_PS <- metadata_df_filtered[keep, ]
mirna_se_PS <- mirna_se_filtered[, keep]

dim(metadata_df_PS)
dim(assay(mirna_se_PS))  # Should have same number of columns as metadata_df_clean rows



# Running DE analysis with no covariates


dds <- DESeqDataSetFromMatrix(
  countData = assay(mirna_se_PS),
  colData = metadata_df_PS,
  design = ~ pathologic_stage
)

dds <- DESeq(dds, test = "LRT", reduced = ~1)

res <- results(dds)
head(res)

res_df <- as.data.frame(res)

# number of significant miRNAs
sign_res <- res_df %>%
  filter(padj <= 0.05)



# Now lets look at covarites
# Inspect potential covarites


summary(metadata_df_filtered$years_to_birth)
table(metadata_df_filtered$vital_status)
table(metadata_df_filtered$pathology_T_stage)
table(metadata_df_filtered$pathology_N_stage)
table(metadata_df_filtered$pathology_M_stage)
table(metadata_df_filtered$date_of_initial_pathologic_diagnosis)
table(metadata_df_filtered$radiation_therapy)
table(metadata_df_filtered$histological_type)
table(metadata_df_filtered$number_of_lymph_nodes)
table(metadata_df_filtered$race)
table(metadata_df_filtered$patient.breast_carcinoma_estrogen_receptor_status)


# 1. Convert covariates to appropriate types
metadata_df_PS$race <- as.factor(metadata_df_PS$race)
metadata_df_PS$radiation_therapy <- as.factor(metadata_df_PS$radiation_therapy)
metadata_df_PS$years_to_birth <- as.numeric(metadata_df_PS$years_to_birth)

metadata_df_PS$ER_status_created <- metadata_df_PS$patient.breast_carcinoma_estrogen_receptor_status
metadata_df_PS$ER_status_created <- droplevels(as.factor(metadata_df_PS$ER_status_created))

# Drop the single indeterminent sample
metadata_df_PS <- metadata_df_PS[metadata_df_PS$ER_status_created %in% c("positive", "negative"), ]


# 2. Define variables to keep and filter out rows with missing values
vars <- c("race", "radiation_therapy", "years_to_birth", "pathologic_stage", "ER_status_created")
keep_rows <- complete.cases(metadata_df_PS[, vars])

# 3. Filter metadata to retain complete cases
metadata_covariates <- metadata_df_PS[keep_rows, ]

# 4. Identify common patient IDs between metadata and expression data
filtered_ids <- metadata_covariates$patientID
common_ids <- intersect(filtered_ids, colnames(mirna_se_PS))

# 5. Subset expression data and metadata to retain only common patients
mirna_cov <- mirna_se_PS[, common_ids]
metadata_covariates <- metadata_covariates[metadata_covariates$patientID %in% common_ids, ]

# 6. Reorder metadata to match miRNA sample order
metadata_covariates <- metadata_covariates[match(common_ids, metadata_covariates$patientID), ]

# 7. Sanity check
stopifnot(all(colnames(mirna_cov) == metadata_covariates$patientID))

# 8. Check dimensions
dim(mirna_cov)
dim(metadata_covariates)


# Create DESeq2 dataset
dds_cov <- DESeqDataSetFromMatrix(
  countData = assay(mirna_cov),
  colData = metadata_covariates,
  design = ~ race + radiation_therapy + years_to_birth + ER_status_created + pathologic_stage
)


dds_cov <- DESeq(dds_cov)
dds_cov <- DESeq(dds_cov, test = "LRT", reduced = ~ race + radiation_therapy + ER_status_created 
                 + years_to_birth)

res_cov <- results(dds_cov)

res_cov_df <- as.data.frame(res_cov)

sign_res_cov <- res_cov_df %>%
  filter(padj <= 0.05)



# Number of DE genes
sum(res_df$padj < 0.05, na.rm = TRUE)        # original
sum(res_cov_df$padj < 0.05, na.rm = TRUE)    # covariate-adjusted





EnhancedVolcano(res_cov_df,
                lab = rownames(res_cov_df),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "DE miRNAs: Covariate-Adjusted",
                subtitle = "Pathologic Stage (Adjusted for Age, Race, Radiation)",
                caption = "Log2FC vs Adjusted p-value")




######################

# lets attempt variance partition

######################


if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}
library(variancePartition)
library(BiocParallel)
library(edgeR)

dge <- DGEList(counts = assay(mirna_cov))
dge <- calcNormFactors(dge)

# Model
form <- ~ race + radiation_therapy + years_to_birth + ER_status_created + pathologic_stage

# Voom transform
v <- voom(dge, plot = FALSE)

metadata_df <- as.data.frame(colData(dds_cov))

# Fit variance partition model
varPart <- fitExtractVarPartModel(v, form, metadata_df)


plotVarPart(varPart)

View(varPart)

# Based on above results I am going to run DE using ER

dds_er <- dds_cov
dds_er$ER_status_created <- relevel(factor(dds_er$ER_status_created), ref = "negative")

# Redefine the deisgn to focus on ER status
design(dds_er) <- ~ race + radiation_therapy + years_to_birth + pathologic_stage + ER_status_created


ddds_er <- DESeq(dds_er)
res_er <- results(dds_er, contrast = c("ER_status_created", "positive", "negative"))

er_res_df <- as.data.frame(res_er)

View(er_res_df)


sign_res_er <- er_res_df %>%
  filter(padj <= 0.05)

nrow(sign_res_er)
print(head(sign_res_cov))
print(head(sign_res_er))



# Basic volcano plot
EnhancedVolcano(res_er,
                lab = rownames(res_er),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "DE miRNAs: ER Status",
                subtitle = "Positive vs Negative (Covariate-Adjusted)",
                caption = "Log2FC vs Adjusted p-value")


