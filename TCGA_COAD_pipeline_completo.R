# =============================================================================
# PIPELINE COMPLETO — ANÁLISE TCGA-COAD + CIBERSORTx
# Caquexia e TILs no Câncer Colorretal
# =============================================================================
# Requisitos: R >= 4.3, RStudio
# Mirror CRAN recomendado para Brasil: https://cran-r.c3sl.ufpr.br/
# =============================================================================


# =============================================================================
# PRÉ-REQUISITOS — Instalação de pacotes
# =============================================================================

# BiocManager (gateway para pacotes Bioconductor)
install.packages("BiocManager")

# Pacotes Bioconductor
BiocManager::install("TCGAbiolinks")
BiocManager::install("GSVA")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")

# Pacote SummarizedExperiment
install.packages("SummarizedExperiment")

# Pacotes CRAN para análise e visualização
install.packages(c(
  "dplyr",
  "data.table",
  "ggplot2",
  "ggpubr",
  "tidyr",
  "survival",
  "survminer"
))


# =============================================================================
# PHASE 1 — DOWNLOAD E PREPARAÇÃO DOS DADOS TCGA-COAD
# =============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# Fixar select/filter para sempre usar versões do dplyr (evita conflito com AnnotationDbi)
select <- dplyr::select
filter <- dplyr::filter

# -----------------------------------------------------------------------------
# STEP 1 — Definir query
# -----------------------------------------------------------------------------

query <- GDCquery(
  project       = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Primary Tumor"
)

# Inspecionar o que será baixado (~460 amostras esperadas)
getResults(query) %>% head()


# -----------------------------------------------------------------------------
# STEP 2 — Download dos dados
# -----------------------------------------------------------------------------

# Criar pasta do projeto e definir diretório de trabalho
dir.create("~/TCGA_CRC_project", recursive = TRUE)
setwd("~/TCGA_CRC_project")
dir.create("GDCdata", showWarnings = FALSE)

# Em caso de erro de SSL, rodar antes:
# options(timeout = 600)

# Download (~3-6 GB; 20-60 min dependendo da conexão)
# Pode ser interrompido e retomado — não re-baixa arquivos já completos
GDCdownload(
  query     = query,
  method    = "api",
  directory = "GDCdata"
)


# -----------------------------------------------------------------------------
# STEP 3 — Carregar dados em objeto SummarizedExperiment
# -----------------------------------------------------------------------------

# Para sessões longas, adicionar save = TRUE para salvar .rda e evitar reprocessamento
coad_data <- GDCprepare(
  query     = query,
  directory = "GDCdata"
)

# Verificar dimensões (~60.000 genes × ~460 amostras)
dim(coad_data)


# -----------------------------------------------------------------------------
# STEP 4 — Extrair matrizes de expressão
# -----------------------------------------------------------------------------

# Verificar assays disponíveis
assayNames(coad_data)

# Extrair TPM (necessário para CIBERSORTx — NÃO usar contagens brutas)
expr_tpm <- assay(coad_data, "tpm_unstrand")

# Extrair contagens brutas (útil para análise de expressão diferencial)
expr_counts <- assay(coad_data, "unstranded")

# Verificar dimensões e preview
dim(expr_tpm)        # ~60.000 genes × ~460 amostras
expr_tpm[1:5, 1:3]


# -----------------------------------------------------------------------------
# STEP 5 — Converter IDs Ensembl para símbolos gênicos
# -----------------------------------------------------------------------------

library(org.Hs.eg.db)
library(AnnotationDbi)

# Remover sufixos de versão dos IDs Ensembl (ex: ".14")
clean_ids <- gsub("\\..*", "", rownames(expr_tpm))

# Mapear Ensembl → símbolo gênico
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = clean_ids,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"   # se múltiplos símbolos, pegar o primeiro
)

# Verificar quantos mapearam com sucesso (~25.000–30.000 de ~60.000)
sum(!is.na(gene_symbols))

# Atribuir símbolos como rownames
rownames(expr_tpm) <- gene_symbols

# Remover linhas sem símbolo (NA) e duplicatas
expr_tpm <- expr_tpm[!is.na(rownames(expr_tpm)), ]
expr_tpm <- expr_tpm[!duplicated(rownames(expr_tpm)), ]

dim(expr_tpm)   # ~25.000–28.000 genes × ~460 amostras


# -----------------------------------------------------------------------------
# STEP 6 — Salvar arquivo de mistura para CIBERSORTx
# -----------------------------------------------------------------------------

# CIBERSORTx requer:
# - Arquivo .txt delimitado por tabulação
# - Primeira coluna nomeada exatamente "GeneSymbol"
# - Uma amostra por coluna restante
# - Valores em TPM (NÃO log-transformado)

mixture_file <- cbind(GeneSymbol = rownames(expr_tpm), as.data.frame(expr_tpm))

write.table(
  mixture_file,
  file      = "CIBERSORTx_mixture_COAD.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# Verificar arquivo gerado
readLines("CIBERSORTx_mixture_COAD.txt", n = 3)
# Primeira linha deve começar: GeneSymbol  TCGA-AA-3489...


# -----------------------------------------------------------------------------
# STEP 7 — Download de dados clínicos
# -----------------------------------------------------------------------------

# Dados clínicos via GDC API
clinical <- GDCquery_clinic(
  project = "TCGA-COAD",
  type    = "clinical"
)

# Verificar colunas disponíveis
colnames(clinical)

# Selecionar colunas relevantes (usando dplyr::select explicitamente para evitar conflito)
clinical_subset <- clinical %>%
  dplyr::select(
    submitter_id,
    age_at_index,
    gender,
    race,
    ethnicity,
    vital_status,
    days_to_last_follow_up,
    days_to_death,
    ajcc_pathologic_stage,
    ajcc_pathologic_t,
    ajcc_pathologic_n,
    ajcc_pathologic_m,
    prior_malignancy,
    synchronous_malignancy,
    alcohol_history,
    icd_10_code,
    treatments_pharmaceutical_therapeutic_agents,
    treatments_pharmaceutical_treatment_or_therapy
  )

write.csv(clinical_subset, "TCGA_COAD_clinical_gdc.csv", row.names = FALSE)

# NOTA: BMI não está disponível no GDC API para TCGA-COAD.
# Dados clínicos complementares (com sobrevida estruturada) foram obtidos via cBioPortal:
# cbioportal.org → "Colon Adenocarcinoma (TCGA, GDC)" → Download → Clinical Data

# Carregar arquivo cBioPortal (após download manual)
file_path <- file.choose()   # abre janela de seleção de arquivo
cbio <- read.table(file_path,
                   sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)

colnames(cbio)

# Selecionar colunas por número (evitar problemas com coluna mal formada)
cbio_clean <- cbio[, c(2, 4, 12, 13, 22, 23, 26, 27, 28, 29, 31, 32, 37, 41, 42, 44)]

colnames(cbio_clean) <- c(
  "submitter_id",
  "age",
  "disease_free_months",
  "disease_free_status",
  "overall_survival_months",
  "overall_survival_status",
  "stage_m",
  "stage_n",
  "stage",
  "stage_t",
  "tumor_site",
  "prior_malignancy",
  "race",
  "sex",
  "tmb",
  "year_of_diagnosis"
)

head(cbio_clean)
nrow(cbio_clean)

write.csv(cbio_clean, "TCGA_COAD_clinical.csv", row.names = FALSE)
cat("Total de pacientes:", nrow(cbio_clean), "\n")

# =============================================================================
# FIM DA PHASE 1
# Arquivos gerados:
#   CIBERSORTx_mixture_COAD.txt  — upload para cibersortx.stanford.edu
#   TCGA_COAD_clinical.csv       — dados clínicos com estágio e sobrevida
# =============================================================================


# =============================================================================
# PHASE 2 — ATRIBUIÇÃO DOS GRUPOS PROXY DE CAQUEXIA
# =============================================================================

library(dplyr)
library(GSVA)

# -----------------------------------------------------------------------------
# STEP 2.1 — Carregar matriz de expressão completa
# -----------------------------------------------------------------------------

# Carregar a matriz completa com gene symbols como rownames
expr_tpm <- read.table("CIBERSORTx_mixture_COAD.txt",
                        sep = "\t", header = TRUE,
                        row.names = 1)   # row.names = 1 define GeneSymbol como rownames

dim(expr_tpm)  # ~36.000 genes × 481 amostras


# -----------------------------------------------------------------------------
# STEP 2.2 — Assinatura de caquexia e scoring com ssGSEA
# -----------------------------------------------------------------------------

# Assinatura transcriptômica de caquexia / atrofia muscular
cachexia_signature <- list(
  cachexia = c(
    "TRIM63",    # MuRF-1 — ubiquitina ligase muscular
    "FBXO32",    # Atrogin-1 — ubiquitina ligase muscular
    "UBB",       # Ubiquitina B
    "ATF4",      # Estresse de ER / atrofia
    "FOXO1",     # Regulador mestre de atrofia
    "FOXO3",     # Fator de transcrição de atrofia
    "MSTN",      # Miostatina
    "IL6",       # Inflamação sistêmica
    "TNF",       # Mediador de caquexia
    "CXCL8",     # IL-8 / cachexina do grupo
    "GDF15",     # Biomarcador de caquexia em câncer
    "PPARGC1A"   # Downregulado na caquexia
  )
)

# Converter matriz para formato numérico (necessário para GSVA)
expr_matrix <- as.matrix(expr_tpm)
expr_log2   <- log2(expr_matrix + 1)   # Log2-transformação para ssGSEA

# Verificações
dim(expr_log2)
class(expr_log2[1, 1])       # deve ser "numeric"
expr_log2[1:3, 1:3]          # deve mostrar números sem NA
rownames(expr_log2)[1:5]     # deve mostrar símbolos gênicos

# Rodar ssGSEA (~5-10 min)
gsva_params     <- ssgseaParam(expr_log2, cachexia_signature)
cachexia_scores <- gsva(gsva_params, verbose = TRUE)

# Converter para dataframe e padronizar IDs
scores_df <- as.data.frame(t(cachexia_scores))
scores_df$submitter_id <- substr(rownames(scores_df), 1, 12)

head(scores_df)


# -----------------------------------------------------------------------------
# STEP 2.3 — Merge com dados clínicos e atribuição de grupos
# -----------------------------------------------------------------------------

# Carregar dados clínicos
clinical <- read.csv("TCGA_COAD_clinical.csv")

# Padronizar formato de IDs (pontos → hífens)
scores_df$submitter_id <- gsub("\\.", "-", scores_df$submitter_id)

# Verificar correspondência
head(clinical$submitter_id)   # formato: TCGA-4N-A93T
head(scores_df$submitter_id)  # deve ser igual após gsub

# Merge
clinical <- merge(clinical, scores_df, by = "submitter_id")
cat("Amostras após merge:", nrow(clinical), "\n")

# Divisão pela mediana — grupos HIGH vs LOW caquexia
median_score <- median(clinical$cachexia, na.rm = TRUE)
cat("Score mediano de caquexia:", median_score, "\n")

clinical <- clinical %>%
  mutate(
    # Proxy transcriptômico (divisão pela mediana)
    cachexia_proxy = ifelse(cachexia >= median_score,
                            "HIGH_cachexia",
                            "LOW_cachexia"),
    # Proxy por estágio clínico (M1 = metastático = maior risco de caquexia)
    stage_proxy = case_when(
      stage_m %in% c("M1", "M1a", "M1b")                                 ~ "HIGH_cachexia",
      stage_m == "M0" & stage %in% c("Stage I", "Stage II",
                                      "Stage IIA", "Stage IIB",
                                      "Stage IIC")                        ~ "LOW_cachexia",
      TRUE ~ NA_character_
    )
  )

# Verificar tamanho dos grupos
cat("\n--- Grupos proxy transcriptômico ---\n")
table(clinical$cachexia_proxy)

cat("\n--- Grupos proxy por estágio ---\n")
table(clinical$stage_proxy, useNA = "always")

# Salvar arquivo clínico anotado
write.csv(clinical, "TCGA_COAD_clinical_annotated.csv", row.names = FALSE)
cat("\nPhase 2 concluída. Arquivo salvo.\n")

# =============================================================================
# FIM DA PHASE 2
# Arquivo gerado: TCGA_COAD_clinical_annotated.csv
#
# PHASE 3 — Upload e execução no CIBERSORTx (manual, via browser)
# Site: cibersortx.stanford.edu
# Configurações utilizadas:
#   Signature matrix:            LM22
#   Mixture file:                CIBERSORTx_mixture_COAD.txt (File Type = Mixture)
#   Batch correction:            Habilitado — B-mode
#   Permutations:                1000
#   Disable quantile norm.:      Marcado
#   Run in absolute mode:        Desmarcado
# =============================================================================


# =============================================================================
# PHASE 4 — ANÁLISE ESTATÍSTICA E FIGURAS
# (Executar após download dos resultados do CIBERSORTx)
# =============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(survival)
library(survminer)

# -----------------------------------------------------------------------------
# STEP 4.1 — Carregar e fazer merge dos resultados CIBERSORTx
# -----------------------------------------------------------------------------

# Baixar resultado do site como CSV e salvar na pasta do projeto
cibersort <- read.csv("CIBERSORTx_Results.csv",
                       stringsAsFactors = FALSE)

# Inspecionar
dim(cibersort)
colnames(cibersort)

# Filtro de qualidade — remover amostras que falharam na deconvolução
cibersort_filtered <- cibersort[cibersort$P.value <= 0.05, ]
cat("Amostras aprovadas no QC:", nrow(cibersort_filtered),
    "de", nrow(cibersort), "\n")

# Padronizar IDs (pontos → hífens, manter 12 caracteres)
cibersort_filtered$submitter_id <- gsub("\\.", "-",
                                   substr(cibersort_filtered$Mixture, 1, 12))

# Carregar arquivo clínico anotado
clinical <- read.csv("TCGA_COAD_clinical_annotated.csv")

# Merge final
merged <- merge(cibersort_filtered, clinical, by = "submitter_id")
cat("Dataset final:", nrow(merged), "amostras\n")
table(merged$cachexia_proxy)


# -----------------------------------------------------------------------------
# STEP 4.2 — Comparação estatística entre grupos
# -----------------------------------------------------------------------------

# Tipos celulares de interesse (LM22)
til_cols <- c(
  "T.cells.CD8",
  "T.cells.CD4.naive",
  "T.cells.CD4.memory.activated",
  "T.cells.regulatory..Tregs.",
  "NK.cells.activated",
  "Macrophages.M1",
  "Macrophages.M2",
  "T.cells.follicular.helper"
)

# Teste de Wilcoxon para cada tipo celular
stats_results <- lapply(til_cols, function(col) {
  high <- merged[merged$cachexia_proxy == "HIGH_cachexia", col]
  low  <- merged[merged$cachexia_proxy == "LOW_cachexia",  col]
  test <- wilcox.test(high, low)
  data.frame(
    cell_type   = col,
    p_value     = test$p.value,
    median_high = median(high, na.rm = TRUE),
    median_low  = median(low,  na.rm = TRUE),
    direction   = ifelse(median(high, na.rm = TRUE) > median(low, na.rm = TRUE),
                         "Higher in cachexia", "Lower in cachexia")
  )
})

stats_df       <- do.call(rbind, stats_results)
stats_df$p_adj <- p.adjust(stats_df$p_value, method = "BH")
stats_df$significant <- ifelse(stats_df$p_adj < 0.05, "Yes", "No")

print(stats_df)
write.csv(stats_df, "statistics_group_comparison.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# FIGURA 1 — CD8+ T cells: violin plot HIGH vs LOW caquexia
# -----------------------------------------------------------------------------

fig1 <- ggplot(
    merged[!is.na(merged$cachexia_proxy), ],
    aes(x = cachexia_proxy, y = T.cells.CD8, fill = cachexia_proxy)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  geom_jitter(width = 0.08, alpha = 0.3, size = 0.8) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c(
    "HIGH_cachexia" = "#C00000",
    "LOW_cachexia"  = "#1F4E79")) +
  scale_x_discrete(labels = c(
    "HIGH_cachexia" = "HIGH Cachexia\n(n=274)",
    "LOW_cachexia"  = "LOW Cachexia\n(n=266)")) +
  labs(
    title    = "CD8+ T Cell Infiltration in Colorectal Cancer",
    subtitle = "TCGA-COAD — CIBERSORTx LM22 Deconvolution (n=540)",
    x        = "",
    y        = "Estimated CD8+ T Cell Fraction",
    caption  = "Mann-Whitney U test; BH-adjusted p=0.0025") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave("Fig1_CD8_violin.pdf", fig1, width = 5, height = 6, dpi = 300)
cat("Fig1 salva!\n")


# -----------------------------------------------------------------------------
# FIGURA 2 — Razão CD8:Treg HIGH vs LOW caquexia
# -----------------------------------------------------------------------------

merged <- merged %>%
  mutate(
    CD8_Treg_ratio = T.cells.CD8 / (T.cells.regulatory..Tregs. + 0.001)
  )

fig2 <- ggplot(
    merged[!is.na(merged$cachexia_proxy), ],
    aes(x = cachexia_proxy, y = CD8_Treg_ratio, fill = cachexia_proxy)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  geom_jitter(width = 0.08, alpha = 0.3, size = 0.8) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c(
    "HIGH_cachexia" = "#C00000",
    "LOW_cachexia"  = "#1F4E79")) +
  scale_x_discrete(labels = c(
    "HIGH_cachexia" = "HIGH Cachexia\n(n=274)",
    "LOW_cachexia"  = "LOW Cachexia\n(n=266)")) +
  labs(
    title    = "CD8+ : Treg Ratio in Colorectal Cancer TME",
    subtitle = "TCGA-COAD — CIBERSORTx LM22 Deconvolution (n=540)",
    x        = "",
    y        = "CD8+ / Treg Fraction Ratio") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave("Fig2_CD8_Treg_ratio.pdf", fig2, width = 5, height = 6, dpi = 300)
cat("Fig2 salva!\n")


# -----------------------------------------------------------------------------
# FIGURA 3 — Composição imune do tumor: stacked bar
# -----------------------------------------------------------------------------

summary_df <- merged %>%
  filter(!is.na(cachexia_proxy)) %>%
  dplyr::select(cachexia_proxy, all_of(til_cols)) %>%
  pivot_longer(cols      = all_of(til_cols),
               names_to  = "Cell_Type",
               values_to = "Fraction") %>%
  group_by(cachexia_proxy, Cell_Type) %>%
  summarise(Mean_Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

fig3 <- ggplot(summary_df,
               aes(x = cachexia_proxy, y = Mean_Fraction, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(
    "HIGH_cachexia" = "HIGH Cachexia",
    "LOW_cachexia"  = "LOW Cachexia")) +
  labs(
    title    = "Tumor Immune Composition by Cachexia Proxy Group",
    subtitle = "TCGA-COAD — CIBERSORTx LM22 (n=540)",
    x        = "",
    y        = "Mean Immune Cell Fraction",
    fill     = "Cell Type") +
  theme_classic(base_size = 14)

ggsave("Fig3_immune_composition.pdf", fig3, width = 9, height = 6, dpi = 300)
cat("Fig3 salva!\n")


# -----------------------------------------------------------------------------
# FIGURA 4 — Curva de sobrevida: caquexia × CD8 (grupo combinado)
# -----------------------------------------------------------------------------

median_cd8 <- median(merged$T.cells.CD8, na.rm = TRUE)

merged <- merged %>%
  mutate(
    CD8_group = ifelse(T.cells.CD8 >= median_cd8, "CD8_High", "CD8_Low"),
    combined_group = case_when(
      cachexia_proxy == "HIGH_cachexia" & CD8_group == "CD8_Low"  ~ "Cachexia+/CD8-Low",
      cachexia_proxy == "HIGH_cachexia" & CD8_group == "CD8_High" ~ "Cachexia+/CD8-High",
      cachexia_proxy == "LOW_cachexia"  & CD8_group == "CD8_High" ~ "Cachexia-/CD8-High",
      cachexia_proxy == "LOW_cachexia"  & CD8_group == "CD8_Low"  ~ "Cachexia-/CD8-Low",
      TRUE ~ NA_character_
    ),
    OS_event = ifelse(overall_survival_status == "1:DECEASED", 1, 0)
  )

surv_obj <- Surv(merged$overall_survival_months, merged$OS_event)
fit      <- survfit(surv_obj ~ combined_group, data = merged)

fig4 <- ggsurvplot(fit,
  data        = merged,
  palette     = c("#C00000", "#FF9999", "#1F4E79", "#6CA0C7"),
  legend.labs = c("Cachexia+/CD8-High", "Cachexia+/CD8-Low",
                  "Cachexia-/CD8-High",  "Cachexia-/CD8-Low"),
  pval        = TRUE,
  risk.table  = TRUE,
  xlab        = "Time (months)",
  ylab        = "Overall Survival Probability",
  title       = "Survival by Cachexia Proxy + CD8 TIL Status (TCGA-COAD)",
  ggtheme     = theme_classic(base_size = 13))

ggsave("Fig4_survival.pdf", fig4$plot, width = 9, height = 7, dpi = 300)
cat("Fig4 salva!\n")

cat("\n✅ Phase 4 concluída — verifique a pasta do projeto para os PDFs das figuras\n")

# =============================================================================
# FIM DO PIPELINE
# Arquivos gerados:
#   statistics_group_comparison.csv  — tabela estatística completa
#   Fig1_CD8_violin.pdf              — CD8+ HIGH vs LOW caquexia
#   Fig2_CD8_Treg_ratio.pdf          — Razão CD8:Treg
#   Fig3_immune_composition.pdf      — Composição imune completa
#   Fig4_survival.pdf                — Curva de sobrevida combinada
# =============================================================================
