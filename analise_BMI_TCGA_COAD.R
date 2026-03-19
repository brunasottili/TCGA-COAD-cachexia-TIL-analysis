################################################################################
# ANÁLISE DE BMI COMO CONFUNDIDOR — TCGA-COAD
# Projeto: TILs em caquexia colorrectal | Orientadora: Prof. Marilia Seelaender
# Data: Março 2026
################################################################################

# ── 1. DOWNLOAD DOS DADOS CLÍNICOS (BCR Biotab) ───────────────────────────────

query_clinical <- GDCquery(
  project       = "TCGA-COAD",
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)

GDCdownload(query_clinical, directory = "GDCdata/")

clinical_data <- GDCprepare(query_clinical, directory = "GDCdata/")

# Verificar tabelas disponíveis no objeto
print(names(clinical_data))
# Output:
# [1] "clinical_follow_up_v1.0_coad"     "clinical_radiation_coad"
# [3] "clinical_nte_coad"                "clinical_omf_v4.0_coad"
# [5] "clinical_follow_up_v1.0_nte_coad" "clinical_patient_coad"
# [7] "clinical_drug_coad"

# ── 2. EXTRAIR TABELA DE PACIENTES E VERIFICAR COLUNAS DE PESO/ALTURA ─────────

patient_tab <- clinical_data$clinical_patient_coad

# Buscar colunas relacionadas a peso, altura e BMI
bmi_cols <- grep("bmi|weight|height|mass", names(patient_tab),
                 ignore.case = TRUE, value = TRUE)
print(bmi_cols)
# Output: [1] "weight_kg_at_diagnosis" "height_cm_at_diagnosis"
# Nota: BMI não está pré-calculado — será calculado a partir de peso e altura

# ── 3. CALCULAR BMI ────────────────────────────────────────────────────────────

# Converter de character para numeric (TCGA usa "[Not Available]" para missings)
patient_tab$weight_kg_at_diagnosis <- as.numeric(patient_tab$weight_kg_at_diagnosis)
patient_tab$height_cm_at_diagnosis <- as.numeric(patient_tab$height_cm_at_diagnosis)

# Calcular BMI: peso(kg) / altura(m)²
patient_tab$bmi <- patient_tab$weight_kg_at_diagnosis /
  (patient_tab$height_cm_at_diagnosis / 100)^2

# Verificar cobertura
cat("Weight disponível:", sum(!is.na(patient_tab$weight_kg_at_diagnosis)), "de", nrow(patient_tab), "\n")
cat("Height disponível:", sum(!is.na(patient_tab$height_cm_at_diagnosis)), "de", nrow(patient_tab), "\n")
cat("BMI calculado:    ", sum(!is.na(patient_tab$bmi)), "de", nrow(patient_tab), "\n")
# Output:
# Weight disponível: 251 de 461
# Height disponível: 233 de 461
# BMI calculado:     233 de 461

# ── 4. MERGE COM GRUPOS DE CAQUEXIA ───────────────────────────────────────────

# Padronizar ID do patient_tab para formato TCGA-XX-XXXX (12 caracteres)
patient_tab$submitter_id <- substr(patient_tab$bcr_patient_barcode, 1, 12)

# Merge com scores_df (contém score ssGSEA e submitter_id)
merged_bmi <- scores_df %>%
  left_join(patient_tab %>% select(submitter_id, weight_kg_at_diagnosis,
                                    height_cm_at_diagnosis, bmi),
            by = "submitter_id")

cat("Total amostras:", nrow(merged_bmi), "\n")
cat("BMI disponível:", sum(!is.na(merged_bmi$bmi)), "\n")
# Output:
# Total amostras: 481
# BMI disponível: 255

# ── 5. TESTE DE CONFUNDIDOR: BMI ~ GRUPO CACHEXIA ─────────────────────────────

# Definir grupos HIGH/LOW pela mediana do score ssGSEA
median_score <- median(merged_bmi$cachexia)
merged_bmi$cachexia_group <- ifelse(merged_bmi$cachexia >= median_score, "HIGH", "LOW")

# Filtrar amostras com BMI disponível
df_bmi <- merged_bmi %>% filter(!is.na(bmi))

# Teste de Wilcoxon (Mann-Whitney U)
wilcox_bmi <- wilcox.test(bmi ~ cachexia_group, data = df_bmi)

cat("Mediana BMI HIGH:", round(median(df_bmi$bmi[df_bmi$cachexia_group == "HIGH"]), 1), "\n")
cat("Mediana BMI LOW: ", round(median(df_bmi$bmi[df_bmi$cachexia_group == "LOW"]), 1), "\n")
cat("p-value:", round(wilcox_bmi$p.value, 4), "\n")
# Output:
# Mediana BMI HIGH: 27.9
# Mediana BMI LOW:  26.7
# p-value: 0.1506

# ── CONCLUSÃO ─────────────────────────────────────────────────────────────────
# BMI balanceado entre grupos HIGH e LOW cachexia (p = 0.151 > 0.05)
# Cobertura: 255/481 amostras (~53%) via BCR Biotab
# O BMI não é um confundidor significativo para a análise de TILs
