# RA Project: Re-Analysis
# Clean and prep counts matrices

library(dplyr)


setwd("/home/ibab/Sem3_NGS/Term_Paper/datasets")

data = read.csv(file = "GSE89408_GEO_count_matrix_rename.txt", sep = "\t", row.names = 1)
data_control = data %>% select(starts_with("n"))
metadata = read.csv(file = "GSE89408_metadata.csv")
col_names = metadata$sex
control_colnames = col_names[1:28]
colnames(data_control) = control_colnames
write.csv(data_control, file = "synovial_control_counts.csv")

data_ra = data %>% select(starts_with("r"))

data_ra_early = data_ra %>% select(RA_tissue_1:RA_tissue_57)
metadata_ra_early = metadata %>% filter(disease == "Rheumatoid arthritis (early)") %>% arrange(Sample.Name)
colnames(data_ra_early) = metadata_ra_early$sex
write.csv(data_ra_early, file = "ra_early_counts.csv")

data_ra_estd = data_ra %>% select(RA_tissue_58:RA_tissue_152)
metadata_ra_estd = metadata %>% filter(disease == "Rheumatoid arthritis (established)") %>% arrange(Sample.Name)
colnames(data_ra_estd) = metadata_ra_estd$sex
write.csv(data_ra_estd, file = "ra_estd_counts.csv")

# segregate males from females for all

controls = read.csv(file = "synovial_control_counts.csv", row.names = 1)
control_male = controls %>% select(starts_with("m"))
write.csv(control_male, file = "synovial_controls_male.csv")

control_female = controls %>% select(starts_with("f"))
write.csv(control_female, file = "synovial_controls_female.csv")

ra_early = read.csv(file = "ra_early_counts.csv",row.names = 1)
ra_early_males = ra_early %>% select(starts_with("m"))
write.csv(ra_early_males, file = "ra_early_males.csv")


ra_early_females = ra_early %>% select(starts_with("f"))
write.csv(ra_early_females, file = "ra_early_females.csv")


ra_estd = read.csv(file = "ra_estd_counts.csv", row.names = 1)
ra_estd_males = ra_estd %>% select(starts_with("m"))
write.csv(ra_estd_males, file = "ra_estd_males.csv")

ra_estd_females = ra_estd %>% select(starts_with("f"))
write.csv(ra_estd_females, file = "ra_estd_females.csv")
