# Process each counts matrix

library(dplyr)

# Control Males:

ctrl_m = read.csv(file = "synovial_controls_male.csv")

# remove ENSG ids
ctrl_m$X <- sub(".*:", "", ctrl_m$X)
ctrl_m_ensg = ctrl_m %>% filter(grepl("ENSG", X))
ctrl_m_ensg_rem = ctrl_m[!grepl("ENSG", ctrl_m$X), ]

# remove genbank accessions
ctrl_m_genbank <- ctrl_m %>% 
  filter(grepl("^[A-Z]{1,2}_?[0-9]+\\.?[0-9]*$", X) & nchar(X) >= 6)
ctrl_m_genbank = ctrl_m_genbank %>% filter(!grepl("PP", X))
ctrl_m_ensg_genbank_rem = ctrl_m_ensg_rem %>% anti_join(ctrl_m_genbank, by = "X")

# At this point we have 3 dfs of importance:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg: this df consists of only rows that contain ensmbl ids
# 3. xx_genbank: this df consists of only rows that contain genbank ids

# load ensmbl conversion mapping df

ensg_conversion = read.csv(file = "ensg_conversion_file.csv")
ensg_conversion = ensg_conversion %>%
  mutate(hgnc_symbol = na_if(hgnc_symbol, "")) %>%  # Convert "" to NA
  mutate(hgnc_symbol = na_if(hgnc_symbol, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(hgnc_symbol))

# map the gene symbols to ensmbl ids 
ctrl_m_ensg = ctrl_m_ensg %>%
  left_join(ensg_conversion, by = c("X" = "ensembl_gene_id")) %>%
  mutate(X = ifelse(is.na(hgnc_symbol), X, hgnc_symbol)) %>%
  dplyr::select(-hgnc_symbol)  # Remove the extra column

# remove empty rows
ctrl_m_ensg = ctrl_m_ensg %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

# remove rows that were not converted
ctrl_m_ensg_converted = ctrl_m_ensg[!grepl("ENSG", ctrl_m_ensg$X), ]
ctrl_m_ensg_converted$X.y = NULL


# load genbank accesion mapping df
genbank_conversion = read.csv(file = "genbank_conversion.csv")

# map the gene symbols to genbank accesions 
ctrl_m_genbank = ctrl_m_genbank %>%
  left_join(genbank_conversion, by = c("X" = "Gnebank_id")) %>%
  mutate(X = ifelse(is.na(Gene_symbol), X, Gene_symbol)) %>%
  dplyr::select(-Gene_symbol)  # Remove the extra column

# remove empty rows
ctrl_m_genbank = ctrl_m_genbank %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

ctrl_m_genbank = ctrl_m_genbank %>% filter(!is.na(X.y))

ctrl_m_genbank$X.y = NULL

# At this point we hav 3 dfs:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg_converted: this df consists of only rows that contain ensmbl ids converted to their corresponding gene symbols
# 3. xx_genbank: this df consists of only rows that contain genbank ids converted to their corresponding gene symbols
# All 3 dfs have the same colnames so we can do a cbind to merge them together.

ctrl_m_final = bind_rows(ctrl_m_ensg_genbank_rem, ctrl_m_ensg_converted, ctrl_m_genbank)

# load list of ncRNA

nc = read.csv(file = "non-coding_RNA.txt", sep = "\t")
nc_rna = nc$symbol

ctrl_m_coding = ctrl_m_final %>% filter(!X %in% nc_rna)
ctrl_m_noncoding = ctrl_m_final %>% filter(X %in% nc_rna)

write.csv(ctrl_m_coding, file = "synovium_controls_coding_final.csv")
write.csv(ctrl_m_noncoding, file = "synovium_controls_noncoding_final.csv")







# Control Females:

ctrl_f = read.csv(file = "synovial_controls_female.csv")

# remove ENSG ids
ctrl_f$X <- sub(".*:", "", ctrl_f$X)
ctrl_f_ensg = ctrl_f %>% filter(grepl("ENSG", X))
ctrl_f_ensg_rem = ctrl_f[!grepl("ENSG", ctrl_f$X), ]

# remove genbank accessions
ctrl_f_genbank <- ctrl_f %>% 
  filter(grepl("^[A-Z]{1,2}_?[0-9]+\\.?[0-9]*$", X) & nchar(X) >= 6)
ctrl_f_genbank = ctrl_f_genbank %>% filter(!grepl("PP", X))
ctrl_f_ensg_genbank_rem = ctrl_f_ensg_rem %>% anti_join(ctrl_f_genbank, by = "X")

# At this point we have 3 dfs of importance:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg: this df consists of only rows that contain ensmbl ids
# 3. xx_genbank: this df consists of only rows that contain genbank ids

# load ensmbl conversion mapping df

ensg_conversion = read.csv(file = "ensg_conversion_file.csv")
ensg_conversion = ensg_conversion %>%
  mutate(hgnc_symbol = na_if(hgnc_symbol, "")) %>%  # Convert "" to NA
  mutate(hgnc_symbol = na_if(hgnc_symbol, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(hgnc_symbol))

# map the gene symbols to ensmbl ids 
ctrl_f_ensg = ctrl_f_ensg %>%
  left_join(ensg_conversion, by = c("X" = "ensembl_gene_id")) %>%
  mutate(X = ifelse(is.na(hgnc_symbol), X, hgnc_symbol)) %>%
  dplyr::select(-hgnc_symbol)  # Remove the extra column

# remove empty rows
ctrl_f_ensg = ctrl_f_ensg %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

# remove rows that were not converted
ctrl_f_ensg_converted = ctrl_f_ensg[!grepl("ENSG", ctrl_f_ensg$X), ]
ctrl_f_ensg_converted$X.y = NULL


# load genbank accesion mapping df
genbank_conversion = read.csv(file = "genbank_conversion.csv")

# map the gene symbols to genbank accesions 
ctrl_f_genbank = ctrl_f_genbank %>%
  left_join(genbank_conversion, by = c("X" = "Gnebank_id")) %>%
  mutate(X = ifelse(is.na(Gene_symbol), X, Gene_symbol)) %>%
  dplyr::select(-Gene_symbol)  # Remove the extra column

# remove empty rows
ctrl_f_genbank = ctrl_f_genbank %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

ctrl_f_genbank = ctrl_f_genbank %>% filter(!is.na(X.y))

ctrl_f_genbank$X.y = NULL

# At this point we hav 3 dfs:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg_converted: this df consists of only rows that contain ensmbl ids converted to their corresponding gene symbols
# 3. xx_genbank: this df consists of only rows that contain genbank ids converted to their corresponding gene symbols
# All 3 dfs have the same colnames so we can do a cbind to merge them together.

ctrl_f_final = bind_rows(ctrl_f_ensg_genbank_rem, ctrl_f_ensg_converted, ctrl_f_genbank)

# load list of ncRNA

nc = read.csv(file = "non-coding_RNA.txt", sep = "\t")
nc_rna = nc$symbol

ctrl_f_coding = ctrl_f_final %>% filter(!X %in% nc_rna)
ctrl_f_noncoding = ctrl_f_final %>% filter(X %in% nc_rna)

write.csv(ctrl_f_coding, file = "synovium_controls_f_coding_final.csv")
write.csv(ctrl_f_noncoding, file = "synovium_controls_f_noncoding_final.csv")








# Early RA Females:

ra_early_f = read.csv(file = "ra_early_females.csv")

# remove ENSG ids
ra_early_f$X <- sub(".*:", "", ra_early_f$X)
ra_early_f_ensg = ra_early_f %>% filter(grepl("ENSG", X))
ra_early_f_ensg_rem = ra_early_f[!grepl("ENSG", ra_early_f$X), ]

# remove genbank accessions
ra_early_f_genbank <- ra_early_f %>% 
  filter(grepl("^[A-Z]{1,2}_?[0-9]+\\.?[0-9]*$", X) & nchar(X) >= 6)
ra_early_f_genbank = ra_early_f_genbank %>% filter(!grepl("PP", X))
ra_early_f_ensg_genbank_rem = ra_early_f_ensg_rem %>% anti_join(ra_early_f_genbank, by = "X")

# At this point we have 3 dfs of importance:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg: this df consists of only rows that contain ensmbl ids
# 3. xx_genbank: this df consists of only rows that contain genbank ids

# load ensmbl conversion mapping df

ensg_conversion = read.csv(file = "ensg_conversion_file.csv")
ensg_conversion = ensg_conversion %>%
  mutate(hgnc_symbol = na_if(hgnc_symbol, "")) %>%  # Convert "" to NA
  mutate(hgnc_symbol = na_if(hgnc_symbol, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(hgnc_symbol))

# map the gene symbols to ensmbl ids 
ra_early_f_ensg = ra_early_f_ensg %>%
  left_join(ensg_conversion, by = c("X" = "ensembl_gene_id")) %>%
  mutate(X = ifelse(is.na(hgnc_symbol), X, hgnc_symbol)) %>%
  dplyr::select(-hgnc_symbol)  # Remove the extra column

# remove empty rows
ra_early_f_ensg = ra_early_f_ensg %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

# remove rows that were not converted
ra_early_f_ensg_converted = ra_early_f_ensg[!grepl("ENSG", ra_early_f_ensg$X), ]
ra_early_f_ensg_converted$X.y = NULL


# load genbank accesion mapping df
genbank_conversion = read.csv(file = "genbank_conversion.csv")

# map the gene symbols to genbank accesions 
ra_early_f_genbank = ra_early_f_genbank %>%
  left_join(genbank_conversion, by = c("X" = "Gnebank_id")) %>%
  mutate(X = ifelse(is.na(Gene_symbol), X, Gene_symbol)) %>%
  dplyr::select(-Gene_symbol)  # Remove the extra column

# remove empty rows
ra_early_f_genbank = ra_early_f_genbank %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

ra_early_f_genbank = ra_early_f_genbank %>% filter(!is.na(X.y))

ra_early_f_genbank$X.y = NULL

# At this point we hav 3 dfs:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg_converted: this df consists of only rows that contain ensmbl ids converted to their corresponding gene symbols
# 3. xx_genbank: this df consists of only rows that contain genbank ids converted to their corresponding gene symbols
# All 3 dfs have the same colnames so we can do a cbind to merge them together.

ra_early_f_final = bind_rows(ra_early_f_ensg_genbank_rem, ra_early_f_ensg_converted, ra_early_f_genbank)

# load list of ncRNA

nc = read.csv(file = "non-coding_RNA.txt", sep = "\t")
nc_rna = nc$symbol

ra_early_f_coding = ra_early_f_final %>% filter(!X %in% nc_rna)
ra_early_f_noncoding = ra_early_f_final %>% filter(X %in% nc_rna)

write.csv(ra_early_f_coding, file = "ra_early_f_coding_final.csv")
write.csv(ra_early_f_noncoding, file = "ra_early_f_noncoding_final.csv")








# Early RA Males:

ra_early_m = read.csv(file = "ra_early_males.csv")

# remove ENSG ids
ra_early_m$X <- sub(".*:", "", ra_early_m$X)
ra_early_m_ensg = ra_early_m %>% filter(grepl("ENSG", X))
ra_early_m_ensg_rem = ra_early_m[!grepl("ENSG", ra_early_m$X), ]

# remove genbank accessions
ra_early_m_genbank <- ra_early_m %>% 
  filter(grepl("^[A-Z]{1,2}_?[0-9]+\\.?[0-9]*$", X) & nchar(X) >= 6)
ra_early_m_genbank = ra_early_m_genbank %>% filter(!grepl("PP", X))
ra_early_m_ensg_genbank_rem = ra_early_m_ensg_rem %>% anti_join(ra_early_m_genbank, by = "X")

# At this point we have 3 dfs of importance:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg: this df consists of only rows that contain ensmbl ids
# 3. xx_genbank: this df consists of only rows that contain genbank ids

# load ensmbl conversion mapping df

ensg_conversion = read.csv(file = "ensg_conversion_file.csv")
ensg_conversion = ensg_conversion %>%
  mutate(hgnc_symbol = na_if(hgnc_symbol, "")) %>%  # Convert "" to NA
  mutate(hgnc_symbol = na_if(hgnc_symbol, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(hgnc_symbol))

# map the gene symbols to ensmbl ids 
ra_early_m_ensg = ra_early_m_ensg %>%
  left_join(ensg_conversion, by = c("X" = "ensembl_gene_id")) %>%
  mutate(X = ifelse(is.na(hgnc_symbol), X, hgnc_symbol)) %>%
  dplyr::select(-hgnc_symbol)  # Remove the extra column

# remove empty rows
ra_early_m_ensg = ra_early_m_ensg %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

# remove rows that were not converted
ra_early_m_ensg_converted = ra_early_m_ensg[!grepl("ENSG", ra_early_m_ensg$X), ]
ra_early_m_ensg_converted$X.y = NULL


# load genbank accesion mapping df
genbank_conversion = read.csv(file = "genbank_conversion.csv")

# map the gene symbols to genbank accesions 
ra_early_m_genbank = ra_early_m_genbank %>%
  left_join(genbank_conversion, by = c("X" = "Gnebank_id")) %>%
  mutate(X = ifelse(is.na(Gene_symbol), X, Gene_symbol)) %>%
  dplyr::select(-Gene_symbol)  # Remove the extra column

# remove empty rows
ra_early_m_genbank = ra_early_m_genbank %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

ra_early_m_genbank = ra_early_m_genbank %>% filter(!is.na(X.y))

ra_early_m_genbank$X.y = NULL

# At this point we hav 3 dfs:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg_converted: this df consists of only rows that contain ensmbl ids converted to their corresponding gene symbols
# 3. xx_genbank: this df consists of only rows that contain genbank ids converted to their corresponding gene symbols
# All 3 dfs have the same colnames so we can do a cbind to merge them together.

ra_early_m_final = bind_rows(ra_early_m_ensg_genbank_rem, ra_early_m_ensg_converted, ra_early_m_genbank)

# load list of ncRNA

nc = read.csv(file = "non-coding_RNA.txt", sep = "\t")
nc_rna = nc$symbol

ra_early_m_coding = ra_early_m_final %>% filter(!X %in% nc_rna)
ra_early_m_noncoding = ra_early_m_final %>% filter(X %in% nc_rna)

write.csv(ra_early_m_coding, file = "ra_early_m_coding_final.csv")
write.csv(ra_early_m_noncoding, file = "ra_early_m_noncoding_final.csv")








# Estd RA Males:

ra_estd_m = read.csv(file = "ra_estd_males.csv")

# remove ENSG ids
ra_estd_m$X <- sub(".*:", "", ra_estd_m$X)
ra_estd_m_ensg = ra_estd_m %>% filter(grepl("ENSG", X))
ra_estd_m_ensg_rem = ra_estd_m[!grepl("ENSG", ra_estd_m$X), ]

# remove genbank accessions
ra_estd_m_genbank <- ra_estd_m %>% 
  filter(grepl("^[A-Z]{1,2}_?[0-9]+\\.?[0-9]*$", X) & nchar(X) >= 6)
ra_estd_m_genbank = ra_estd_m_genbank %>% filter(!grepl("PP", X))
ra_estd_m_ensg_genbank_rem = ra_estd_m_ensg_rem %>% anti_join(ra_estd_m_genbank, by = "X")

# At this point we have 3 dfs of importance:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg: this df consists of only rows that contain ensmbl ids
# 3. xx_genbank: this df consists of only rows that contain genbank ids

# load ensmbl conversion mapping df

ensg_conversion = read.csv(file = "ensg_conversion_file.csv")
ensg_conversion = ensg_conversion %>%
  mutate(hgnc_symbol = na_if(hgnc_symbol, "")) %>%  # Convert "" to NA
  mutate(hgnc_symbol = na_if(hgnc_symbol, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(hgnc_symbol))

# map the gene symbols to ensmbl ids 
ra_estd_m_ensg = ra_estd_m_ensg %>%
  left_join(ensg_conversion, by = c("X" = "ensembl_gene_id")) %>%
  mutate(X = ifelse(is.na(hgnc_symbol), X, hgnc_symbol)) %>%
  dplyr::select(-hgnc_symbol)  # Remove the extra column

# remove empty rows
ra_estd_m_ensg = ra_estd_m_ensg %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

# remove rows that were not converted
ra_estd_m_ensg_converted = ra_estd_m_ensg[!grepl("ENSG", ra_estd_m_ensg$X), ]
ra_estd_m_ensg_converted$X.y = NULL


# load genbank accesion mapping df
genbank_conversion = read.csv(file = "genbank_conversion.csv")

# map the gene symbols to genbank accesions 
ra_estd_m_genbank = ra_estd_m_genbank %>%
  left_join(genbank_conversion, by = c("X" = "Gnebank_id")) %>%
  mutate(X = ifelse(is.na(Gene_symbol), X, Gene_symbol)) %>%
  dplyr::select(-Gene_symbol)  # Remove the extra column

# remove empty rows
ra_estd_m_genbank = ra_estd_m_genbank %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

ra_estd_m_genbank = ra_estd_m_genbank %>% filter(!is.na(X.y))

ra_estd_m_genbank$X.y = NULL

# At this point we hav 3 dfs:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg_converted: this df consists of only rows that contain ensmbl ids converted to their corresponding gene symbols
# 3. xx_genbank: this df consists of only rows that contain genbank ids converted to their corresponding gene symbols
# All 3 dfs have the same colnames so we can do a cbind to merge them together.

ra_estd_m_final = bind_rows(ra_estd_m_ensg_genbank_rem, ra_estd_m_ensg_converted, ra_estd_m_genbank)

# load list of ncRNA

nc = read.csv(file = "non-coding_RNA.txt", sep = "\t")
nc_rna = nc$symbol

ra_estd_m_coding = ra_estd_m_final %>% filter(!X %in% nc_rna)
ra_estd_m_noncoding = ra_estd_m_final %>% filter(X %in% nc_rna)

write.csv(ra_estd_m_coding, file = "ra_estd_m_coding_final.csv")
write.csv(ra_estd_m_noncoding, file = "ra_estd_m_noncoding_final.csv")







# Estd RA Females:

ra_estd_f = read.csv(file = "ra_estd_females.csv")

# remove ENSG ids
ra_estd_f$X <- sub(".*:", "", ra_estd_f$X)
ra_estd_f_ensg = ra_estd_f %>% filter(grepl("ENSG", X))
ra_estd_f_ensg_rem = ra_estd_f[!grepl("ENSG", ra_estd_f$X), ]

# remove genbank accessions
ra_estd_f_genbank <- ra_estd_f %>% 
  filter(grepl("^[A-Z]{1,2}_?[0-9]+\\.?[0-9]*$", X) & nchar(X) >= 6)
ra_estd_f_genbank = ra_estd_f_genbank %>% filter(!grepl("PP", X))
ra_estd_f_ensg_genbank_rem = ra_estd_f_ensg_rem %>% anti_join(ra_estd_f_genbank, by = "X")

# At this point we have 3 dfs of importance:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg: this df consists of only rows that contain ensmbl ids
# 3. xx_genbank: this df consists of only rows that contain genbank ids

# load ensmbl conversion mapping df

ensg_conversion = read.csv(file = "ensg_conversion_file.csv")
ensg_conversion = ensg_conversion %>%
  mutate(hgnc_symbol = na_if(hgnc_symbol, "")) %>%  # Convert "" to NA
  mutate(hgnc_symbol = na_if(hgnc_symbol, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(hgnc_symbol))

# map the gene symbols to ensmbl ids 
ra_estd_f_ensg = ra_estd_f_ensg %>%
  left_join(ensg_conversion, by = c("X" = "ensembl_gene_id")) %>%
  mutate(X = ifelse(is.na(hgnc_symbol), X, hgnc_symbol)) %>%
  dplyr::select(-hgnc_symbol)  # Remove the extra column

# remove empty rows
ra_estd_f_ensg = ra_estd_f_ensg %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

# remove rows that were not converted
ra_estd_f_ensg_converted = ra_estd_f_ensg[!grepl("ENSG", ra_estd_f_ensg$X), ]
ra_estd_f_ensg_converted$X.y = NULL


# load genbank accesion mapping df
genbank_conversion = read.csv(file = "genbank_conversion.csv")

# map the gene symbols to genbank accesions 
ra_estd_f_genbank = ra_estd_f_genbank %>%
  left_join(genbank_conversion, by = c("X" = "Gnebank_id")) %>%
  mutate(X = ifelse(is.na(Gene_symbol), X, Gene_symbol)) %>%
  dplyr::select(-Gene_symbol)  # Remove the extra column

# remove empty rows
ra_estd_f_genbank = ra_estd_f_genbank %>%
  mutate(X = na_if(X, "")) %>%  # Convert "" to NA
  mutate(X = na_if(X, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(X))

ra_estd_f_genbank = ra_estd_f_genbank %>% filter(!is.na(X.y))

ra_estd_f_genbank$X.y = NULL

# At this point we hav 3 dfs:
# 1. xx_ensg_genbank_rem: this df consists of only rows which were already gene symbols
# 2. xx_ensg_converted: this df consists of only rows that contain ensmbl ids converted to their corresponding gene symbols
# 3. xx_genbank: this df consists of only rows that contain genbank ids converted to their corresponding gene symbols
# All 3 dfs have the same colnames so we can do a cbind to merge them together.

ra_estd_f_final = bind_rows(ra_estd_f_ensg_genbank_rem, ra_estd_f_ensg_converted, ra_estd_f_genbank)

# load list of ncRNA

nc = read.csv(file = "non-coding_RNA.txt", sep = "\t")
nc_rna = nc$symbol

ra_estd_f_coding = ra_estd_f_final %>% filter(!X %in% nc_rna)
ra_estd_f_noncoding = ra_estd_f_final %>% filter(X %in% nc_rna)

write.csv(ra_estd_f_coding, file = "ra_estd_f_coding_final.csv")
write.csv(ra_estd_f_noncoding, file = "ra_estd_f_noncoding_final.csv")
