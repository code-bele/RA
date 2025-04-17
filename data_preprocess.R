# Mapping ensmbl ids and genbank accessions to gene symbols

library(dplyr)
library(biomaRt)

setwd("/home/ibab/Sem3_NGS/Term_Paper/datasets/GSE89408")
male_early_ra = read.csv(file = "ra_early_males.csv", row.names = 1)

male_early_ra$gene = rownames(male_early_ra)
male_early_ra$gene <- sub(".*:", "", male_early_ra$gene)
male_early_ra_ensg = male_early_ra %>% filter(grepl("ENSG", gene))
rownames(male_early_ra) = male_early_ra$gene
male_early_ra_ensg$gene <- sub(".*:", "", male_early_ra_ensg$gene)
write.csv(male_early_ra_ensg, file = "male_early_ra_ensg.csv")
male_early_ra_ensg_rem = male_early_ra[!grepl("ENSG", male_early_ra$gene), ]

# Initialize Ensembl mart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query biomaRt to get gene symbols
conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                    filters = "ensembl_gene_id", 
                    values = male_early_ra_ensg$gene, 
                    mart = ensembl)

conversion <- conversion %>%
  mutate(hgnc_symbol = na_if(hgnc_symbol, "")) %>%  # Convert "" to NA
  mutate(hgnc_symbol = na_if(hgnc_symbol, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(hgnc_symbol))
write.csv(conversion, file = "ensg_conversion_file.csv")
# Merge with original dataframe
male_early_ra_ensg <- merge(male_early_ra_ensg, conversion, by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE)



male_early_ra_genbank <- male_early_ra %>% 
  filter(grepl("^[A-Z]{1,2}_?[0-9]+\\.?[0-9]*$", gene) & nchar(gene) >= 6)
male_early_ra_genbank = male_early_ra_genbank %>% filter(!grepl("PP", gene))
write.csv(male_early_ra_genbank, file = "male_early_ra_genbank.csv")
male_early_ra_ensg_genbank_rem = male_early_ra_ensg_rem %>% anti_join(male_early_ra_genbank, by = "gene")


# Processing of Genbank accesions to gene symbols was done in Python

# Process master counts matrix

master_counts = read.csv("GSE89408_GEO_count_matrix_rename.txt", sep = "\t")
master_counts$X <- sub(".*:", "", master_counts$X)

# master_counts <- as_tibble(master_counts)

library(dplyr)

male_early_ra <- male_early_ra %>%
  left_join(conversion, by = c("gene" = "ensembl_gene_id")) %>%
  mutate(Gene = ifelse(is.na(hgnc_symbol), gene, hgnc_symbol)) %>%
  dplyr::select(-hgnc_symbol)  # Explicitly calling select() from dplyr
write.csv(male_early_ra, file = "male_early_ra_new.csv")

sum(is.na(male_early_ra$Gene))  # Count missing values in the column



genbank_symbols = read.csv("genbank_with_symbol.csv")

male_early_ra <- male_early_ra %>%
  left_join(conversion, by = c("gene" = "ensembl_gene_id")) %>%
  mutate(Gene = ifelse(is.na(hgnc_symbol), gene, hgnc_symbol)) %>%
  dplyr::select(-hgnc_symbol)  # Explicitly calling select() from dplyr

male_early_ra <- male_early_ra %>%
  mutate(Gene = na_if(Gene, "")) %>%  # Convert "" to NA
  mutate(Gene = na_if(Gene, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(Gene))  # Remove rows with NA in the column


male_early_ra$gene = NULL

male_early_ra <- male_early_ra %>%
  left_join(genbank_symbols, by = c("Gene" = "gene")) %>%
  mutate(Gene = ifelse(is.na(gene_symbol), Gene, gene_symbol)) %>%
  dplyr::select(-gene_symbol)

# write.csv(master_counts, "master_counts_matrix.csv")

# Now, counts_df will have a new column "Gene"
genbank_conversion = data.frame(Gnebank_id = genbank_symbols$gene, 
                                Gene_symbol = genbank_symbols$gene_symbol)

genbank_conversion <- genbank_conversion %>%
  mutate(Gene_symbol = na_if(Gene_symbol, "")) %>%  # Convert "" to NA
  mutate(Gene_symbol = na_if(Gene_symbol, " ")) %>% # Convert " " to NA (optional)
  filter(!is.na(Gene_symbol))

write.csv(genbank_conversion, file = "genbank_conversion.csv")


