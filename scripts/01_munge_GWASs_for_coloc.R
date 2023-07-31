# Run coloc on Immune expression qtls and PD GWAS
# Aaron Wagen
# Jan 2023

# This script prepares the following GWAS's for use with coloc:



# Load libraries -------------------------------------------------------------------

library(here) # For project-specific paths
library(tidyverse) # For tidy manipulation of data
library(coloc)
library(colochelpR)
library(data.table)
library(plyranges)
library(tidytable) # For faster separation of columns of eqtl data
library(GenomicRanges) 
library(biomaRt)
library(qdapTools)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major)
library(MafDb.1Kgenomes.phase3.GRCh38) # If using GRCh38, if not then use GRCh37.
library(rutils)
library(GenomicScores)

# Note: here using the 1000 genomes project east asian population as it gave more results for MAF than other databases (eg gnomad 2. gnomad3 MafH5.gnomAD.v3.1.2.GRCh38 did not have an east asian population)
mafdb1K <- MafDb.1Kgenomes.phase3.GRCh38




# Define paths ------------------------------------------------------------


here::i_am("scripts/01_munge_GWASs_for_coloc.R")

args <- list(
  pat_discovery_gwas = c(path = here::here("raw_data",
                                                    "PAT_discovery_gwas.txt"),
                                  n_total = 18774,
                                  cc_or_quant = "quant"),
  pat_replication_gwas = c(path = here::here("raw_data",
                                                    "PAT_replication_gwas.txt"),
                                  n_total = 9387,
                                  cc_or_quant = "quant"),
  pat_meta_gwas = c(path = here::here("raw_data",
                                               "PAT_meta_gwas.txt"),
                             ntotal = 28161,
                             cc_or_quant = "quant"),
  hg19tohg38_liftover = "/nemo/lab/gandhis/home/users/wagena/references/GRCh38/tools/liftover/hg19ToHg38.over.chain",
  visceral_adipose_eqtl = here::here("txt"),
  gwas_results_path = here::here("raw_data")
  )


# Create GWAS summary dataframe from list above (choosing only top 3 args to convert to dataframe)
gwas_summary <- as.data.frame(t(do.call(cbind, args[1:3]))) %>%
  rownames_to_column(var = "gwas") %>% 
  dplyr::mutate(gwas = str_replace(gwas, "_gwas", ""))


# Define functions --------------------------------------------------------

source(here::here("scripts",
                  "colochelpr_functions.R"))


# Pericardial fat GWAS - discovery ----------------------------------------



# The outcome measure of the GWAS was pericardial fat volume from cardiac MRI (ie a continueous metric), including 28,126 participants: 18,774 utilised in the discovery set and 9,387 in the replication dataset.
# The data is from UK Brain bank: it is European participants, mapped to GRCh37.
# Significant results were found in the following loci: TBX15, WARS2, and EBF2.
# Note: there are 3 GWAS: the discovery and replication GWAS that have the same columns. Allele 1 is the major/prevalent effect allele, and Allele 2 is the minor prevalence/effect allele (ie, the columns are not coded by effect and reference against the genome build). 
# The allele frequency is specific to each discovery and replicaiton dataset, and is quoted for allele 2 (this is the minor allele frequency)
# There is a 3rd GWAS of the meta-analysis. This does not have allele frequency quoted, and the nucleotides ascribed to each allele are not consistent with the discovery and replication datasets.
# In the first instance I will process only the discovery and replication datasets

# GWAS is tidied by:
# - Renaming the column headings  and then to change them to colochelpr convention for future processing.
# - Liftover GRCh37 to 38 (to match GTEx8)
# - Removing SNPs with beta 0 and multiallelic SNPs
# - The minor allele frequency is already present
# - Adding in varbeta, the number of participants, and the gwas name, n_total and cc_or_quant columns
# - Recoding SNPs with MAF >0.5, and indels as I and D.

for (gwas_number in 1:2) {

gwas_name = gwas_summary$gwas[gwas_number] 
gwas <-  fread(gwas_summary$path[gwas_number])

print(gwas_name)
gwas <- gwas %>% 
  rename_gwas_columns() %>% 
  dplyr::rename(maf = Allele2Freq,
                p.value = P_LINREG) %>% 
  dplyr::mutate(GWAS = gwas_name,
                n_total = gwas_summary$n_total[gwas_number],
                cc_or_quant = gwas_summary$cc_or_quant[gwas_number]) %>% 
  dplyr::select(GWAS, CHR, BP, Al1, Al2, maf, beta, se, p.value, rs_number = SNP, n_total, cc_or_quant) %>%
  rutils::liftover_coord(df = .,  # Liftover from GRCh37 to 38, note this function removes the "chr" prefix, which is what we want to match our code for coloc
                         path_to_chain = args$hg19tohg38_liftover) %>% 
  relocate(., GWAS) %>% 
    generate_location_based_SNP_column(.) # If using location based SNP column for coloc, generate it using the command below from CHR and BP columns

# Generate varbeta
gwas <- colochelpR::get_varbeta(gwas) %>% 
  drop_na(varbeta) # Remove SNPs that do not have a SE/varbeta

# Change indels to I and D
gwas <- change_indels(gwas)

useful_gwas_numbers(gwas, maf_present = TRUE)

# Remove NAs and 0 from maf and beta
gwas <- remove_snps_beta_0(gwas) 

# Remove multiallelic snps by taking the snp with the lowest p value
gwas <- remove_duplicated_snps(gwas)

# Switch mafs
gwas <- correct_mafs_over05(gwas)

colochelpR::check_coloc_data_format(gwas, beta_or_pval = "beta", check_maf = T) # Check columns are the correct format

print(str_c("writing",
            args$gwas_results_path,
            "/",
            gwas_name,
            "_tidy_varbeta.txt" ))

write_delim(gwas,
            file = str_c(args$gwas_results_path,
                         "/",
                         gwas_name,
                         "_tidy_varbeta.txt"))


}



