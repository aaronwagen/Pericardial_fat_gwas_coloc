# Run coloc on Immune expression qtls and PD GWAS
# Aaron Wagen
# Jan 2023

# Run this script on a slurm cluster with:
# ml R/4.2.0-foss-2021b
# sbatch --time=1-00:00:00 --output logs/coloc_pericard-fat_gtex8-adipose.log --open-mode=truncate --mem=128G --wrap="Rscript scripts/02_run_coloc_pericard-fat_gtex-adipose.R"

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
library(MafDb.1Kgenomes.phase3.GRCh38)
library(GenomicScores)



here::i_am("scripts/02_run_coloc_pericard-fat_gtex-adipose.R")



# Datasets

# Pericardial fat GWAS - discovery and replication

# The outcome measure of the GWAS was pericardial fat volume from cardiac MRI (ie a continueous metric), including 28,126 participants: 18,774 utilised in the discovery set and 9,387 in the replication dataset.
# The data is from UK Brain bank: it is European participants, mapped to GRCh37.
# Significant results were found in the following loci: TBX15, WARS2, and EBF2.
# Note: there are 3 GWAS: the discovery and replication GWAS that have the same columns. Allele 1 is the major/prevalent effect allele, and Allele 2 is the minor prevalence/effect allele (ie, the columns are not coded by effect and reference against the genome build). 
# The allele frequency is specific to each discovery and replicaiton dataset, and is quoted for allele 2 (this is the minor allele frequency)
# There is a 3rd GWAS of the meta-analysis. This does not have allele frequency quoted, and the nucleotides ascribed to each allele are not consistent with the discovery and replication datasets.
# In the first instance I will process only the discovery and replication datasets


# GTEx8 eQTL

## Visceral adipose tissue
## This is mapped to GRCh38 using gencode v26.
## 



# Assign arguments --------------------------------------------------------

args <- list(
    pat_discovery_gwas = c(path = here::here("raw_data",
                                             "pat_discovery_tidy_varbeta.txt"),
                           n_total = 18774,
                           cc_or_quant = "quant"),
    pat_replication_gwas = c(path = here::here("raw_data",
                                               "pat_replication_tidy_varbeta.txt"),
                             n_total = 9387,
                             cc_or_quant = "quant"),
    visceral_adipose_eqtl = here::here(""),
    raw_data_path = here::here("raw_data"),
    output_path = here::here("processed_data",
                             "coloc"),
    gtex_sample_numbers = here::here("raw_data",
                                     "GTEx_v8_sample_counts_by_tissues.csv"),
    eqtl_snp_reference_file = here::here("raw_data",
                                         "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"),
    eqtl_gtf_gencode_26_file = "~/references/GRCh38/annotation/gencode.v26.annotation.gtf.gz"
  )
  


# Load functions ----------------------------------------------------------

source(here::here("scripts",
                  "colochelpr_functions.R"))

munge_gtex8_eqtl <- function(eqtl_df){
  #' Function that munges the raw gtex8 eqtl for use with colochelpR
  #' @param eqtl_df loaded eqtl dataframe
  #' @return dataframe munged for coloc
  eqtl_df <- eqtl_df %>% 
    dplyr::mutate(eqtl_dataset = eqtl_name,
                  N = eqtl_sample_no,
                  gene = str_replace(gene_id, "\\..*$", "")) %>%  #remove suffix digit from gene_id%>% 
    dplyr::left_join(ensembl_genes_to_test,
                     by = "gene_id") %>% 
    separate(variant_id,
             into = c("CHR", "BP", "Al2", "Al1", "GRCh"), #split variant column into snp identifier columns. The first column is the ref column (Al2), the second is the alt/effect column (AL1)
             sep = "_") %>% 
    change_CHR_gencode_to_ensembl(.) %>%    # change gencode Chr format to ensembl
    generate_location_based_SNP_column() %>% 
    dplyr::select(eqtl_dataset, gene, gene_name, SNP, CHR, BP, beta = slope, se = slope_se, p.value = pval_nominal, Al1, Al2, maf, N, gene_name)  %>% 
    get_varbeta() %>% 
    change_indels() %>% 
    remove_duplicated_snps()
  
  return(eqtl_df)
}


# Create summary dataframes for GWAS and eQTL -----------------------------

# GWAS

gwas_summary <- as.data.frame(t(do.call(cbind, args[grep("_gwas", names(args))]))) %>%
  rownames_to_column(var = "gwas") %>% 
  as_tibble() %>% 
  dplyr::mutate(gwas = str_replace(gwas, "_gwas", ""),
                n_total = as.numeric(n_total))

# Create eqtl summary dataframe

# In gtex8 dataset:
## Ref allele is the allele in the reference genome - this is coloc Al2
## Alt allele is the effect allele - this is coloc Al1

eqtl_list <- list.files(args$raw_data_path) %>% 
  .[grepl("GTEx_Analysis_v8_QTLs-GTEx_Analysis_v8_eQTL_all_associations-", .)]  %>% 
  as_tibble()

eqtl_numbers <- fread(args$gtex_sample_numbers) # File of sample numbers from GTEx. The column of interest is `#RNAseq and genotyped samples`.

eqtl_numbers <- eqtl_numbers %>% 
  dplyr::select(name = Tissue,
                total_n = `# RNASeq and Genotyped samples`) %>% 
  dplyr::mutate(name = str_replace_all(name, "-", "_")) # Replace hypehns with underscore to match the file names


eqtl_summary <- eqtl_list %>% 
  dplyr::mutate(path = here::here(args$raw_data_path,
                                  value),
                name = sub("GTEx_Analysis_v8_QTLs-GTEx_Analysis_v8_eQTL_all_associations-(.*?)\\.allpairs\\.txt",
                           "\\1",
                           value),
                cc_or_quant = "quant") %>% 
  dplyr::left_join(.,
                   eqtl_numbers,
                   by = "name")

# Create output directory

dir.create(args$output_path)

# Load gtf from gencode v26 to use to lookup gene names from gene id

gtf <- rtracklayer::readGFF(args$eqtl_gtf_gencode_26_file) %>% 
  as_tibble() %>% 
  dplyr::select(gene_id, gene_name) %>% 
  dplyr::distinct(gene_name, gene_id)

 
# In this case only need to explore these 6 loci
ensembl_genes_to_test <- tibble(gene_name = c("TBX15", "WARS2", "EBF2", "RP11-418J17.1", "RP4-712E4.1", "RPS3AP12"))
ensembl_genes_to_test <- ensembl_genes_to_test %>% 
  dplyr::left_join(gtf,  # Add in gene_id
                   by = "gene_name")

# GWAS loop ---------------------------------------------------------------
    
for (gwas_number in 1:nrow(gwas_summary)) {
  # gwas_number = 1
  
  print(str_c("Loading GWAS: ", gwas_number))
  gwas_name <- gwas_summary$gwas[gwas_number]
  gwas_tidy <- read_delim(gwas_summary$path[gwas_number])
  
  
  # Find all genes within +/-1Mb of significant GWAS hits -------------------------
  
  # Usually would use this code which finds all genes within 1Mb of significant GWAS hits, and then filters the eQTL to these individual genes for munging within the loop.
  # Given there are only 3 gene-loci of interest and these are known a piori from the GWAS results, I will filter the entire eQTL for these genes and munge to start.
  
  print(str_c("GWAS: ", gwas_name))
  print("Extracting genes of interest")
  
  # Use the following code to extract all genes within +/- 1Mb of significant GWAS hits with the following code. In this case we do not require this because we are only looking in specific loci
  # ensembl_gene_ids_overlapping_1Mb_window_hit <- gwas_tidy %>%
  #   colochelpR::get_genes_within_1Mb_of_signif_SNPs(pvalue_column = "p.value",
  #                                                   CHR_column = "CHR",
  #                                                   BP_column = "BP",
  #                                                   mart = 38)
  
    
  ## Create nested results directory for gwas
  results_path_GWAS <- colochelpR::make_results_dir(results_path = args$output_path, folder_name = gwas_name)
  

 
  
  
  # eQTL loop --------------------------------------------------------------
  
  
  for(eqtl_number in 1:nrow(eqtl_summary)) {
    
    # eqtl_number <- 1
    
    print("Load eqtl")
    print(eqtl_number)
    
    eqtl_name <- eqtl_summary$name[eqtl_number]
    eqtl_path <- eqtl_summary$path[eqtl_number]
    eqtl_sample_no <- eqtl_summary$total_n[eqtl_number]
    
    print(str_c("Tissue Name is:   ", eqtl_name))
    print(str_c("Tissue Path is:   ", eqtl_path))
    print(str_c("Number of samples per tissue is: ", eqtl_sample_no))
    
    ## Load eqtl
    eqtl <- fread(eqtl_path) # Load eqtl
    
    eqtl2 <- eqtl %>% # Given only looking at 3 predefined loci, initial filter eqtl at this stage to reduce size.
    dplyr::filter(gene_id %in% ensembl_genes_to_test$gene_id) 
    
    eqtl <- munge_gtex8_eqtl(eqtl) # Use gtrex eqtl specific script to munge eqtl
    
    useful_gwas_numbers(eqtl)  
    
    # Find overlapping genes from qtl and gwas hits to save processing time in the next loop
    eqtl_genes <- unique(eqtl$gene_name)
    cis_genes_matched_gwas_eqtl <- intersect(eqtl_genes, unlist(ensembl_genes_to_test$gene_name))
    
    ## Create nested results within GWAS directory for gwas, eqtl and priors
    
    results_path_GWAS_region_eqtl <- colochelpR::make_results_dir(results_path = results_path_GWAS, folder_name = eqtl_name)
    
    ## Within results directory, create a folder for "liberal" and "robust" coloc p12 prior
    results_path_priors <- setNames(vector(mode = "list", length = 2),
                                    c("liberal", "robust"))
    
    results_path_priors$liberal <- make_results_dir(results_path = results_path_GWAS_region_eqtl, folder_name = "liberal")
    results_path_priors$robust <- make_results_dir(results_path = results_path_GWAS_region_eqtl, folder_name = "robust")
    
    
    ## Set p value cuttofs
    p12 <- setNames(c(1e-05, 5e-06),
                    c("liberal", "robust"))
    
    
    for(j in seq_along(cis_genes_matched_gwas_eqtl)){
      # j <-2 
      ensembl_gene_id_to_search <-  cis_genes_matched_gwas_eqtl[j]
      
      print(str_c(Sys.time(), " - ", gwas_name, "_gwas - ", eqtl_name, "_eqtl - ", ensembl_gene_id_to_search))
      
      # Filter eqtls for snps matching a single editing site
      eqtl_tidy_gene_filtered <- eqtl %>% 
        dplyr::filter(gene_name == ensembl_gene_id_to_search) %>% 
        colochelpR::check_coloc_data_format(.,
                                            beta_or_pval = "beta",
                                            check_maf = T) 
      
      if (nrow(eqtl_tidy_gene_filtered) == 0) {
        print(str_c("No QTLs overlapping: ", ensembl_gene_id_to_search))
        next
      }
      
      # Run coloc
      p12 <- setNames(c(1e-05, 5e-06),
                      c("liberal", "robust"))
      
      for(k in 1:length(p12)){
      # k <- 1  
        print(str_c("Results for '", names(p12[k]), "' p12 prior;  p12 = ", p12[k]))
        
        coloc_results_annotated <-
          colochelpR::get_coloc_results(df1 = gwas_tidy, df2 = eqtl_tidy_gene_filtered, 
                                        harmonise = T,  # Harmonise set to true as is will flip b-values so alleles match, allowing exploration of directionality of results
                                        df1_type = gwas_summary$cc_or_quant[gwas_number], df2_type = eqtl_summary$cc_or_quant[eqtl_number], 
                                        df1_beta_or_pval = "beta", df2_beta_or_pval = "beta",
                                        #df_1_propor_cases = args$gwas_details$prop_cases, #df_1 proportion is only useful with case control datasets
                                        df1_N = gwas_summary$n_total[gwas_number], # df_N not required for case-control datasets
                                        df2_N = eqtl_summary$total_n[eqtl_number],
                                        annotate_signif_SNP_df1_df2 = T, 
                                        key_cols = c("GWAS_1", "eqtl_dataset_2", "gene_2"), 
                                        df_1_name = "GWAS", df_2_name = "eqtl", 
                                        df1_path = gwas_summary$path[gwas_number], df2_path = eqtl_summary$path[eqtl_number],
                                        p1 = 1e-04, p2 = 1e-04, p12 = as.numeric(p12[k]))
        
        colochelpR::save_coloc_results(coloc_results_annotated, results_dir_path = results_path_priors[[names(p12[k])]])
        
        
      }
      
    }
    
#    rm(eqtl) # To save memory before loading the next eqtl # Comment out this given loading eqtl and munging at once, not in loop
  }
  
}

print("Done")

