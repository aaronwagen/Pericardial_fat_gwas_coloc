#
# Title: "Analyse colocalisation of pericardial fat GWAS with GTEx8 tissue specific eQTLs"
# Author: "Aaron Wagen"
# This script parses the output folders of the run_coloc_gwas_immune.R script and generates a summary table of results

## Load libraries
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
library(LDlinkR)
library(pheatmap)
library(DBI)
library(patchwork)
library(viridis)
library(sciRmdTheme)
library(paletteer)



# Set paths
here::i_am("scripts/03_process_coloc_results_pericard-fat_gtex-adipose.R")

args <- list(
  coloc_path = here::here("processed_data",
                           "coloc"),
  coloc_results_summary_file = here::here("processed_data",
                                          "coloc",
                                          "coloc_summary.Rds"),
  eqtl_gtf_gencode_26_file = "~/references/GRCh38/annotation/gencode.v26.annotation.gtf.gz"
)



# Load functions
source(here::here("scripts",
                  "colochelpr_functions"))

# Methods


# Process Coloc results


# Create table with file paths, datasets, prior type and editing site


results_dir <- find_coloc_results_paths(directory_to_search = args$coloc_path, pattern_to_search = ".rda")


# Split into liberal/robust results
dataset_dir <- setNames(results_dir %>% dplyr::group_split(prior_type),
                        c(results_dir %>% .[["prior_type"]] %>% unique() %>% sort()))
# Priors df
priors_df <- tibble(prior_type = c("liberal", "robust"),
                    p12 = c(1e-05, 5e-06))
all_results <- setNames(vector(mode = "list", length = length(dataset_dir)),
                        names(dataset_dir))


# List of datasets analysed

dataset_variable <- "eqtl" # To explore the difference between using eqtl and gwas

dataset_list <- unique(results_dir[dataset_variable]) %>% unlist()



for(i in 1:length(dataset_dir)){
  
  priors_df_filtered <- 
    priors_df %>% 
    dplyr::filter(prior_type == names(dataset_dir[i]))
  
  dataset_file_paths <- dataset_dir[[i]]
  
  dataset_names <- dataset_file_paths[dataset_variable] %>% unique() %>% unlist()
  
  results_list <- vector(mode = "list", length = length(dataset_names))
  
  for(j in 1:length(dataset_names)){
    
    dataset <- dataset_names[j]
    
    print(dataset)
    
    dir_to_load <- 
      dataset_file_paths %>% 
      dplyr::filter(!!sym(dataset_variable) == dataset_names[j]) %>% 
      .[["dir"]] %>% 
      unique()
    
    dir_to_load <- str_c(dir_to_load, "/", priors_df_filtered$prior_type)
    
    for(k in 1:length(dir_to_load)){
      
      print(dir_to_load[k])

      if(dataset %in% dataset_list){
        
        results <- 
          colochelpR::merge_coloc_summaries(dir_to_load[k], add_signif_SNP = F, recursive = T, pattern = ".rda") %>% 
          dplyr::select(GWAS_1, gene_2, everything(), -eqtl_dataset_2)
        
      }
      
      results <- 
        results %>% 
        dplyr::mutate(eqtl = dataset,
                      p12 = priors_df_filtered$p12) %>% 
        dplyr::select(GWAS_1, eqtl,  gene_2, everything()) #
      
      if(k == 1){
        
        results_list[[j]] <- results 
        
      } else{
        
        results_list[[j]] <- 
          results_list[[j]] %>% 
          dplyr::bind_rows(results)
        
      }
      
    }
    
  }
  
  all_results[[i]] <- results_list
  
}


# Given biomart often fails, use gtf from gencode v26 (the one used in gtex8) to find the gene names from ensembl ids.

gtf <- rtracklayer::readGFF(args$eqtl_gtf_gencode_26_file) %>% 
  as_tibble() %>% 
  dplyr::select(gene_id, gene_name) %>% 
  dplyr::distinct(gene_name, gene_id) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

results <- 
  all_results %>% 
  lapply(., function(x){
    
    x %>% 
      qdapTools::list_df2df() %>%
      dplyr::left_join(., gtf,
                       by = c("gene_2" = "gene_id")) %>% 
      # biomart_df(columnToFilter = "gene_2",
      #            mart = 38,
      #            mirror = "asia",
      #            attributes = c("ensembl_gene_id", "hgnc_symbol"),
      #            filter = c("ensembl_gene_id")) %>%
      dplyr::select(GWAS_1, eqtl, gene_2, gene_name, everything(), -X1)

  }) 

saveRDS(results, file = args$coloc_results_summary_file)
  



