# Functions for coloc
# Author: Aaron Wagen
# July 2023

# This script contains additional helper functions for colochelper, useful in munging GWAS and coloc datasets and analysing results


# Munging GWAS ------------------------------------------------------------


rename_gwas_columns <- function(gwas_df) {
  #' Function that renames GWAS columns to all be in the same format, allowing other functions to run. Note, it will not check whether Al1 and Al2 are the correct reference
  #' @param gwas_df dataframe from a GWAS with columns to rename. Specifically will use the following labels: beta, p.value, Al1, Al2, 
  #' @return dataframe with renamed columns
  
  
  gwas_df <- gwas_df %>% 
    dplyr::rename(beta = any_of(c("BETA", "beta", "b", "Effect")),
                  p.value = any_of(c( "p.value", "p", "P", "P.VALUE", "P-value")),
                  se = any_of(c( "SE", "se", "standard.error", "StdErr")),
                  Al1 = any_of(c( "A1", "Al1", "Allele1")),
                  Al2 = any_of(c( "A2", "Al2", "Allele2"))
    )
  
  return(gwas_df)
}


switch_al1_al2 <- function(gwas_df) {
  #' Function that switched names of Al1 and Al2 columns to match the requirements of colochelpR:
  #' Al1 as the effect allele
  #' Al2 as the Alternate/Reference allele (from the genome reference).
  #' @param gwas_df dataframe from a GWAS with columns to rename. 
  #' @return dataframe with renamed columns
  gwas_df <- gwas_df %>% 
    dplyr::rename(Al1 = Al2,
                  Al2 = Al1)
  return(gwas_df)
}


change_al_columns_to_upper_case <- function(gwas_df) {
  #' Function that changes the lowercase nucleotides used in Al1 and Al2 to uppercase:
  #' @param gwas_df dataframe from a GWAS with column names Al1 and Al2. 
  #' @return dataframe with re-cased columns
  gwas_df <- gwas_df %>% 
    dplyr::mutate(Al1 = toupper(Al1),
                  Al2 = toupper(Al2))
}

generate_location_based_SNP_column <- function(gwas_df) {
  #' Function that takes a dataframe with columns for CHR and BP, and generates a location based SNP column based on the format CHR_BP
  #' @param gwas_df dataframe from a GWAS containing the columns CHR and BP
  #' @return dataframe with the combined column in the 'SNP' column
  
  gwas_df <- gwas_df %>% 
    dplyr::mutate(SNP = str_c(CHR, "_", BP))
  
  return(gwas_df)
  
}


useful_gwas_numbers <- function(gwas_df, maf_present = T) {
  #' Function that reports useful qc numbers while munging a GWAS
  #' @param gwas_df dataframe from a GWAS containing the columns CHR and BP
  #' @param maf_present, TRUE or FALSE value whether a maf column is present in the dataframe
  #' @return List of useful qc numbers
  
  
  
  print(stringr::str_c("Total number of SNPs: ",
                       nrow(gwas_df)))
  
  
  print(stringr::str_c("Number of multiallelic SNPs: ",
                       duplicated(gwas_df$SNP) %>%
                         sum()))
  
  print(stringr::str_c("Number of SNPs with beta == 0: ",
                       gwas_df %>%
                         dplyr::filter(beta == 0) %>%
                         nrow()))
  
  print(stringr::str_c("Number of SNPs with multiple nucleotides (Indels): ",
                       sum(nchar(gwas_df$Al1)>1, nchar(gwas_df$Al2)>1)))
  
  
  if (maf_present == TRUE) {
    
    print(stringr::str_c("Number of SNPs with MAF = NA: ",
                         gwas_df %>%
                           dplyr::filter(is.na(maf)) %>%
                           nrow()))
    
    print(stringr::str_c("Number of SNPs with MAF > 0.5: ",
                         gwas_df %>%
                           dplyr::filter(maf > 0.5) %>%
                           nrow()))
  }
}




remove_snps_beta_0 <- function(gwas_df) {
  #' Function that removes values where beta = 0 or NA
  #' @param gwas_df dataframe from a GWAS containing the columns beta
  #' @return dataframe with SNPs removed where beta = 0
  
  gwas_df <- gwas_df %>%
    dplyr::filter(beta != 0) 
  
  return(gwas_df)
}



remove_duplicated_snps <- function(gwas_df) {
  #' Function that removes duplicated snps by chosing the one with the lowest p value
  #' @param gwas_df dataframe from a GWAS, with the column SNP and p value
  #' @return dataframe with duplicated SNPs removed
  
  # Find duplicated SNPs
  duplicated_snps <- gwas_df %>%
    dplyr::filter(duplicated(SNP)) %>%
    dplyr::select(SNP) %>%
    unlist()
  
  # Choose duplicate with lowest p value
  filtered_duplicates <- gwas_df %>%
    dplyr::filter(SNP %in% duplicated_snps) %>%
    dplyr::group_by(SNP) %>%
    dplyr::slice_min(order_by = p.value, n=1, with_ties = F) %>% #If both have the same p value shoose at random
    dplyr::ungroup()
  
  # Add in filtered duplicates to non-duplicated snps
  gwas_tidy <- gwas_df %>%
    dplyr::filter(!SNP %in% duplicated_snps) %>% # Remove all duplicates previously identified
    dplyr::bind_rows(filtered_duplicates) %>%  # Add the rows of filtered duplicates
    dplyr::arrange(CHR, BP) 
  
  return(gwas_tidy)
  
}


change_indels <- function(gwas_df) {
  #' Function that recodes indels (SNPs with multiple nucleotides in Al1 or Al2) as Insertions and Deletions
  #' @param gwas_df dataframe from a GWAS, with the column Al1 and Al2 for Allele 1 and 2
  #' @return dataframe with duplicated SNPs removed
  
  gwas_df$Al2[nchar(gwas_df$Al1)>1] = "D"
  gwas_df$Al1[nchar(gwas_df$Al1)>1] = "I"
  gwas_df$Al1[nchar(gwas_df$Al2)>1] = "D"
  gwas_df$Al2[nchar(gwas_df$Al2)>1] = "I"
  
  return(gwas_df)
  
}


retrieve_mafs_GRCh38 <- function(gwas_df, population = c("AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF")) {
  #' Function that retrieves the minor allele frequence for a SNP based on the 1000 genomes project phase 3.
  #' @param gwas_df dataframe from a GWAS, with the columns CHR and BP from the genome build matching the one used in the maf database loaded
  #' @param population which of the [geographic populations](https://www.ensembl.org/Help/Faq?id=532) including: African (AFR_AF);  Ad Mixed American (AMR_AF); East Asian (EAS_AF); European (EUR_AF); and south asian (SAS_AF)
  #' @return dataframe with mafs, and those rows that do not have a maf removed.

  library(MafDb.1Kgenomes.phase3.GRCh38) # If using GRCh38, if not then use GRCh37.

  gwas <- gwas %>% 
    as_granges(seqnames = CHR,    # Make genomic range to get mafs.
               start = BP,
               end = BP) %>% 
    gscores(x = mafdb1K, ranges = ., pop = population) %>% #Use east asian population as the data is from Japanese donors
    as_tibble() %>% 
    drop_na(EAS_AF) %>% 
    dplyr::rename(maf = EAS_AF,
                  CHR = seqnames,
                  BP = start) %>% 
    dplyr::select(!c(end, width, strand))
}


correct_mafs_over05 <- function(gwas_df) {
  #' Function that recodes any maf values >0.5 to less than 0.5
  #' @param gwas_df dataframe from a GWAS, with the columns maf
  #' @return dataframe with recoded mafs
  #' 
gwas_df %>% 
    dplyr::mutate(
          maf =
            dplyr::case_when(
              maf > 0.5 ~ 1-maf,
              TRUE ~ maf
            )
        )
}


change_CHR_gencode_to_ensembl <- function(gwas_df) {
  #' Function that converts formating of gencode CHR to ensembl
  #' @param gwas_df dataframe  with the column CHR to transfor
  #' @return dataframe with reformated columns
 gwas_df <- gwas_df %>% 
    dplyr::mutate(CHR = str_remove(CHR, "chr")) %>% 
    dplyr::mutate(CHR = str_replace(CHR, "M", "MT")) # remove chr and recode M to MT
  return(gwas_df)
}


change_CHR_ensembl_to_gencode <- function(gwas_df) {
    #' Function that converts formating of gencode CHR to ensembl
    #' @param gwas_df dataframe  with the column CHR to transfor
    #' @return dataframe with reformated columns
  gwas_df %>% 
    dplyr::mutate(CHR = str_replace(CHR, "MT", "M")) %>%  #change MT to M
    dplyr::mutate(CHR = str_c("chr", CHR)) # Add in chr prefix
  return(gwas_df)
    
}



find_coloc_results_paths <- function(directory_to_search, pattern_to_search = ".rda") {
  #' Function that takes a directory output by colochelpR and searches within that directory and sub-directories for all files matching a pattern. It will output a dataframe with columns for: directory, file_name, gwas, eqtl, prior_type and gene.
  #' @param directory_to_search parent directory to search within
  #' @param pattern a string pattern to search in the parent directory and subdirectories
  #' @return dataframe with columns for: directory, file_name, dataset, prior_type, and gene.
  
  results_dir <- tibble(file_path =  list.files(directory_to_search, recursive = T, full.names = T, pattern = pattern_to_search)) %>% 
    dplyr::mutate(dir = file_path %>% 
                    str_replace(., "/[^/]*$", "") %>% 
                    str_replace(., "/[^/]*$", ""),
                  file_name = basename(file_path) %>% 
                    str_replace(., ".rda", ""),
                  eqtl = basename(dir),
                  gwas_path = dir %>% 
                    str_replace(., "/[^/]*$", ""),
                  gwas = basename(gwas_path),
                  prior_type = file_path %>% 
                    str_replace(., "/[^/]*$", "") %>% 
                    basename(),
                  gene = str_replace(file_name, ".+_(.*)", "\\1")
    ) 
  
  return(results_dir)
  
}









