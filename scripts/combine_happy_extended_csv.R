##########################################
## This script is to summarize metrics from the hap.py results_extended.csv
## and outputs metrics for all stratitications in stratification .tsv and *
###########################################

library(tidyverse)
library(fs)
library(here)

## directory path for results_extended.csv files.  If multiple benchmarks are in the same directory script will output single file with all benchmark metrics.

data_dir <- here("giab-tibanna-runs/defrabb-HPRC-44")

############################################
## Function to read in  metrics
############################################

read_ext_results <- function(dir) {
  extended_files <- fs::dir_ls(dir, regexp = "extended.csv$", recurse = TRUE)
  ext_results <- extended_files %>%  
          set_names(.,fs::path_file(.)) %>% 
          map_dfr(read_csv, show_col_types = FALSE, .id="source")
  return(ext_results)
}

extended_results <- read_ext_results(data_dir)

############################################
## Create list of stratifcations to subset to for summary
############################################

## Including all subsets
subsets <- extended_results$Subset %>% unique()

############################################
## Filter the orginal data frame for Type (SNP/INDEL) and desired stratifications (Subset). Rename colnames to note
## metric Type in larger table once combined.
############################################
## genrally the function does the following
##    1) subset to desired stratifications and subtype = "*"
##    2) rename colnames to specify Type for the metric (e.g. SNP or INDEL)
##    3) combine metric types into single table where all Type.metric cols are present.
##    4) Add additional columns for subtypes I16 and D16 (insertions and deletions 16+ bp in length) for specific
##       metrics for the same subsets (stratifications) as prior columns
##    5) combines all subtables into single output table (4 rows x 16 cols)
##    6) 8/4 added additional metrics and new calculated metrics to explore ancestry strat stats
############################################

type_rename <- function(combined_results) {
  # create subtable for SNPs and INDELs with Subtype = *
  SNP_ext_results<- combined_results %>%
    filter(Type == "SNP" & Subset %in% subsets & Subtype == "*" ) %>%
    rename( SNP.Recall = METRIC.Recall,
            SNP.Precision = METRIC.Precision,
            SNP.Frac_NA = METRIC.Frac_NA,
            SNP.TRUTH.TOTAL.het_hom_ratio = TRUTH.TOTAL.het_hom_ratio,
            SNP.TRUTH.TP = TRUTH.TP,
            SNP.TRUTH.TP.het = TRUTH.TP.het,
            SNP.TRUTH.TP.homalt = TRUTH.TP.homalt,
            SNP.TRUTH.FN = TRUTH.FN,
            SNP.TRUTH.FN.het = TRUTH.FN.het,
            SNP.TRUTH.FN.homalt = TRUTH.FN.homalt,
            SNP.QUERY.TOTAL.het_hom_ratio = QUERY.TOTAL.het_hom_ratio,
            SNP.QUERY.FP = QUERY.FP,
            SNP.QUERY.FP.het = QUERY.FP.het,
            SNP.QUERY.FP.homalt = QUERY.FP.homalt,
            SNP.FP.gt = FP.gt,
            SNP.FP.al = FP.al,
            SNP.TRUTH.TOTAL = TRUTH.TOTAL,
            SNP.QUERY.TOTAL = QUERY.TOTAL,
            SNP.Subset.IS_CONF.Size = Subset.IS_CONF.Size,
            SNP.Subset.Size = Subset.Size) %>% 
    mutate("SNP.Truth/kb" = (as.numeric(SNP.TRUTH.TOTAL)/as.numeric(SNP.Subset.IS_CONF.Size))*1000) %>%
    mutate("SNP.Query/kb" = (as.numeric(SNP.QUERY.TOTAL)/as.numeric(SNP.Subset.Size))*1000) %>%
    select(-Type,-Subtype, -Filter, -SNP.TRUTH.TOTAL, -SNP.QUERY.TOTAL, -SNP.Subset.IS_CONF.Size, -SNP.Subset.Size)
  
  INDEL_ext_results <- combined_results %>%
    filter(Type == "INDEL" & Subset %in% subsets & Subtype == "*") %>%
    rename(INDEL.Recall = METRIC.Recall,
           INDEL.Precision = METRIC.Precision,
           INDEL.Frac_NA = METRIC.Frac_NA,
           INDEL.TRUTH.TOTAL.het_hom_ratio = TRUTH.TOTAL.het_hom_ratio,
           INDEL.TRUTH.TP = TRUTH.TP,
           INDEL.TRUTH.TP.het = TRUTH.TP.het,
           INDEL.TRUTH.TP.homalt = TRUTH.TP.homalt,
           INDEL.TRUTH.FN = TRUTH.FN,
           INDEL.TRUTH.FN.het = TRUTH.FN.het,
           INDEL.TRUTH.FN.homalt = TRUTH.FN.homalt,
           INDEL.QUERY.TOTAL.het_hom_ratio = QUERY.TOTAL.het_hom_ratio,
           INDEL.QUERY.FP = QUERY.FP,
           INDEL.QUERY.FP.het = QUERY.FP.het,
           INDEL.QUERY.FP.homalt = QUERY.FP.homalt,
           INDEL.FP.gt = FP.gt,
           INDEL.FP.al = FP.al,
           INDEL.TRUTH.TOTAL = TRUTH.TOTAL,
           INDEL.QUERY.TOTAL = QUERY.TOTAL,
           INDEL.Subset.IS_CONF.Size = Subset.IS_CONF.Size,
           INDEL.Subset.Size = Subset.Size) %>% 
    mutate("INDEL.Truth/kb" = (as.numeric(INDEL.TRUTH.TOTAL)/as.numeric(INDEL.Subset.IS_CONF.Size))*1000) %>%
    mutate("INDEL.Query/kb" = (as.numeric(INDEL.QUERY.TOTAL)/as.numeric(INDEL.Subset.Size))*1000) %>%
    select(-Type,-Subtype, -Filter, -INDEL.TRUTH.TOTAL, -INDEL.QUERY.TOTAL, -INDEL.Subset.IS_CONF.Size, -INDEL.Subset.Size)

  # join subtables for SNPs and INDELs where Subtype = *
  star_join <- left_join(SNP_ext_results,
                         INDEL_ext_results,
                         by = c("source", "Subset"))

  # create subtables for insertions and deletions 16+ bp in length
  ext_results_I16_ext_results <- combined_results %>%
    filter(Type == "INDEL"& Subset %in% subsets & Subtype == "I16_PLUS") %>%
    rename(I16_PLUS.Recall = METRIC.Recall,
           I16_PLUS.Precision = METRIC.Precision,
           I16_PLUS.TP = TRUTH.TP,
           I16_PLUS.FP = QUERY.FP,
           I16_PLUS.FN = TRUTH.FN) %>%
    select(source, Subset, I16_PLUS.Recall,I16_PLUS.Precision,I16_PLUS.TP, I16_PLUS.FP, I16_PLUS.FN)

  ext_results_D16_ext_results <- combined_results %>%
    filter(Type == "INDEL"& Subset %in% subsets & Subtype == "D16_PLUS") %>%
    rename(D16_PLUS.Recall = METRIC.Recall,
           D16_PLUS.Precision = METRIC.Precision,
           D16_PLUS.TP = TRUTH.TP,
           D16_PLUS.FP = QUERY.FP,
           D16_PLUS.FN = TRUTH.FN) %>%
    select(source, Subset, D16_PLUS.Recall, D16_PLUS.Precision, D16_PLUS.TP, D16_PLUS.FP, D16_PLUS.FN)

  # join subtables where insertions and deletions 16+ bp in length
  plus_join <- left_join(ext_results_I16_ext_results,
                         ext_results_D16_ext_results,
                         by = c("source", "Subset"))

  # output table created by joining metrics for Subtype = * and 16+bp
  ext_results_combined_renamed <-left_join(star_join,
                                           plus_join,
                                           by = c("source" , "Subset"))

return(ext_results_combined_renamed)
}

ext_results_combined_renamed <- type_rename(extended_results)
write_tsv(ext_results_combined_renamed, 
          here("scratch/defrabb-HPRC44_combined-happy-extended_2022-02-28.tsv"))


