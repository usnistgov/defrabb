################################################################################
####
#### Utility functions for working with extended csv files generate by the
#### happy benchmarking pipeline
####
################################################################################

## Reading and tidying happy benchmarking extended csv files ###################
load_happy_extended <-
  function(path,
           tidy = FALSE,
           filter = FALSE,
           subsets,
           metrics) {
    ## Reading extended csv
    happy_df <- read_happy_extended(path)

    ## Filter df based on metrics and subsets
    if (filter == TRUE) {
      happy_df <- filter_happy_extended(happy_df, subsets, metrics)
    }

    ## Tidying df
    if (tidy == TRUE) {
      happy_df <- tidy_happy_extended(happy_df)
    }

    return(happy_df)
  }


## Reading happy benchmark extended csv files ##################################
read_happy_extended <- function(path) {
  ## Defining col types to avoid ti/tiv logical issues
  ctypes <-
    stringr::str_c(c(rep("c", 7), rep("d", 9), rep("dccddcd", 7)),
      sep = "",
      collapse = ""
    )

  ## Reading as csv file
  happy_df <- readr::read_csv(path, col_types = ctypes)

  return(happy_df)
}

## Tidying happy benchmark extended csv files ##################################
filter_happy_extended <-
  function(happy_df,
           subsets = NULL,
           metrics = NULL,
           subtypes = NULL,
           min.subset.size = 0,
           filter_pass = TRUE) {
    ## Subsetting to defined stratification set
    if (!is.null(subsets)) {
      happy_df <- dplyr::filter(happy_df, Subset %in% subsets)
    }

    ## Only Keeping passing variants
    if (filter_pass == TRUE) {
      happy_df <- dplyr::filter(happy_df, Filter == "PASS")
    } else if (!is.null(metrics)) {
      ## Adding Filter to list of metrics to keep
      metrics <- c("Filter", metrics)
    }

    ## Subsetting to defined metrics set
    if (!is.null(metrics)) {
      ## Keeping main cols - might want to add Genotype, QQ.Field, and QQ values
      cols_to_keep <-
        c("Type", "Subtype", "Subset", metrics)

      happy_df <-
        dplyr::select(happy_df, all_of(cols_to_keep))
    }

    ## Subsetting to defined variant subtypes
    if (!is.null(subtypes)) {
      happy_df <- dplyr::filter(happy_df, Subtype %in% subtypes)
    }

    ## Excluding subsets with a limited number of variants
    if (min.subset.size > 0) {
      happy_df <- dplyr::filter(happy_df, Subset.Size > min.subset.size)
    }

    ## Return filtered and subsetted data frame
    return(happy_df)
  }

tidy_happy_extended <- function(happy_df, scale_phred = TRUE) {
  ## Convert to a long df
  long_happy_df <- happy_df %>%
    select(
      Type,
      Subtype,
      Subset,
      METRIC.Recall,
      METRIC.Precision,
      TRUTH.TP,
      TRUTH.FN,
      QUERY.TP,
      QUERY.FP
    ) %>%
    pivot_longer(
      cols = c("METRIC.Recall", "METRIC.Precision"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(metric = str_remove(metric, "METRIC."))

  ##
  ci_happy_df <- long_happy_df %>%
    group_by(Type, Subtype, Subset, metric) %>%
    mutate(metric_ci = list(
      get_metric_ci(metric,
        alpha = 0.05,
        TRUTH.TP, TRUTH.FN, QUERY.TP, QUERY.FP
      )
    )) %>%
    unnest(cols = c(metric_ci)) %>%
    select(-TRUTH.TP, -TRUTH.FN, -QUERY.TP, -QUERY.FP, -PointEst) %>%
    rename(
      lci = Lower,
      uci = Upper
    )

  if (scale_phred == TRUE) {
    ## Converting metrics to phred
    ci_happy_df <- mutate(
      ci_happy_df,
      value_phred = convert_to_phred(value),
      lci_phred = convert_to_phred(lci),
      uci_phred = convert_to_phred(uci)
    )
  }

  return(ci_happy_df)
}

## Calculating Benchmark Metrics Geometric Means ###############################
calc_geometric_mean <- function(value) {
  exp(mean(log(value)))
}

## Converting phred scale benchmarking metrics #################################
convert_to_phred <- function(metric) {
  -10 * log(1 - metric) / log(10)
}

## Calculate Metric Confidence Intervals #######################################

#' Calculate Metric Confidence Intervals
#'
#' Calculate confidence intervals for benchmarking performance metrics using
#' output from happy small variant benchmarking tool. Metrics described in
#' Krusche et al. 2019 (10.1038/s41587-019-0054-x). Confidence intervals
#' calculated using binconf assuming binomial probability distribution.
#' Currently calculates confidence intervals for recall and precision. All four
#' TP, TN, and FP counts are not needed for all metrics. Only TRUTH.TP and
#' TRUTH.FN needed for Recall and QUERY.TP and QUERY.FP needed for Precision.
#' Future extensions include CIs for Frac_NA, and F_Score.
#'
#' @param metric character string either RECALL or PRECISION
#' @param alpha probability of a type I error, passed to binconf
#' @param TRUTH.TP number of Truth TPs
#' @param TRUTH.FN number of Truth FNs
#' @param QUERY.TP number of Query TPs
#' @param QUERY.FP number of Query FPs.
#'
#' @return data frame with confidence intervals
#' @export
#'
#' @examples get_metric_ci("RECALL", 0.05, 999, 1)
get_metric_ci <- function(metric,
                          alpha,
                          TRUTH.TP,
                          TRUTH.FN,
                          QUERY.TP,
                          QUERY.FP) {
  ## Converting to upper to allow lower case metric names
  metric <- toupper(metric)
  ## Recall defined as TRUTH.TP /(TRUTH.TP + TRUTH.FN)
  if (metric == "RECALL") {
    x <- TRUTH.TP
    n <- TRUTH.TP + TRUTH.FN
    ## Precision defined as QUERY.TP /(QUERY.TP + QUERY.FP)
  } else if (metric == "PRECISION") {
    x <- QUERY.TP
    n <- QUERY.TP + QUERY.FP
  } else {
    stop("metric must be RECALL or PRECISION")
  }
  ## Calculating confidence intervals
  ci_df <- Hmisc::binconf(
    x = x,
    n = n,
    alpha = alpha,
    return.df = TRUE
  )
  return(ci_df)
}
