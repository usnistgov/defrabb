library(tidyverse)

xlog <- "../../logs/exclusions/testC/subtract_GRCh38_asm17aChr21_dipcall-z2k.log" %>%
    read.csv(sep = "\t") %>%
    as_tibble()

init <- xlog[1, ]
after <- xlog[-1, ]

init_length <- init$resulting_length

sprintf("Starting length before exclusions: %ibp", init_length)

after %>%
    mutate(exclusion = basename(exclusion),
           percent_decrease = 100 * resulting_length / init_length,
           bp_diff = lag(resulting_length) - resulting_length)
          
