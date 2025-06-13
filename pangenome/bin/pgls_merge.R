# title: pgls_merge
# author: Dean Pettinga
# date: 2024-12-26

# libraries ----

library(tidyverse) # 2.0.0

# args ----
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("This script requires exactly 4 input arguments.", call.=FALSE)
}

# arguments
phene <- args[1]
in_dir <- args[2]
cutoff <- args[3]
output <- args[4]

# phene <- "de"
# in_dir <- "analysis/pgls2/de"
# cutoff <- 0.25
# output <- "analysis/pgls2/results/de.tsv"

# function ----

read_phene <- function(){

  # read each table and combine 
  res <- lapply(list.files(in_dir, full.names = T), function(x) read_tsv(x ,show_col_types = F)) %>%
    do.call(rbind, .) %>%
    mutate(test = ifelse(!is.na(lambda) & (abs(estimate) > cutoff) & (se < cutoff) & (shapiro.p > 0.05), "yes", "no")) # & (estimate_change < 0.05) 
  
  res.test <- res %>% 
    # identify OGs to test & run p.adjust we do this to limit the p.adjust penalty 
    # by including only the group we think is biologically relevant. since the adjustment
    # is higher for larger number of tests, we reduce the penalty 
    filter(test == "yes") %>%
    # p.adjust with each group and remove values for the un-tested group to avoid confusion
    mutate(FDR = p.adjust(p, method = "BH"))
  
  res.notest <- res %>% 
    filter(test == "no") %>%
    mutate(FDR = NA)
  
  return(rbind(res.test, res.notest))
}

# read results
result <- read_phene()

# save result
write_tsv(result, file = output)
# lapply(names(result), function(x) write_tsv(result[[x]], file = paste0("analysis/pgls/results/",x,".tsv")))