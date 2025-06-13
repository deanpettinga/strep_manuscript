# title: orthogroup_split
# author: Dean Pettinga
# date: 2024-12-18
# adapted from: https://github.com/isaisg/gfobap/blob/master/enrichment_tests/oh.split_matrix_for_pagelphyloglm.R

# libraries ----

library(tidyverse)  # 2.0.0
library(reshape2)   # 1.4.4
library(phytools)   # 2.3-0
library(nlme)       # 3.1-162

# args ----
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("This script requires exactly 4 input arguments.", call.=FALSE)
}

og_file <- args[1] # "symlink_data/Orthogroups.GeneCount.tsv"
col <- args[2]     # "symlink_data/colonization.tsv"
plant <- args[3]   # "symlink_data/microbox-plant-phenes_mean.tsv"
outdir <- args[4]  # "analysis/pgls/split/"

# import orthogroups ----

# OG <- read_tsv("../analysis/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv") %>% select(-Total)
OG <- read_tsv(og_file) %>% select(-Total)
vars <- matrixStats::rowVars(as.matrix(OG[,-1]))
OG <- cbind(OG, vars) %>%
  filter(vars > 0.25) %>%
  select(-vars)
rm(vars)

# binarize the data and exclude OGs without orthologues present in at least 1/2 of individuals
OG.bin <- ifelse(OG[2:ncol(OG)] > 0, 1,0) %>%
  as.data.frame() %>%
  cbind(OG[,1], .) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  # exclude OGs without orthologues present in at least 1/2 of individuals
  filter(sum >= ((ncol(OG)-1)/4)) %>%
  select(-sum)

# now filter OG to include the same OGs as in OG.bin
OG <- OG %>% 
  filter(Orthogroup %in% OG.bin[,1])

print("Orthogroups imported.")

# join w/ phenotypes----

col <- read_tsv(col) %>%
  dplyr::rename(mean_col_drought = mean_col.drought,
         mean_col_water = mean_col.water,
         enrichment = Enrichment,
         # Drought Enrichment score
         de = `log2(Drought/Water)` )

print("Colonization data imported.")

plant <- inner_join(# just the drought data
  read_tsv(plant) %>%
    dplyr::filter(isolate != "Mock") %>%
    reshape2::dcast(isolate + treatment ~ variable ,value.var = "value") %>%
    dplyr::filter(treatment == "drought") %>% select(-treatment) %>%
    dplyr::rename(pwc_drought = pwc,
           shoot_dryweight_g_drought = shoot_dryweight_g,
           shoot_freshweight_g_drought = shoot_freshweight_g),
  # just the water data
  read_tsv(plant) %>%
    dplyr::filter(isolate != "Mock") %>%
    reshape2::dcast(isolate + treatment ~ variable ,value.var = "value") %>%
    dplyr::filter(treatment == "water") %>% select(-treatment) %>%
    dplyr::rename(pwc_water = pwc,
                  shoot_dryweight_g_water = shoot_dryweight_g,
                  shoot_freshweight_g_water = shoot_freshweight_g),
  by = "isolate")

# join all data to one frame
phene <- inner_join(col, plant, by = "isolate") %>%
  select(-ANI_group, -biolog)

print("Plant phenotype data imported.")

# delete originals for clarity
rm(col, plant)

# split and save df ----

print("Starting split loop.")

# Create directory to contain the split files
if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}

# set up variables for df splitting in the for-loop
n_splits <- 100
span <- n_splits-1
num_paral <- 0
start <- seq(from = 1, to = nrow(OG), by = n_splits)

# Loop over the matrix splitting it
for(index in start){
        num_paral = num_paral + 1
        out_file = paste0(outdir, "/split-", num_paral, ".tsv")
        top <- span + index
        
        # for the last split including <100 orthogroups:
        if(top > nrow(OG)){
                # reasssign the last index for this split as the last OG
                top <- nrow(OG)
                
                # subset the Orthogroups df
                subset <- OG %>% column_to_rownames("Orthogroup") %>% 
                  t() %>% .[,index:top] %>% as.data.frame() %>%
                  rownames_to_column("isolate") %>%
                  # now join with the phenotype data
                  left_join(phene, ., by = "isolate")
                
                # save the file
                write_tsv(subset, file = out_file)
                
        # for all other splits of 100 orthogroups each:
        }else{
                # subset the Orthogroups df
                subset <- OG %>% column_to_rownames("Orthogroup") %>% 
                  t() %>% .[,index:top] %>% as.data.frame() %>%
                  rownames_to_column("isolate") %>%
                  # join with the phenotype data
                  left_join(phene, ., by = "isolate")
                
                # save the file
                write_tsv(subset, file = out_file)
        }
        
        print(paste0("split ", num_paral," complete."))
}

# sessionInfo() ----
sessionInfo()

