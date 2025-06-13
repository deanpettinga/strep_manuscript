# title: pgls_annotate
# author: Dean Pettinga
# date: 2024-12-26

# libraries ----

library(tidyverse) # 2.0.0
library(farrowandball)
# args ----
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("This script requires exactly 5 input arguments.", call.=FALSE)
}

# arguments
pgls_in <- args[1]
egg_in <- args[2]
fun_in <- args[3]
outfile <- args[4]

# pgls_in <- "analysis/pgls2/results/de.tsv"
# egg_in <- "analysis/emapper/out.emapper.annotations"
# fun_in <- "symlink_data/cog-2020/fun-20.tab"
# file_out <- "analysis/pgls/results/de-anno.tsv"

pgls <- read_tsv(pgls_in)

egg <- read_tsv(egg_in, skip = 4, comment = "##") %>%
  rename(query = "#query") %>%
  mutate(query = gsub("-consensus","",query)) 

cog.fun <- read_tsv(fun_in,
                    col_names = c("cog_category","HEX","cog_category_description"),
                    col_select = c("cog_category","cog_category_description"))

cog <- egg %>%
  dplyr::select(query,eggNOG_OGs,COG_category,Description) %>% 
  mutate(COG = str_extract(eggNOG_OGs, "^COG[0-9]{4}")) %>%
  dplyr::select(query,COG,Description,COG_category) %>%
  rename("orthogroup" = "query",
         "cog" = "COG",
         "cog_description" = "Description",
         "cog_category" = "COG_category")

ec <- egg %>% 
  dplyr::select(query,EC) %>%
  rename(orthogroup = query)

# cog.melt <- cog %>%
#   tidyr::separate(col = cog_category, sep = "", into = c("X","1","2","3","4","5")) %>%
#   select(-X) %>%
#   reshape2::melt(id.vars = c("orthogroup","cog","cog_description"), value.name = "cog_category", na.rm = T) %>%
#   dplyr::select(-variable) %>%
#   left_join(cog.fun, by = "cog_category")

anno <- pgls %>%
  left_join(cog, by = "orthogroup") %>%
  left_join(ec, by = "orthogroup")

write_tsv(anno, outfile)
