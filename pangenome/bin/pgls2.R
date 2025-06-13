# title: pgls
# author: Dean Pettinga
# date: 2024-12-19
# adapted from: https://github.com/isaisg/gfobap/blob/master/enrichment_tests/oh.phyloglm.R

# libraries ----

library(tidyverse) # 2.0.0
library(reshape2) # 1.4.4
library(phytools) # 2.3-0
library(nlme) # 3.1-162
# require(caper)
# require(ape)
# require(nlme)
# require(dplyr)

# args ----
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("This script requires exactly 4 input arguments.", call.=FALSE)
}

# arguments
split_file <-  args[1] # "analysis/pgls/split/split-1.tsv"
tree_file <-   args[2] # "analysis/gtotreeN/gtotreeN.tre"
phenotype <-   args[3] # "de"
output <-      args[4] # "analysis/pgls/mean_col/de-1.tsv"

# #setwd()
# split_file <-  "analysis/pgls_split/split-1.tsv"
# tree_file <-   "analysis/gtotreeN/gtotreeN.tre"
# phenotype <-   "de"
# output <-      "analysis/pgls/de/de-1.tsv"

run_pgls <- function(DF, COL, TREE){
  
  # prep df----
  
  # print what is being tested
  print(paste("column", COL, ":", phenotype, "~", colnames(DF[,COL])))
  # current orthogroup string
  orthogroup <- colnames(DF)[COL]
  # setup the table to be input into PGLS
  dat <- DF %>% dplyr::select(all_of(c("isolate", phenotype, orthogroup)))
  colnames(dat) <- c("isolate","phenotype","OG")
  
  # choose lambda ----
  PGLS_est <- tryCatch(gls(phenotype ~ OG,
                           correlation = corPagel(value = 0.4, phy = TREE, fixed = FALSE,form = ~isolate),
                           data = dat,
                           method = "ML"),
                       error = function(e) NULL)
  
  PGLS_1 <- tryCatch(gls(phenotype ~ OG,
                         correlation = corPagel(value = 1, phy = TREE, fixed = T, form = ~isolate),
                         data = dat,
                         method = "ML"),
                     error = function(e) NULL)
  
  PGLS_0 <- tryCatch(gls(phenotype ~ OG,
                         correlation = corPagel(value = 0, phy = TREE, fixed = T, form = ~isolate),
                         data = dat,
                         method = "ML"),
                     error = function(e) NULL)
  
  print(paste0("PGLS_est lambda =", PGLS_est$modelStruct$corStruct[1]))
  print(paste0("PGLS_0 lambda =", PGLS_0$modelStruct$corStruct[1]))
  print(paste0("PGLS_1 lambda =", PGLS_1$modelStruct$corStruct[1]))
  
  # if lambda estimation with corPagel works and estimated lambda is positive
  # go ahead and compare all 3 models.
  if (!is.null(PGLS_est) && PGLS_est$modelStruct$corStruct[1] > 0){
      print("testing all three lambda...")
      # test all 3 models for lowest AIC and return model name
      best_model <- anova(PGLS_est,PGLS_0,PGLS_1) %>% 
        as.data.frame() %>% 
        arrange(AIC) %>%
        top_n(1,-AIC) %>% rownames
      print(paste0("choosing from estimate, lambda=0, or lambda=1: ",best_model))
  } else {
    print("testing lambda=0,1...")
    # test lambda=0 and lambda =1 for lowest AIC, return best
    best_model <- anova(PGLS_0,PGLS_1) %>% 
      as.data.frame() %>% 
      arrange(AIC) %>%
      top_n(1,-AIC) %>% rownames
    print(paste0("choosing from lambda=0 or lambda=1: ",best_model))
  }
  
  # save results----
  # note the model we're using
  if (best_model == "PGLS_0"){
    print("0")
    PGLS <- PGLS_0
  }
  if (best_model == "PGLS_1"){
    print("1")
    PGLS <- PGLS_1
  }
  if (best_model == "PGLS_est"){
    print("est")
    PGLS <- PGLS_est
  }
  
  # now save results
  if(is.null(PGLS)){
    print("gls() failed, saving results as NA...")
    res <- data.frame(orthogroup = orthogroup,
                      phenotype = phenotype,
                      lambda = NA,
                      estimate = NA,
                      intercept = NA,
                      se = NA,
                      t = NA,
                      p = NA,
                      shapiro.p = NA)
  } else{
    print("saving results...")
    res <- data.frame(orthogroup = orthogroup,
                      phenotype = phenotype,
                      lambda = PGLS$modelStruct$corStruct[1],
                      estimate = coef(summary(PGLS))["OG","Value"],
                      intercept = coef(summary(PGLS))["(Intercept)","Value"],
                      se = coef(summary(PGLS))["OG","Std.Error"],
                      t = coef(summary(PGLS))["OG","t-value"],
                      p = coef(summary(PGLS))["OG","p-value"],
                      shapiro.p = shapiro.test(PGLS$residuals) %>% .$p.value)
    
    # save RData
    if(!dir.exists("analysis/pgls/RData/")){dir.create("analysis/pgls/RData/")}
    save.image(file = paste0("analysis/pgls/RData/",phenotype,"-",orthogroup,".RData"))
    
    return(res)
    }
}
# IMPORT ----

# read in data from args
df <- read_tsv(split_file)
tree <- read.tree(tree_file)

# root the tree first using MRCA from outgroup Nocardioides
root_node <- phytools::findMRCA(tree=tree, tips=c("N_dokdonensis","N_humi"))
tree <- reroot(tree, node.number=root_node)
# Subset the phylogenetic tree to exclude tips without data
subtree <- drop.tip(phy = tree, tip=which(!(tree$tip.label %in% df$isolate)))
# add small value to length 0 tree branches to prevent issues in the PGLS
subtree$edge.length[subtree$edge.length == 0] <- 0.0000001

# # print what is being tested
# print(paste("column", COL, ":", phenotype, "~", colnames(DF[,COL])))
# # current orthogroup string
# orthogroup <- colnames(DF)[COL]
# # setup the table to be input into PGLS
# dat <- DF %>% dplyr::select(all_of(c("isolate", phenotype, orthogroup)))
# colnames(dat) <- c("isolate","phenotype","OG")

# RUN FUNCTION ----
table <- do.call(rbind, lapply(grep("OG",colnames(df)), function(x) run_pgls(DF = df, COL = x, TREE = subtree)))

# SAVE ----
write_tsv(table, file = output)

# sessionInfo()
sessionInfo()
