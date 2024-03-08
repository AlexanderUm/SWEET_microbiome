#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "broom", "qiime2R", 
               "MicrobiomeStat")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")


################################################################################
# MicrobiomeStats

#-------------------------------------------------------------------------------
# Variables grid for loops 
#-------------------------------------------------------------------------------
var.grid <- expand.grid(da.vars.ls$data_set, 
                        da.vars.ls$used_ps) %>% 
                  mutate(across(everything(), as.character))


#-------------------------------------------------------------------------------
# Linda: Trend Analysis 
#-------------------------------------------------------------------------------
da.long <- list()

linda.df <- NULL  

for(i.grid in 1:nrow(var.grid)) {
  
  i.set <- var.grid[i.grid, 1]
  i.lvl <- var.grid[i.grid, 2]
  
  # Convert to mstat object
  mstat <- mStat_convert_phyloseq_to_data_obj(pss.ls[[i.set]][[i.lvl]])
  
  # Change to class 
  mstat$meta.dat <- meta.ls[[i.set]]
  
  linda.trend.res <- generate_taxa_trend_test_long(
                            data.obj = mstat,
                            subject.var = da.vars.ls$subj_var,
                            time.var = da.vars.ls$time_var,
                            group.var = da.vars.ls$group_var,
                            adj.vars = da.vars.ls$adjst_var,
                            prev.filter = da.vars.ls$min_prev,
                            feature.mt.method = "fdr",
                            feature.level = "original",
                            feature.dat.type = "count")
  
  da.long[["LinDA_trand"]][[i.set]][[i.lvl]] <- linda.trend.res$original[[1]]
  
  linda.df <- linda.trend.res$original[[1]] %>% 
                      mutate(Tax_level = i.lvl, 
                             Subset = i.set, 
                             Refference = "Sug",
                             Method ="LD_tr", 
                             Prevalence_CutOff = da.vars.ls$min_prev) %>% 
                      bind_rows(linda.df, .)
  
}

da.long.df <- list(LinDA = linda.df) 

#-------------------------------------------------------------------------------
# Write results 
#-------------------------------------------------------------------------------
save(list = c("da.long.df", "da.long"), 
     file = "out/supp/da_res3.1_long.Rdata")


