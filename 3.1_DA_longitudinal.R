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


#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------
taxa.lvls <- c("Genus", "Family")

adj.vars.ld <- var.use["Country"]


################################################################################
# MicrobiomeStats

da.long <- list()

linda.df <- NULL  

for(i.ps in names(pss.ls)) {
  
  for(i.lvl in taxa.lvls) {
    
    ps.inst <- pss.ls[[i.ps]][[i.lvl]]
    
    #---------------------------------------------------------------------------
    # Linda: Trend Analysis 
    #---------------------------------------------------------------------------
    mstat <- mStat_convert_phyloseq_to_data_obj(ps.inst)
    
    # Change to class 
    mstat$meta.dat <- mstat$meta.dat %>% 
              mutate(!!var.use["Time"] := as.numeric(.data[[var.use["Time"]]]), 
                      across(var.use[c("Country", "Subject")], 
                             as.factor), 
                     !!var.use["Group"] := factor(.data[[var.use["Group"]]], 
                                                  levels = c("Sugar", "S&SEs")))
    
    linda.trend.res <- generate_taxa_trend_test_long(
                                  data.obj = mstat,
                                  subject.var = var.use["Subject"],
                                  time.var = var.use["Time"],
                                  group.var = var.use["Group"],
                                  adj.vars = adj.vars.ld,
                                  prev.filter = prev.da.cut.off,
                                  feature.mt.method = "fdr",
                                  feature.level = "original",
                                  feature.dat.type = "count")
    
    da.long[["LinDA_trand"]][[i.ps]][[i.lvl]] <- linda.trend.res$original[[1]]
    
    linda.df <- linda.trend.res$original[[1]] %>% 
                          mutate(Tax_level = i.lvl, 
                                 Subset = i.ps, 
                                 Refference = "Sug",
                                 Method ="LD_tr", 
                                 Prevalence_CutOff = prev.da.cut.off) %>% 
                          bind_rows(linda.df, .)
 
}}


da.long.df <- list(LinDA = linda.df) 

#-------------------------------------------------------------------------------
# Write results 
#-------------------------------------------------------------------------------
save(list = c("prev.da.cut.off", "da.long.df", "da.long", "adj.vars.ld"), 
     file = "out/supp/da_res3.1_long.Rdata")


