#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "ggsignif", "broom", "vegan", 
               "metagenomeSeq", "qiime2R", "picante", "phangorn", "FSA", 
               "knitr", "MicrobiomeStat", "ggvegan", "ggpmisc", 
               "Maaslin2", "ComplexHeatmap", "circlize", "splinectomeR")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")
load("out/supp/5_Responders.Rdata")



################################################################################
# Outline
################################################################################
resp.da.vars.ls <- list(resp_cols = c("Resp_Body_weight"), 
                     samp_gr = c("Sugar", "S&SEs"), 
                     samp_set = c("all"), 
                     used_ps = c("Genus", "Family"), 
                     min_prev = 0.5,
                     subj_var = "Subject", 
                     time_var = "Time", 
                     group_var = "Group",
                     adjst_var = c("Country")) 

var.grid <- expand.grid("Samp_Set" = resp.da.vars.ls$samp_set, 
                        "Taxa_lvl" = resp.da.vars.ls$used_ps, 
                        "Sub_Group" = resp.da.vars.ls$samp_gr, 
                        "Test_Var" = resp.da.vars.ls$resp_cols) %>% 
                 mutate(across(everything(), as.character))

#-------------------------------------------------------------------------------
# Add Response data to original metadata - will have to be moved to data prep.
#-------------------------------------------------------------------------------
resp.df <- resp.data.ls$Resp_data %>% 
                  bind_rows() %>% 
                  select(starts_with(c("Subject", "Resp")))

for(i in names(meta.ls)) {
  
  meta.ls[[i]] <- meta.ls[[i]] %>% 
                    left_join(., resp.df, by = "Subject") %>% 
                    'rownames<-'(.[["SeqID"]])
  
  
}


################################################################################
# LinaDA 
################################################################################
# Subjects with undetermiend response are removed from the dataset 

da.resp <- list()

linda.df.resp <- NULL  

for(i in 1:nrow(var.grid)) { 
  
  prm <- var.grid[i, ] %>% unlist()
 
  meta.inst <- meta.ls[[prm[["Samp_Set"]]]] %>% 
                  filter(.data[[resp.da.vars.ls$group_var]] == prm[["Sub_Group"]]) %>% 
                  filter(!is.na(.data[[prm[["Test_Var"]]]]))
  
  ps.inst <- pss.ls[[prm[["Samp_Set"]]]][[prm[["Taxa_lvl"]]]] %>% 
                  prune_samples(rownames(meta.inst), .)
  
  
  # Convert to mstat object
  mstat <- mStat_convert_phyloseq_to_data_obj(ps.inst)
  
  # Change to class 
  mstat$meta.dat <- meta.inst
  
  ld.res <- generate_taxa_trend_test_long(
                                    data.obj = mstat,
                                    subject.var = resp.da.vars.ls$subj_var,
                                    time.var = resp.da.vars.ls$time_var,
                                    group.var = prm[["Test_Var"]],
                                    adj.vars = resp.da.vars.ls$adjst_var,
                                    prev.filter = resp.da.vars.ls$min_prev,
                                    feature.mt.method = "fdr",
                                    feature.level = "original",
                                    feature.dat.type = "count")
  
  res.out <- list("Results" = ld.res$original[[1]], 
                  "Parameters" = prm, 
                  "Method" = "LinDA trend")
  
  da.resp[[paste0("DA_par", i)]] <- res.out
  
  linda.df.resp <- ld.res$original[[1]] %>% 
                          mutate(Tax_level = prm["Taxa_lvl"], 
                                 Group_level = prm["Sub_Group"],
                                 Subset = prm["Samp_Set"], 
                                 Response_value = prm["Test_Var"],
                                 Method ="LD_tr", 
                                 Prevalence_CutOff = resp.da.vars.ls$min_prev) %>% 
                          bind_rows(linda.df.resp, .)
                        
}

resp.da.res <- list("LinDA_ls" = da.resp, 
                    "LinDA_long_df" = linda.df.resp, 
                    "Param_table" = var.grid)

save(list = c("resp.da.res"), 
     file = "out/supp/5.1_Resp_DA.Rdata")
