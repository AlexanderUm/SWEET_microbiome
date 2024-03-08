#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "ggsignif", "broom", "vegan", 
               "metagenomeSeq", "qiime2R", "picante", "phangorn", "FSA", 
               "knitr", "MicrobiomeStat", "ggvegan", "ggpmisc", 
               "Maaslin2", "ComplexHeatmap", "circlize", "splinectomeR", "rfPermute")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")
load("out/supp/5_Responders.Rdata")

source("R/phy_taxa_filter.R")


################################################################################
# Outline
################################################################################
resp.vars.ls <- list(resp_cols = c("Resp_Body_weight"), 
                     samp_gr = c("Sugar", "S&SEs"), 
                     samp_set = c("all"), 
                     used_ps = c("Genus", "Family"), 
                     min_prev = 0.5,
                     subj_var = "Subject", 
                     time_var = "Time", 
                     tp_var = "CID",
                     group_var = "Group",
                     adjst_var = c("Country"), 
                     RF_min_prev = 0.25,
                     RF_sets = c("all", "CIDs_3"), 
                     RF_ps = c("css_Genus", "css_Family"), 
                     RF_time_p = c("CID_1", "CID_2", "CID_3", "CID_4"), 
                     RF_time_p_set = c("all"),
                     RF_ntree = 1001, 
                     RF_mtry = 12, 
                     RF_nperm = 99) 


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
# RF for combined dataset
################################################################################
var.grid.comb <- expand.grid("Samp_Set" = resp.vars.ls$RF_sets, 
                             "Taxa_lvl" = resp.vars.ls$RF_ps, 
                             "Sub_Group" = resp.vars.ls$samp_gr, 
                             "Test_Var" = resp.vars.ls$resp_cols) %>% 
                mutate(across(everything(), as.character))

rf.res.comb.ls <- list()

for(i in 1:nrow(var.grid.comb)) { 
  
  prm <- var.grid.comb[i, ] %>% unlist()
  
  #-----------------------------------------------------------------------------
  # Extract data for RF 
  #-----------------------------------------------------------------------------
  meta.inst <- meta.ls[[prm["Samp_Set"]]] %>% 
                  filter(.data[[resp.vars.ls$group_var]] == prm["Sub_Group"]) %>% 
                  filter(!is.na(.data[[prm[["Test_Var"]]]]))
  
  rf.df.inst <- pss.ls[[prm["Samp_Set"]]][[prm["Taxa_lvl"]]] %>% 
                  prune_samples(rownames(meta.inst), .) %>% 
                  phy_taxa_filter(., prev_fraction = resp.vars.ls$RF_min_prev) %>% 
                  otu_table() %>% 
                  as.matrix() %>% 
                  t() %>% 
                  as.data.frame() %>% 
                  mutate(!!prm["Test_Var"] := 
                           as.factor(meta.inst[[prm["Test_Var"]]])) %>% 
                  setNames(., gsub("\\(|\\)", "_", colnames(.)))
  
  s.size.b <- balancedSampsize(meta.inst[[prm["Test_Var"]]])
  
  #-----------------------------------------------------------------------------
  # Run rundom forest 
  #-----------------------------------------------------------------------------
  set.seed(349675)
  
  rf.perm.res <- rfPermute(as.formula(paste0(prm["Test_Var"], " ~ .")), 
                           data = rf.df.inst, 
                           ntree = resp.vars.ls$RF_ntree, 
                           # mtry= resp.vars.ls$RF_mtry, 
                           num.cores = 4, 
                           num.rep = resp.vars.ls$RF_nperm, 
                           sampsize = s.size.b, 
                           proximity = TRUE)
  
  res.ls <- list("RF_res" = rf.perm.res, 
                 "meta"= meta.inst, 
                 "RF_data" = rf.df.inst, 
                 "Parameters" = prm)
  
  rf.res.comb.ls[[paste0("comb_par", i)]] <- res.ls
  
                
}

################################################################################
# Per time point 
################################################################################
var.grid.tp <- expand.grid("Samp_Set" = resp.vars.ls$RF_time_p_set, 
                           "Taxa_lvl" = resp.vars.ls$RF_ps, 
                           "Sub_Group" = resp.vars.ls$samp_gr, 
                           "Test_Var" = resp.vars.ls$resp_cols, 
                           "CID" = resp.vars.ls$RF_time_p) %>% 
               mutate(across(everything(), as.character))

rf.res.tp.ls <- list()

for(i in 1:nrow(var.grid.tp)) { 
  
  prm <- var.grid.tp[i, ] %>% unlist()
  
  #-----------------------------------------------------------------------------
  # Extract data for RF 
  #-----------------------------------------------------------------------------
  meta.inst <- meta.ls[[prm["Samp_Set"]]] %>% 
                  filter(.data[[resp.vars.ls$group_var]] == prm["Sub_Group"], 
                         .data[[resp.vars.ls$tp_var]] == prm["CID"]) %>% 
                  filter(!is.na(.data[[prm[["Test_Var"]]]]))
  
  rf.df.inst <- pss.ls[[prm["Samp_Set"]]][[prm["Taxa_lvl"]]] %>% 
                    prune_samples(rownames(meta.inst), .) %>% 
                    phy_taxa_filter(., prev_fraction = resp.vars.ls$RF_min_prev) %>% 
                    otu_table() %>% 
                    as.matrix() %>% 
                    t() %>% 
                    as.data.frame() %>% 
                    mutate(!!prm["Test_Var"] := 
                             as.factor(meta.inst[[prm["Test_Var"]]])) %>% 
                    setNames(., gsub("\\(|\\)", "_", colnames(.)))
  
  s.size.b <- balancedSampsize(meta.inst[[prm["Test_Var"]]])
  
  set.seed(349675)
  
  rf.perm.res <- rfPermute(as.formula(paste0(prm["Test_Var"], " ~ .")), 
                           data = rf.df.inst, 
                           ntree = resp.vars.ls$RF_ntree, 
                           # mtry= resp.vars.ls$RF_mtry, 
                           num.cores = 4, 
                           num.rep = resp.vars.ls$RF_nperm, 
                           sampsize = s.size.b, 
                           proximity = TRUE)
  
  res.ls <- list("RF_res" = rf.perm.res, 
                 "meta"= meta.inst, 
                 "RF_data" = rf.df.inst, 
                 "Parameters" = prm)
  
  rf.res.tp.ls[[paste0("tp_par", i)]] <- res.ls
  
}

resp.RF.res <- list("RF_comb" = rf.res.comb.ls, 
                    "RF_per_tp" = rf.res.tp.ls, 
                    "RF_Param_table_comb" = var.grid.comb, 
                    "RF_Param_table_tp" = var.grid.tp)

save(list = c("resp.RF.res"), 
     file = "out/supp/5.2_Resp_RF.Rdata")
