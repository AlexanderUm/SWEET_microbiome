#-------------------------------------------------------------------------------
# List of parameter to use 
#-------------------------------------------------------------------------------
PRM <- list()

PRM[["general"]] <- list("cols" = list("GroupFac"= "Group", 
                                       "TimeNum" = "Time", 
                                       "CountryFac" = "Country", 
                                       "SubjectFac" = "Subject"), 
                         "libs" = c("phyloseq", "tidyverse", "broom", "vegan", 
                                    "metagenomeSeq", "qiime2R",
                                    "MicrobiomeStat", "ggvegan", 
                                    "ComplexHeatmap", 
                                    "rfPermute", "caret", "randomForest", 
                                    "pROC", "cowplot", "lme4"), 
                         "seed" = 394573, 
                         "round_to" = 4)


# Data preparation
PRM[["data"]] <- list("out_dir" = "out/Rdata",
                      "qiime_path" = "data/", 
                      "meta_path" = "data/meta_data.csv", 
                      "picrust_path" = c("MetaCys" = "data/path_abun_unstrat.tsv"),
                      "min_reads_per_taxa" = 50, 
                      "glom_lvls" = c("Genus", "Family"), 
                      "glom_NArm" = TRUE,
                      "rare_depth" = 40000,
                      "sample_subs" = list("no_M0" = c("CID_2", "CID_3", "CID_4")),
                      "count_norm" = c("CSS_log2", "TSS_log2", "Rare"))


# Alpha diversity 
PRM[["alpha"]] <- list("out_dir" = "out/alpha_diversity", 
                       "alpha_ind" = c("shannon", "simpson", 
                                       "observed_species", "chao1"), 
                       "taxa_lvl" = "ASV_alpha", 
                       "data_set" = c("all"),
                       "count_norm" = c("Rare"), 
                       "time_var" = "Time", 
                       "subject_var" = "Subject", 
                       "group_var" = "Group", 
                       "adjust_var" = "Country")


PRM[["beta"]] <- list("out_dir" = "out/beta_diversity", 
                      "taxa_lvl" = c("ASV"), 
                      "data_set" = c("all"),
                      "count_norm" = c("TSS_log2"),
                      "dists" = c("Unweighted UniFrac" = "unifrac", 
                                  "Weighted UniFrac" = "wunifrac", 
                                  "Jaccard" = "jaccard", 
                                  "Bray-Curtis" = "bray"),
                      "n_perm" = 999,
                      "formula_RDA" = "Time*Group + Condition(Country)",
                      "p_color" = "Group", 
                      "test_var" = "Group",
                      "p_shape" = "CID", 
                      "p_group" = "Subject")


PRM[["DA"]] <- list("out_dir" = "out/DA", 
                    "taxa_lvl" = c("Genus", "Family", "MetaCys"),
                    "data_set" = c("all"),
                    "count_norm" = c("count"),
                    "count_norm_plot" = c("CSS_log2"),
                    "time_var" = "Time", 
                    "subject_var" = "Subject", 
                    "group_var" = "Group", 
                    "adjust_var" = "Country", 
                    "tax_min_prev" = 0.5, 
                    "max_q-val" = 0.1)


# Responders/non-responders
PRM[["resp"]] <- list("out_dir" = "out/Response")

PRM[["resp"]][["calc"]] <- list("Body_weight" = 
                                  list("index_fun" = "index_wms_fun", 
                                       "resp_fun" = "resp_5tile_fun"), 
                                "Fasting_glucose" = 
                                  list("index_fun" = "index_cid4_cid2_fun", 
                                       "resp_fun" = "resp_5tile_fun"), 
                                "HbA1c" = 
                                  list("index_fun" = "index_cid4_cid2_fun", 
                                       "resp_fun" = "resp_morethan0_fun"))

PRM[["resp"]][["RF"]] <- list("data_set_cross" = c("CID_1", "CID_2", 
                                                   "CID_3", "CID_4"),
                              "taxa_lvl" = c("Genus"), 
                              "group_lvls" = c("S&SEs", "Sugar"), 
                              "count_norm" = c("TSS_log2"),
                              "tax_min_prev" = 0.5, 
                              "n_random" = 10,
                              "time_var" = "Time", 
                              "subject_var" = "Subject", 
                              "group_var" = "Group", 
                              "n_trees" = 1999, 
                              "n_cv" = 5, 
                              "cv_repeats" = 25, 
                              "n_rf_permute" = 199, 
                              "imp_max_pval" = 0.05)

save(list = c("PRM"), file = "PRM.Rdata")

rm(list = ls())
gc()
