#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
load("PRM.Rdata")

load(paste0(PRM$data$out_dir, "/0_data.Rdata"))
load(paste0(PRM$data$out_dir, "/4_Resp.Rdata"))

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE)}

# Custom functions
source("R/phy_taxa_filter.R")
source("R/RF_supp_functions.R")

# Create directory 
DirOut <- paste0(PRM$resp$out_dir, "/RF")

dir.create(paste0(DirOut, "/plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(DirOut, "/tabs"), recursive = TRUE, showWarnings = FALSE)


################################################################################
# Cross-sectional (per CID)
################################################################################
# Variables grid 
RespGrid <- expand.grid("taxa_lvl" = PRM$resp$RF$taxa_lvl, 
                        "group_lvls" = PRM$resp$RF$group_lvls, 
                        "resp_vars" = paste0("Resp_", 
                                             names(PRM$resp$calc), 
                                             "_index"),
                        "count_norm" = PRM$resp$RF$count_norm, 
                        "samples_set" = PRM$resp$RF$data_set_cross,
                        stringsAsFactors = FALSE) 

ResRfCombDf <- NULL

RocDataDf <- NULL

for(i in 1:nrow(RespGrid)) {  
  
  iSampleSet <- RespGrid[i, "samples_set"]
  
  iTaxaLvl <- RespGrid[i, "taxa_lvl"]
  
  iCountNorm <- RespGrid[i, "count_norm"]
  
  iGroupLvl <- RespGrid[i, "group_lvls"]
  
  iRespVar <- RespGrid[i, "resp_vars"]
  
  
  # Extract corresponding metadata and filter out NAs 
  Meta <- RespData$Resp_meta_comb$all %>% 
             filter(!is.na(.data[[iRespVar]]))
  
  # Leave one focus group 
  MetaFilt <- Meta %>% 
               filter(.data[[PRM$resp$RF$group_var]] == iGroupLvl)
  
  # Extract phyloseq, remove NA group samples filter out taxa with less than x prevalence 
  Ps <- DataComb[["all"]][["ps"]][[iTaxaLvl]][[iCountNorm]] %>% 
            prune_samples(rownames(Meta), .) %>% 
            phy_taxa_filter(prev_fraction = PRM$resp$RF$tax_min_prev, 
                            group_col = PRM$resp$RF$group_var) %>% 
            prune_samples(rownames(MetaFilt), .)
  
  # Extract OTU table & filter
  DataRf <- Ps %>% 
            otu_table() %>% 
            as.matrix() %>% 
            t() %>% 
            as.data.frame() %>% 
            bind_cols(MetaFilt) %>% 
            filter(CID %in% iSampleSet) %>% 
            select(all_of(c(taxa_names(Ps), iRespVar))) %>% 
            rename(Response = iRespVar)
  
  # Vector with the true response
  RespTrue <- DataRf$Response
  
  
  #-----------------------------------------------------------------------------
  # RF with CV
  #-----------------------------------------------------------------------------
  # Mtry
  SeedMtry <- round(sqrt((ncol(DataRf)-1)), 0) 
  
  MtryGrid <- expand.grid(.mtry = SeedMtry)
  
  
  # Control
  Control <- trainControl(method = "repeatedcv", 
                          repeats = PRM$resp$RF$cv_repeats,
                          number = PRM$resp$RF$n_cv, 
                          classProbs = TRUE, 
                          savePredictions = "all", 
                          summaryFunction = twoClassSummary) 
  
  # Model
  RfCvRes <- train(Response ~ .,
                   data = DataRf,
                   method = "rf",
                   metric = "ROC",
                   importance=TRUE, 
                   tuneGrid=MtryGrid, 
                   trControl = Control,
                   ntree = PRM$resp$RF$n_trees)
  
  # Data for ROC
  RocDataDf <- roc(RfCvRes$pred$obs, 
                   RfCvRes$pred$Resp)[c("sensitivities", "specificities")] %>% 
                as.data.frame() %>% 
                mutate(Random = "True", 
                       Correspondance = 100) %>% 
                arrange(-row_number()) %>% 
                bind_cols(., RespGrid[i, ]) %>% 
                bind_rows(RocDataDf, .) %>% 
                suppressMessages()
  
  # Results table
  ResRfCombDf <- RfCvRes$results %>% 
                    bind_cols(., RespGrid[i, ]) %>% 
                    mutate(Random = "True", 
                           Correspondance = 100, 
                           Mtry = SeedMtry) %>% 
                    bind_rows(ResRfCombDf, .)
  
  
  #!!!--------------------------------------------------------------------------
  #!!! Randomize response
  #!!!--------------------------------------------------------------------------
  for(j in 1:PRM$resp$RF$n_random) {
    
    # Randomize response with sample
    DataRf <- DataRf %>% 
                mutate(Response = sample(Response))
    
    # Correspondence between randomized and true response 
    CoresPrec <- sum(DataRf$Response == RespTrue)/length(RespTrue)*100 %>% 
                    round(2)
    
    # Run the model (control from above)
    RfCvRes <- train(Response ~ .,
                     data = DataRf,
                     method = "rf",
                     metric = "ROC",
                     importance=TRUE, 
                     tuneGrid=MtryGrid, 
                     trControl = Control,
                     ntree = PRM$resp$RF$n_trees)
    
    # Add to ROC data frame 
    RocDataDf <- roc(RfCvRes$pred$obs, 
                     RfCvRes$pred$Resp)[c("sensitivities", "specificities")] %>% 
                  as.data.frame() %>% 
                  mutate(Random = paste0("Random_", j), 
                         Correspondance = CoresPrec) %>% 
                  bind_cols(., RespGrid[i, ]) %>% 
                  bind_rows(RocDataDf, .) %>% 
                  suppressMessages()
    
    # Add to results data frame 
    ResRfCombDf <- RfCvRes$results %>% 
                    bind_cols(., RespGrid[i, ]) %>% 
                    mutate(Random = paste0("Random_", j), 
                           Correspondance = CoresPrec, 
                           Mtry = SeedMtry) %>% 
                    bind_rows(ResRfCombDf, .)
    
  }
  
}


################################################################################
# Visualization
################################################################################
#-------------------------------------------------------------------------------
# Plot ROC curves
#-------------------------------------------------------------------------------
RocDataDf <- RocDataDf %>% 
                mutate(Order = gsub("_.*", "", Random), 
                       resp_vars = gsub("Resp_", "", resp_vars)) %>% 
                mutate(resp_vars = gsub("_", " ", resp_vars), 
                       x_grid_var = paste0(group_lvls, " [", samples_set, "]"))

RocPlot <- ggplot(RocDataDf, aes(x = sensitivities, 
                                 y = specificities)) + 
                geom_line(aes(group = Random, 
                              color = Order, 
                              alpha = Correspondance), 
                          linewidth = 0.75) + 
                facet_grid(c("resp_vars", "x_grid_var")) + 
                scale_x_reverse() + 
                theme_bw() +
                theme(panel.grid = element_blank(),
                      strip.text = element_text(size = 7, face = "italic"),
                      axis.title = element_blank(), 
                      axis.text.x = element_text(angle = 45, hjust = 1), 
                      legend.position = "bottom") + 
                geom_abline(intercept = 1) + 
                scale_alpha_continuous(guide = "none")

ggsave(paste0(DirOut, "/plots/ROC_comb.png"), plot = RocPlot, 
       width = 12, height = 5)


#-------------------------------------------------------------------------------
# Plot AUC for original and random
#-------------------------------------------------------------------------------
ResRfCombDf <- ResRfCombDf %>% 
                  mutate(Order = gsub("_.*", "", Random)) %>% 
                  mutate(resp_vars = gsub("Resp_", "", resp_vars)) %>% 
                  mutate(resp_vars = gsub("_", " ", resp_vars))

AucPlot <- ggplot(ResRfCombDf, aes(y = ROC, x = group_lvls)) + 
                  geom_jitter(aes(colour = Order), 
                              height = 0, width = 0.05, 
                              size = 2.5, alpha = 0.5) + 
                  facet_grid(c("resp_vars", "samples_set")) + 
                  theme_bw() +
                  theme(panel.grid = element_blank(),
                        strip.text = element_text(size = 8, face = "italic"),
                        axis.title = element_blank(), 
                        axis.text.x = element_text(angle = 45, hjust = 1), 
                        legend.position = "bottom")

ggsave(paste0(DirOut, "/plots/AUC_comb.png"), 
       plot = AucPlot, 
       width = 7, height = 6)


################################################################################
# Taxa significance
################################################################################
# Sets with the reasonable prediction 
SetsForImp <- list(c("resp_vars" = "Resp_HbA1c_index", 
                     "samples_set" = "CID_1", 
                     "group_lvls" = "S&SEs"), 
                   c("resp_vars" = "Resp_HbA1c_index", 
                     "samples_set" = "CID_2", 
                     "group_lvls" = "S&SEs"), 
                   c("resp_vars" = "Resp_Body_weight_index", 
                     "samples_set" = "CID_2", 
                     "group_lvls" = "S&SEs"))

ImpRfPlotsLs <- list()

ImpRfTabsLs <- list()

ImpRfResLs <- list()

for(i in 1:length(SetsForImp)) {
  
  iSet <- SetsForImp[[i]]
  
  iSetName <- paste0("set", i)
  
  # Extract metadata
  Meta <- RespData$Resp_meta_comb$all %>%
            filter(!is.na(.data[[iSet[["resp_vars"]]]]),
                   CID == iSet[["samples_set"]]) %>%
            droplevels()
  
  # Leave one focus group
  MetaFilt <- Meta %>%
                filter(Group == iSet[["group_lvls"]])
  
  
  # Extract phyloseq, remove NA group samples filter out taxa with less than x prevalence
  DataRf <- DataComb[["all"]][["ps"]][["Genus"]][[PRM$resp$RF$count_norm]] %>%
                prune_samples(rownames(Meta), .) %>%
                phy_taxa_filter(prev_fraction = PRM$resp$RF$tax_min_prev,
                                group_col = PRM$resp$RF$group_var) %>%
                prune_samples(rownames(MetaFilt), .) %>%
                otu_table() %>%
                as.matrix() %>%
                t() %>%
                as.data.frame() %>%
                mutate(Response = MetaFilt[[iSet[["resp_vars"]]]])
  
  # Define mtry
  SeedMtry <- round(sqrt((ncol(DataRf)-1)), 0)
  
  #-----------------------------------------------------------------------------
  # Run RF permute
  ImpRfResLs[[iSetName]] <- rfPermute(Response ~ .,
                                      data = DataRf,
                                      ntree = PRM$resp$RF$n_trees,
                                      mtry=SeedMtry,
                                      num.cores = 4,
                                      num.rep = PRM$resp$RF$n_rf_permute,
                                      proximity = TRUE)
  
  #-----------------------------------------------------------------------------
  # Extract results and plot taxa
  ImpRfTabsLs[[iSetName]] <- RF_extract_importance(ImpRfResLs[[iSetName]]) %>%
                  select(all_of(c("feature",
                                  "Resp_importance",
                                  "NonResp_importance",
                                  "MeanDecreaseAccuracy_importance",
                                  "NonResp.unscaled_pval",
                                  "Resp.unscaled_pval",
                                  "MeanDecreaseAccuracy.unscaled_pval"))) %>%
                  filter(if_any(ends_with("_pval"), 
                                function(x){x <= PRM$resp$RF$imp_max_pval})) %>%
                  setNames(gsub("MeanDecreaseAccuracy", "MDA", colnames(.)))
  
  
  # Importance plots with a custom funciton
  ImpRfPlotsLs[[iSetName]] <- RF_plot_importance(ImpRfTabsLs[[iSetName]], 
                                                 cols_to_plot = c("Resp_importance", 
                                                                  "NonResp_importance", 
                                                                  "MDA_importance"), 
                                                 col_arrange_by = "MDA_importance") + 
                                      theme(legend.position = "bottom") +
                                      scale_y_discrete(position = "right")
}


#-------------------------------------------------------------------------------
# Combine Importance and AUC plots 
#-------------------------------------------------------------------------------
# Keep ROC data only for the important models 
RocDataDfFilt <- NULL

for(i in 1:length(SetsForImp)) {
  
  iSet <- SetsForImp[[i]]
  
  RocDataDfFilt <- RocDataDf %>% 
                    filter(.data[[names(iSet)[1]]] == gsub("_", " ", 
                                                           gsub("Resp_", 
                                                                "", iSet[[1]])), 
                           .data[[names(iSet)[2]]] == iSet[[2]],
                           .data[[names(iSet)[3]]] == iSet[[3]],) %>% 
                    bind_rows(RocDataDfFilt, .)
}

# Rearrange factors 
RocDataDfFilt <- RocDataDfFilt %>% 
                  mutate(samples_set = factor(samples_set, 
                                              levels = unique(samples_set)), 
                         resp_vars = factor(resp_vars, 
                                            levels = unique(resp_vars)))

# ROC curve plots 
AucPlotStrong <- ggplot(RocDataDfFilt, 
                        aes(x = sensitivities, y = specificities)) + 
                  geom_line(aes(group = Random, 
                                color = Order, 
                                alpha = Correspondance), 
                            linewidth = 1) + 
                  facet_wrap(c("resp_vars", "samples_set"), scales = "free") + 
                  scale_x_reverse() + 
                  theme_bw() +
                  theme(panel.grid = element_blank(),
                        strip.text = element_text(size = 10, face = "italic"),
                        axis.title = element_blank(), 
                        axis.text.x = element_text(angle = 45, hjust = 1), 
                        legend.position = "right", 
                        panel.spacing = unit(100, "pt")) + 
                  geom_abline(intercept = 1) + 
                  scale_alpha_continuous(guide = "none") + 
                  scale_y_continuous(position = "right")


# Remove legend 
ImpRfPlotsLsNl <- lapply(ImpRfPlotsLs, 
                         function(x){x + theme(legend.position = "none")})

# Get legend 
LegendP <- get_legend((ImpRfPlotsLs[[1]] + 
                         theme(legend.position = "right")))

# Arrange ROC and Importance plots in the grid
ImpPlotGrid <- plot_grid(plotlist = ImpRfPlotsLsNl, nrow = 1)

ImpPlotGridL <- plot_grid(ImpPlotGrid, LegendP, rel_widths = c(0.9, 0.1))

AucPlotGrid <- plot_grid(AucPlotStrong , NULL, rel_widths = c(0.9, 0.1))

ImpFinalFig <- plot_grid(AucPlotGrid, ImpPlotGridL, ncol = 1, 
                         rel_heights = c(0.35, 0.65), labels = "AUTO", scale = 0.975)


ggsave(paste0(DirOut, "/plots/FinalComb.png"), plot = ImpFinalFig, 
       width = 15, height = 8)


#-------------------------------------------------------------------------------
# Collect and back up results 
#-------------------------------------------------------------------------------
ResRFRes_04 <- list("RocData" = RocDataDf, 
                    "RfData" = ResRfCombDf, 
                    "RfPermuteRes" = ImpRfResLs,
                    "RfPermuteSigRes" = ImpRfTabsLs,
                    "ImpSetsInfo" = SetsForImp)


save(list = c("ResRFRes_04"),
     file = paste0(PRM$data$out_dir, "/4_RespRF.Rdata"))

# rm(list = ls())
# gc()

