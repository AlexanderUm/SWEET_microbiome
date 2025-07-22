#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
load("PRM.Rdata")

load(paste0(PRM$data$out_dir, "/0_data.Rdata"))

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE)}

# Create directory 
DirOut <- PRM$beta$out_dir

dir.create(paste0(DirOut, "/plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(DirOut, "/tabs"), recursive = TRUE, showWarnings = FALSE)

# Custom functions
source("R/phy_dists_ls.R")

#-------------------------------------------------------------------------------
# Variables grid for loops 
#-------------------------------------------------------------------------------
BetaRes <- list()

GridBeta <- expand.grid("taxa_lvl" = PRM$beta$taxa_lvl, 
                        "count_norm" = PRM$beta$count_norm, 
                        stringsAsFactors = FALSE) 

GridRda <- expand.grid("dist" = names(PRM$beta$dists),
                       "formula" = PRM$beta$formula_RDA, 
                       "samples_set" = PRM$beta$data_set,
                       stringsAsFactors = FALSE) 


for(i in 1:nrow(GridBeta)) {
  
  # Create sub-directories 
  for(dir in c("/plots/", "/tabs/")) {
    
    dir.create(paste0(DirOut, dir, paste(GridBeta[i, ], collapse = "-")), 
               recursive = TRUE, showWarnings = FALSE)
    
  }
  
  # Variables 
  iTaxaLvl <- GridBeta[i, "taxa_lvl"]
  
  iCountNorm <- GridBeta[i, "count_norm"]
  
  # Extract data for all samples
  Ps <- DataComb[["all"]]$ps[[iTaxaLvl]][[iCountNorm]] 
  
  #-----------------------------------------------------------------------------
  # Calculate distances - custom function 
  #-----------------------------------------------------------------------------
  DistsLs <- phy_dist_ls(Ps, dists = PRM$beta$dists) %>% 
              setNames(names(PRM$beta$dists))
  
  
  #-----------------------------------------------------------------------------
  # Build RDAs & Plot 
  #-----------------------------------------------------------------------------
  StatResComb <- NULL
  
  for(j in 1:nrow(GridRda)) {
    
    jSampleSet <- GridRda[j, "samples_set"]
    
    jDist <- GridRda[j, "dist"]
    
    jFormula <- GridRda[j, "formula"]
    
    jTitle <- paste0("Distance: ", jDist)
    
    Meta <- DataComb[[jSampleSet]][["meta"]]
    
    Dist <- DistsLs[[jDist]] %>% 
              usedist::dist_subset(idx = rownames(Meta))
    
    # Built and test RDA
    RdaObj <- dbrda(as.formula(paste0("Dist ~ ", jFormula)), 
                    data = Meta)
    
    RdaStat <- anova.cca(RdaObj,
                         permutations = how(nperm = PRM$beta$n_perm),
                         by = "terms")
    
    
    #---------------------------------------------------------------------------
    # Plot RDA
    #---------------------------------------------------------------------------
    # Summarize RDA object 
    RdaObjSummary <- summary(RdaObj)
    
    AxisExpl <- round(RdaObjSummary$cont$importance[2, 1:2]*100, 2)
    
    # Extract scores
    AxisScore <- scores(RdaObj, display = "sites")[, 1:2] %>%
                            as.data.frame() %>%
                            bind_cols(Meta) %>% 
                            arrange(Subject, Time)
    
    # Extract centroids and vectors
    RdaAutoPlot <- autoplot(RdaObj) %>%
                        ggplot_build()
    
    CentrDf <- RdaAutoPlot$data[[length(RdaAutoPlot$data)]] %>%
                        mutate(Centroids = gsub(PRM$beta$test_var, "", label))
    
    VectorDf <- RdaAutoPlot$data[[(length(RdaAutoPlot$data)-1)]] %>%
                        mutate(Vector = gsub(PRM$beta$test_var, "", label)) %>% 
                        mutate(Vector = factor(Vector, levels = rev(Vector)))
    
    # Plot 
    BetaPlot <- ggplot() +
                geom_path(data = AxisScore,
                          aes(x = dbRDA1, y = dbRDA2,
                              color = .data[[PRM$beta$p_color]],
                              group = .data[[PRM$beta$p_group]]),
                          alpha = 0.1) +
                geom_point(data = AxisScore,
                           aes(x = dbRDA1, y = dbRDA2,
                               color = .data[[PRM$beta$p_color]],
                               shape = .data[[PRM$beta$p_shape]]),
                           alpha = 0.5) +
                geom_segment(data = VectorDf,
                             aes(x = 0, y = 0, 
                                 xend = x, yend = y, 
                                 linetype = Vector),
                             lineend = "round", 
                             linejoin = "round",
                             arrow = arrow(length = unit(0.2, "cm"),  
                                           type = "closed"), 
                             linewidth = 0.75) +
                scale_color_manual(values = PlotsAttr$color_gr) +
                scale_shape_manual(values = PlotsAttr$shape_cid) +
                scale_linetype_manual(values = c(1, 4)) +
                coord_fixed() +
                theme_bw() +
                xlab(paste0(colnames(AxisScore)[1], " [", AxisExpl[1], "%]")) +
                ylab(paste0(colnames(AxisScore)[2], " [", AxisExpl[2], "%]")) +
                ggtitle(jTitle) + 
                guides(color = guide_legend(order=1),
                       fill = guide_legend(order=2),
                       shape = guide_legend(order=3)) + 
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
    
    
    #---------------------------------------------------------------------------
    # Collect results 
    #---------------------------------------------------------------------------
    BetaRes[[jSampleSet]][[iTaxaLvl]][[iCountNorm]][[jDist]][["plot"]] <- BetaPlot
    
    RdaStatDf <- RdaStat %>% 
                    as.data.frame() %>% 
                    mutate(Explained = SumOfSqs/sum(SumOfSqs)) %>% 
                    mutate(across(everything(), function(x){round(x, 5)})) %>% 
                    bind_cols(GridBeta[i, ], GridRda[j, ]) %>% 
                    rownames_to_column(var = "Terms")
    
    StatResComb <- RdaStatDf %>% 
                    add_row() %>% 
                    bind_rows(StatResComb, .)
    
    write.csv(RdaStatDf, 
              file = paste0(DirOut, "/tabs/", 
                            paste(GridBeta[i, ], collapse = "-"), 
                                  "/", gsub(" ", "_", jDist), "--", 
                                  jSampleSet, "--stat.csv"), 
              row.names = FALSE, na = "") 
    
    BetaRes[[jSampleSet]][[iTaxaLvl]][[iCountNorm]][[jDist]][["stat"]] <- RdaStatDf
    
    
    ggsave(filename = paste0(DirOut, "/plots/", 
                             paste(GridBeta[i, ], collapse = "-"), 
                             "/", gsub(" ", "_", jDist), 
                             "--", jSampleSet, "--dbRDA.svg"), 
           plot = BetaPlot, width = 5, height = 3.5)
    
  }
  
  write.csv(StatResComb, 
            file = paste0(DirOut, "/tabs/", 
                          paste(GridBeta[i, ], collapse = "-"), 
                          "/Combined_stat.csv"), 
            row.names = FALSE, na = "") 
  
}

#-------------------------------------------------------------------------------
# Plot for publication 
#-------------------------------------------------------------------------------
BetaPlotLs <- list(BetaRes$all$ASV$TSS_log2$Jaccard$plot,
                    BetaRes$all$ASV$TSS_log2$`Bray-Curtis`$plot,
                    BetaRes$all$ASV$TSS_log2$`Unweighted UniFrac`$plot,
                    BetaRes$all$ASV$TSS_log2$`Weighted UniFrac`$plot)


BetaPlotLs <- lapply(BetaPlotLs, 
                     function(x){x + theme(legend.position = "none")})

PlotLegend <- get_legend(BetaRes$all$ASV$TSS_log2$Jaccard$plot)

PlotGrid0 <- plot_grid(plotlist = BetaPlotLs)

FullGrid <- plot_grid(PlotGrid0, NULL,
                      PlotLegend, 
                      ncol = 3, 
                      rel_widths = c(1, 0.001, 0.1))

save_plot(filename = paste0(DirOut, "/plots/Figure_3.tiff"), 
          compression = "lzw",
          dpi = 300,
          plot = FullGrid, 
          base_width = 10.5, 
          base_height = 7.5) 

#-------------------------------------------------------------------------------
# Save and clean up 
#-------------------------------------------------------------------------------
save(list = c("BetaRes"), 
     file = paste0(PRM$data$out_dir, "/2_beta.Rdata"))

rm(list = ls())
gc()
