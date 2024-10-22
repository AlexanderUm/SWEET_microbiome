#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
load("PRM.Rdata")

load(paste0(PRM$data$out_dir, "/0_data.Rdata"))

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE)}

# Create directory 
DirOut <- PRM$DA$out_dir

dir.create(paste0(DirOut, "/plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(DirOut, "/tabs"), recursive = TRUE, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# Variables grid for loops 
#-------------------------------------------------------------------------------
DaRes <- list()

DaResDf <- NULL

SigTaxaLs <- list()

GridDa <- expand.grid("taxa_lvl" = PRM$DA$taxa_lvl, 
                      "count_norm" = PRM$DA$count_norm, 
                      "samples_set" = PRM$DA$data_set,
                       stringsAsFactors = FALSE) 


for(i in 1:nrow(GridDa)) {
  
  # Variables 
  iTaxaLvl <- GridDa[i, "taxa_lvl"]
  
  iCountNorm <- GridDa[i, "count_norm"]
  
  iSampleSet <- GridDa[i, "samples_set"]
  
  
  # Extract data 
  Ps <- DataComb[[iSampleSet]][["ps"]][[iTaxaLvl]][[iCountNorm]] 
  
  Meta <- DataComb[[iSampleSet]][["meta"]]
  
  
  #-----------------------------------------------------------------------------
  # LinDA
  #-----------------------------------------------------------------------------
  # Convert to msata object 
  MstatObj <- mStat_convert_phyloseq_to_data_obj(Ps)
  
  # Adjust metadata 
  MstatObj$meta.dat <- Meta 
  
  # Run trend test 
  LinDaRes <- generate_taxa_trend_test_long(
                      data.obj = MstatObj,
                      subject.var = PRM$DA$subject_var,
                      time.var = PRM$DA$time_var,
                      group.var = PRM$DA$group_var,
                      adj.vars = PRM$DA$adjust_var,
                      prev.filter = PRM$DA$tax_min_prev,
                      feature.level = "original",
                      feature.dat.type = "count")
  
  # Format & collect results 
  LinDaResDf <- LinDaRes$original[[1]] %>% 
                  bind_cols(GridDa[i, ]) %>% 
                       mutate(MinPrev = PRM$DA$tax_min_prev)
  
  DaRes[[iTaxaLvl]][[iSampleSet]][[iCountNorm]] <- LinDaResDf
  
  DaResDf <- bind_rows(DaResDf, LinDaResDf)
  
  LinDaResDfSig <- LinDaResDf %>% 
                      filter(Adjusted.P.Value <= PRM$DA$`max_q-val`) %>% 
                      select(Variable, P.Value, Adjusted.P.Value, Coefficient)
  
  SigTaxaLs[[iTaxaLvl]][[iSampleSet]][[iCountNorm]] <- LinDaResDfSig
  
  # Write CSVs
  dir.create(paste0(DirOut, "/tabs/", iTaxaLvl, "/"), 
                    recursive = TRUE, showWarnings = FALSE)
  
  write.csv(LinDaResDf, 
            paste0(DirOut, "/tabs/", iTaxaLvl, "/", 
                   iSampleSet, "--", iCountNorm, ".csv"))
  
  write.csv(LinDaResDfSig, 
            paste0(DirOut, "/tabs/", iTaxaLvl, "/",
                   iSampleSet, "--", iCountNorm, "--Sig.csv"))
  
}


################################################################################
#-------------------------------------------------------------------------------
# Plot results 
#-------------------------------------------------------------------------------
PubSigTaxLs <- list()

PubSigTaxLs[["Genus"]] <- read.csv(file = "data/sig_taxa_pub.csv")$x

PubSigTaxLs[["Family"]] <- read.csv(file = "data/sig_taxa_pub_supp.csv")$Variable


GridDaPlot <- expand.grid("taxa_lvl" = PRM$DA$taxa_lvl, 
                          "count_norm_plot" = PRM$DA$count_norm_plot, 
                          "samples_set" = PRM$DA$data_set,
                          stringsAsFactors = FALSE) 

SpaghPlotsLs <- list()

for(i in 1:nrow(GridDaPlot)) {
  
  # Variables 
  iTaxaLvl <- GridDaPlot[i, "taxa_lvl"]
  
  iCountNorm <- GridDaPlot[i, "count_norm_plot"]
  
  iSampleSet <- GridDaPlot[i, "samples_set"]
  
  iDirOut <- paste0(DirOut, "/plots/", iTaxaLvl)
  
  # Make a directory for plots 
  dir.create(iDirOut, recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(iDirOut, "/", iSampleSet, "--", iCountNorm), 
             recursive = TRUE, showWarnings = FALSE)
  
  # Extract data 
  Meta <- DataComb[[iSampleSet]][["meta"]] 
  
  # Significant taxa
  SigTaxa <- DaResDf %>% 
              filter(taxa_lvl == iTaxaLvl, 
                     samples_set == iSampleSet, 
                     Adjusted.P.Value <= PRM$DA$`max_q-val`)
  
  
  if(nrow(SigTaxa) > 0 ) {
    
    #---------------------------------------------------------------------------
    # Spaghetti plots
    #---------------------------------------------------------------------------
    TaxaLongDf <- DataComb[[iSampleSet]][["ps"]][[iTaxaLvl]][[iCountNorm]] %>% 
                      prune_taxa(SigTaxa$Variable, .) %>% 
                      otu_table() %>% 
                      as.matrix() %>% 
                      t() %>% 
                      as.data.frame() %>% 
                      bind_cols(Meta) %>% 
                      pivot_longer(cols = all_of(SigTaxa$Variable), 
                                   names_to = "Taxa", 
                                   values_to = "Taxa_Abundance") %>% 
                      mutate(Mean = mean(Taxa_Abundance), 
                             Median = median(Taxa_Abundance),
                             .by = c(PRM$DA$group_var, "Taxa", PRM$DA$time_var))
    
    FullGraph <- ggplot(TaxaLongDf, 
                   aes(x = .data[[PRM$DA$time_var]])) + 
                    geom_smooth(aes(x = .data[[PRM$DA$time_var]], 
                                    y = Taxa_Abundance, 
                                    fill = .data[[PRM$DA$group_var]], 
                                    color = .data[[PRM$DA$group_var]]), 
                                method = "lm") +
                   geom_line(aes(y = Mean, 
                                  group = .data[[PRM$DA$group_var]], 
                                  linetype = .data[[PRM$DA$group_var]])) + 
                   geom_point(aes(y = Mean, 
                                   shape = .data[[PRM$DA$group_var]])) + 
                   scale_x_continuous(breaks = sort(unique(TaxaLongDf[[PRM$DA$time_var]]))) +
                   scale_y_continuous(labels = function(x) {sprintf("%.2f", x)}) +
                   facet_wrap("Taxa", scales = "free_y", ncol = 5) + 
                   theme_bw() + 
                   theme(panel.grid = element_blank(), 
                         strip.text = element_text(size = 5, face = "italic"), 
                         legend.position = "bottom", 
                         axis.title = element_blank()) + 
                   scale_color_manual(values = PlotsAttr$color_gr, 
                                      name = "Group(lm)") +
                   scale_fill_manual(values = PlotsAttr$color_gr, 
                                       name = "Group(lm)") + 
                   scale_linetype_discrete(name = "Group(mean)") + 
                   scale_shape_discrete(name = "Group(mean)")
    
    ggsave(filename = paste0(iDirOut, "/", 
                             iSampleSet, "--", iCountNorm, ".png"), 
           plot = FullGraph, 
           width = 10, 
           height = (ceiling(length(SigTaxa$Variable)/5)* 1.25 + 1)) 
            
    
    IndPlotsLs <- list()

    # Individual Plots
    for(tax in SigTaxa$Variable) {

      IndPlotsLs[[tax]] <- TaxaLongDf %>%
                            filter(Taxa == tax) %>%
                            ggplot(aes(x = .data[[PRM$DA$time_var]])) +
                            geom_smooth(aes(x = .data[[PRM$DA$time_var]], 
                                            y = Taxa_Abundance, 
                                            fill = .data[[PRM$DA$group_var]], 
                                            color = .data[[PRM$DA$group_var]]), 
                                        method = "lm") +
                              geom_line(aes(y = Mean,
                                            group = .data[[PRM$DA$group_var]],
                                            linetype = .data[[PRM$DA$group_var]])) +
                              geom_point(aes(y = Mean,
                                             shape = .data[[PRM$DA$group_var]])) +
                              scale_x_continuous(breaks = sort(unique(TaxaLongDf[[PRM$DA$time_var]]))) +
                              scale_y_continuous(labels = function(x) {sprintf("%.2f", x)}) +
                              facet_wrap("Taxa", scales = "free_y", ncol = 5) +
                              theme_bw() +
                              theme(panel.grid = element_blank(),
                                    strip.text = element_text(size = 5, face = "italic"),
                                    legend.position = "bottom",
                                    axis.title = element_blank()) + 
                              scale_color_manual(values = PlotsAttr$color_gr, 
                                                 name = "Group(lm)") +
                              scale_fill_manual(values = PlotsAttr$color_gr, 
                                                name = "Group(lm)") + 
                              scale_linetype_discrete(name = "Group(mean)") + 
                              scale_shape_discrete(name = "Group(mean)")
      
      ggsave(filename = paste0(iDirOut, "/", iSampleSet, "--", iCountNorm, "/", 
                               tax, ".svg"), 
             plot = IndPlotsLs[[tax]], 
             width = 2.5, height = 2.25)

    }
    
    SpaghPlotsLs[[iTaxaLvl]][[iSampleSet]][[iCountNorm]] <- IndPlotsLs
    
    
    #-------------------------------------------------------------------------
    # Heat maps 
    #-------------------------------------------------------------------------
    ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
    
    HeatDf <- DataComb[[iSampleSet]][["ps"]][[iTaxaLvl]][[iCountNorm]] %>% 
                        prune_taxa(SigTaxa$Variable, .) %>% 
                        otu_table() %>% 
                        as.matrix() %>% 
                        as.data.frame()
    
    # Heat map plot size
    HeatHight <- nrow(HeatDf)/8
    
    HeatWidth <- HeatHight * 1.5
    
    HeatPlot <- Heatmap(HeatDf, 
                        name = paste0("Abundance (", iCountNorm, ")"), 
                        width = unit(HeatWidth, "in"), 
                        height = unit(HeatHight, "in"), 
                        show_column_names = FALSE, 
                        show_row_dend = FALSE, 
                        show_column_dend = FALSE,
                        column_split = Meta$GroupTime,
                        border_gp = gpar(col = "black"),
                        cluster_column_slices = FALSE, 
                        column_title_gp = gpar(fill = c("white")),
                        row_names_gp = gpar(fontsize = 6.5, 
                                            fontface = "italic"), 
                        heatmap_legend_param = list(
                          legend_direction = "horizontal", 
                          legend_width = unit(3.5, "cm")))
    
    
    png(filename = paste0(iDirOut, "/", 
                          iSampleSet, "--", iCountNorm, "--heat.png"), 
        width = HeatWidth + 2.5, 
        height = HeatHight + 1, units = "in", res = 300)
    draw(HeatPlot, heatmap_legend_side="bottom")
    dev.off()
  }
  
  
  #-----------------------------------------------------------------------------
  # Sig taxa used in publication (SCFA & Methane producers)
  #-----------------------------------------------------------------------------
  if(iTaxaLvl %in% c("Genus", "Family")) {
    
    PubSigTax <- PubSigTaxLs[[iTaxaLvl]]
    
    PresentTax <- intersect(PubSigTax, names(IndPlotsLs))
    
    if(length(PresentTax) > 0) {
      
      PlotsList <- lapply(IndPlotsLs[PresentTax], 
                          function(x){ x + theme(legend.position = "none")}) 
      
      RotLegend <- IndPlotsLs[[1]] + 
                      theme(legend.position = "right")
      
      PlotLegend <- get_legend(RotLegend)
    
      PlotGrid0 <- plot_grid(plotlist = PlotsList, ncol = 3, labels = "AUTO")
      
      FullGrid <- plot_grid(PlotGrid0, NULL,
                            PlotLegend, 
                            ncol = 3, 
                            rel_widths = c(1, 0.02, 0.15))
      
      
      save_plot(filename = paste0(iDirOut, "/", 
                               iSampleSet, "--", iCountNorm, "--SCFA.png"), 
             plot = FullGrid, 
             base_width = 7.5, 
             base_height = (ceiling(length(PresentTax)/3)* 1.5 + 1)) 
      
    }

  }
 
}


#-------------------------------------------------------------------------------
# Write results 
#-------------------------------------------------------------------------------
DaRes_03 <- list("da_res_ls" = DaRes, 
                 "da_res_df" = DaResDf, 
                 "da_sig_tax" = SigTaxaLs)

save(list = c("DaRes_03"), 
     file = paste0(PRM$data$out_dir, "/3_DA.Rdata"))

rm(list = ls())
gc()


