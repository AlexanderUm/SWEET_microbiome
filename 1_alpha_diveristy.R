#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
load("PRM.Rdata")

load(paste0(PRM$data$out_dir, "/0_data.Rdata"))

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE)}

# Create directory 
DirOut <- PRM$alpha$out_dir

dir.create(paste0(DirOut, "/plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(DirOut, "/tabs"), recursive = TRUE, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# Alpha diversity
#-------------------------------------------------------------------------------
AlphaRes <- list()

GridAlpha <- expand.grid("samples_set" = PRM$alpha$data_set,
                         "taxa_lvl" = PRM$alpha$taxa_lvl, 
                         "count_norm" = PRM$alpha$count_norm, 
                         stringsAsFactors = FALSE) 


for(i in 1:nrow(GridAlpha)) {
  
  iSampleSet <- GridAlpha[i, "samples_set"]
  
  iTaxaLvl <- GridAlpha[i, "taxa_lvl"]
  
  iCountNorm <- GridAlpha[i, "count_norm"]
  
  # Extract data
  Ps <- DataComb[[iSampleSet]]$ps[[iTaxaLvl]][[iCountNorm]]
  
  if(!taxa_are_rows(Ps)) {
    
    Ps <- phyloseq::t(Ps)
    
  }
    
  OtuTab <- Ps %>% 
              otu_table() %>% 
              as.matrix()
  
  
  #-----------------------------------------------------------------------------
  # Test differences
  #-----------------------------------------------------------------------------
  # Calculate alpha diversity with mStat
  AlphaObj <- mStat_calculate_alpha_diversity(OtuTab, 
                                               PRM$alpha$alpha_ind) 
    
  # Mstat object
  MstatObj <- mStat_convert_phyloseq_to_data_obj(Ps)
    
  # Change to class 
  MstatObj$meta.dat <- DataComb[[iSampleSet]][["meta"]]
    
  # Compare Trends
  AlphaTestRes <- generate_alpha_trend_test_long(
                                    data.obj = MstatObj, 
                                    alpha.obj = AlphaObj,
                                    alpha.name = PRM$alpha$alpha_ind,
                                    time.var = PRM$alpha$time_var,
                                    subject.var = PRM$alpha$subject_var,
                                    group.var = PRM$alpha$group_var,
                                    adj.vars = PRM$alpha$adjust_var)
    
  
  #-----------------------------------------------------------------------------
  # Plot Alpha diversity 
  #-----------------------------------------------------------------------------
  TimePoints <- MstatObj$meta.dat[[PRM$alpha$time_var]] %>% 
                  unique() %>% 
                  sort()
  
  AlphaDf <- as.data.frame(AlphaObj) %>% 
                      bind_cols(., MstatObj$meta.dat) %>% 
                      pivot_longer(cols = PRM$alpha$alpha_ind, 
                                   names_to = "Index", 
                                   values_to = "Value")
  
  AlphaDfMean <- AlphaDf %>% 
                    reframe(Mean = mean(Value), 
                              .by = c(PRM$alpha$time_var, 
                                      PRM$alpha$group_var, 
                                      "Index"))
  

  # Add statistics to plots 
  StatTab <- lapply(AlphaTestRes, 
                     function(x){ x <- x[nrow(x), ] 
                                 paste0(gsub("Group", "", x[[1]]),
                                        " [Est=", round(x[["Estimate"]], 2), 
                                        "; p = ", round(x[["P.Value"]], 3), "]")
                                 }) %>% 
                        unlist() %>% 
                        as.data.frame() %>% 
                        setNames("Text") %>% 
                        rownames_to_column(var = "Index")
            
    
  StatTabText <- AlphaDf %>% 
                summarise(y_text = max(Value)*1.1, 
                          .by = "Index") %>% 
                left_join(., StatTab, by = "Index") %>% 
                mutate(x_text = TimePoints[1])
    
    
    
  AlphaPlot <- ggplot() + 
                    geom_line(data = AlphaDf, 
                              aes(x = .data[[PRM$alpha$time_var]], 
                                  y = Value, 
                                  group = .data[[PRM$alpha$subject_var]],
                                  color = .data[[PRM$alpha$group_var]]), 
                                  alpha = 0.15) + 
                        geom_line(data = AlphaDfMean, 
                                  aes(x = .data[[PRM$alpha$time_var]], 
                                      y = Mean, 
                                      group = .data[[PRM$alpha$group_var]],
                                      color = .data[[PRM$alpha$group_var]]), 
                                  linewidth = 1) +
                        geom_point(data = AlphaDfMean, 
                                   aes(x = .data[[PRM$alpha$time_var]], 
                                       y = Mean), 
                                   size = 1) +
                        geom_text(data = StatTabText, 
                                  aes(label = Text,
                                      y = y_text, 
                                      x = x_text), 
                                  size =2.25, 
                                  hjust = 0) +
                        facet_wrap(~Index, scales = "free_y") + 
                        theme_bw() + 
                        scale_color_manual(values = PlotsAttr$color_gr) + 
                        scale_x_continuous(breaks = TimePoints) + 
                        xlab("Time(Month)") + 
                        ylab("Index value") + 
                        theme(panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank())
  
  
  #-----------------------------------------------------------------------------  
  # Collect and write results 
  #-----------------------------------------------------------------------------
  AlphaRes[[iSampleSet]][[iTaxaLvl]][[iCountNorm]][["Stat"]] <- AlphaTestRes
  
  AlphaRes[[iSampleSet]][[iTaxaLvl]][[iCountNorm]][["Plot"]] <- AlphaPlot
  
  Map(cbind, AlphaTestRes, Index=names(AlphaTestRes)) %>% 
    bind_rows() %>% 
    write.csv(., file = paste0(DirOut, "/tabs/AlphaLinDa_", 
                               iSampleSet, "_", iTaxaLvl, "_", 
                               iCountNorm, ".csv"))
  
  ggsave(filename = paste0(DirOut, "/plots/AlphaLinDa_", 
                           iSampleSet, "_", iTaxaLvl, "_", 
                           iCountNorm, ".svg"), 
         plot = AlphaPlot, 
         width = 7, height = 4)
  
    
}


#-------------------------------------------------------------------------------
# Save and clean up 
#-------------------------------------------------------------------------------
save(list = c("AlphaRes"), 
     file = paste0(PRM$data$out_dir, "/1_alpha.Rdata"))

rm(list = ls())
gc()
