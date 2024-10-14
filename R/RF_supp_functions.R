################################################################################
# Extract features importance 
#-------------------------------------------------------------------------------
RF_extract_importance <- function(rfPermut_results, 
                                  col_arrange_by = "MeanDecreaseAccuracy_importance") {
  
  require(tidyverse)
  
  imp.df <- rfPermut_results$rf$importance %>% 
                as.data.frame() %>% 
                setNames(paste0(names(.), "_importance")) 
  
  imp.sd.df <- rfPermut_results$rf$importanceSD %>% 
                as.data.frame() %>% 
                setNames(paste0(names(.), "_SD")) 
              
  p.val.df <- rfPermut_results$pval %>% 
                as.data.frame() %>% 
                setNames(paste0(names(.), "_pval")) 
  
  comb.df <- bind_cols(imp.df, imp.sd.df, p.val.df) %>% 
                arrange(desc(.data[[col_arrange_by]])) %>% 
                rownames_to_column(var = "feature")
  
  return(comb.df) 
  
}


################################################################################
# Plot Features importance 
#-------------------------------------------------------------------------------
RF_plot_importance <- function(importance_tab, 
                               cols_to_plot = c("LIR_importance",
                                                "MIR_importance",
                                                "MeanDecreaseAccuracy_importance",
                                                "MeanDecreaseGini_importance"), 
                               col_arrange_by = "MeanDecreaseAccuracy_importance", 
                               pvals_to_use = "unscaled",
                               feature_col = "feature", 
                               remove_from_names = "_importance", 
                               extend_y_axis = 1.75, 
                               y_text_face_italic = TRUE) {
  
  require(tidyverse)
  
  pval.cols <- gsub(remove_from_names, "", cols_to_plot) %>% 
                paste0(., ".*", pvals_to_use) %>% 
                paste(., collapse = "|") %>% 
                grep(., colnames(importance_tab), value = TRUE)
  
  DfToPlot0 <- importance_tab %>% 
                  select(all_of(c(feature_col, cols_to_plot, pval.cols))) %>% 
                  setNames(gsub("\\.", "_", names(.))) %>% 
                  setNames(gsub(paste0("_", pvals_to_use), "", names(.))) %>% 
                  arrange(.data[[col_arrange_by]])
  
  DfToPlot <- DfToPlot0 %>% 
                  pivot_longer(-all_of(feature_col), 
                               names_to = c("Index", ".value"), 
                               names_sep = "_") %>% 
                  arrange(importance) %>% 
                  mutate(PvalShort = sprintf("%.3f", round(pval, 3)), 
                         Index = factor(Index, 
                                        levels = gsub(remove_from_names, "", 
                                                      cols_to_plot)), 
                         feature = factor(feature, levels = DfToPlot0$feature)) %>% 
                  mutate(PvalText = ifelse(PvalShort == 0, 
                                           "P<0.001", 
                                           paste0("P=", PvalShort)), 
                         !!"P-value" := ifelse(pval <= 0.05, pval, NA), 
                         PvalLabPos = ifelse(importance < 0, 0, importance))
  
  PlotOut <- ggplot(DfToPlot, aes(y = .data[[feature_col]], 
                                  x = importance, 
                                  fill = .data[["P-value"]])) + 
                geom_bar(stat = "identity") + 
                # geom_point(aes(x = importance*extend_y_axis), shape = NA) +
                # geom_text(aes(x = PvalLabPos, label = PvalText), 
                #           hjust = -0.1, size = 3) +
                facet_grid(~Index) + 
                theme_bw() + 
                theme(axis.title = element_blank(), 
                      axis.text.x = element_text(angle = 72.5, hjust = 1))
  
  
  if(y_text_face_italic) {
    
    PlotOut <-  PlotOut + 
                   theme(axis.text.y = element_text(face = "italic"))
    
  }
  
  return(PlotOut)
  
}



