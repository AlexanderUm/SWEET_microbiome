#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------
qvalue.cutoff <- 0.05

#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "ggsignif", "broom", "vegan", 
               "MicrobiomeStat", "ggvegan", "ggpmisc", "ComplexHeatmap", 
               "circlize", "RColorBrewer")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")
load("out/supp/da_res3.1_long.Rdata")
load("out/supp/da_res3.2_per_tp.Rdata")



methods.names <- c("Maas_TP" = "MaAsLin2 per time point",
                   "LD_tr" = "LinDA trend", 
                   "LD_TP" = "LinDA per time point")

#-------------------------------------------------------------------------------
# Combine data for plotting 
#-------------------------------------------------------------------------------
df.linda.f <- bind_rows(da.long.df$LinDA, 
                        da.per.tp.df$LinDA) %>% 
                        rename(p_value = P.Value, 
                               q_value = Adjusted.P.Value, 
                               Taxa = Variable) %>%
                        mutate(ColID = paste0(Method, "|set_", 
                                              gsub("CID_", "c", Subset), 
                                              "|ref_Sug|cont_S&SEs"))

# MaasLin results
df.maas.f <- da.per.tp.df$MaasLin %>% 
                      rename(Coefficient = coef, 
                             p_value = pval,
                             q_value = qval,
                             Taxa = feature)  %>% 
                  mutate(ColID = paste0(Method, 
                                        "|set_", gsub("CID_", "c", Subset), 
                                        "|ref_", Refference,
                                        "|cont_", value)) 


df.heat.sig.all <- bind_rows(df.maas.f, df.linda.f) %>% 
                    mutate(Subset = gsub("CID_", "c", Subset)) %>% 
                    filter(q_value <= qvalue.cutoff)


#-------------------------------------------------------------------------------
# Heat maps aesthetics
#-----------------------------------------------------------------------
# Heatmap fill 
col_fun = colorRamp2(c(-4, -0.0001, 0, 0.0001, 4), 
                     c("blue", "white", "black", "white", "gold1"))

# Color annotation - Method
col.method <- setNames(brewer.pal(length(unique(df.heat.sig.all$Method)), "Set1"), 
                       unique(df.heat.sig.all$Method))

col.samp.set <- c(setNames(brewer.pal(8, "YlGnBu")[c(8, 6, 4)], 
                           c("all", "CIDs_3", "CIDs_4")), 
                  setNames(brewer.pal(8, "BuGn")[c(8, 6, 4, 2)], 
                           c("c1", "c2", "c3", "c4")))

################################################################################
# Plot results per method 
#-------------------------------------------------------------------------------
# Extract data and plot heat maps 
heat.sum.ls <- list()

for (i.lvl in unique(df.heat.sig.all$Tax_level)) {
  
    filt.res.wide  <- df.heat.sig.all %>% 
                        filter(Tax_level == i.lvl) %>% 
                        select(Taxa, Coefficient, ColID) %>% 
                        pivot_wider(names_from = ColID, 
                                    values_from = Coefficient) %>% 
                        as.data.frame() %>% 
                        column_to_rownames(var = "Taxa") %>% 
                        relocate(order(colnames(.)))
         
    filt.res.wide[is.na(filt.res.wide)] <- 0

  
  
  
  # Top annotation 
  meth.v <- gsub("\\|.*", "", colnames(filt.res.wide))
  
  cid.v <- gsub(".*set_", "", colnames(filt.res.wide)) %>% 
              gsub("\\|.*", "", .)
            
  
  top.anot <- HeatmapAnnotation(Method = meth.v, 
                                Subset = cid.v, 
                                col = list(Method = col.method, 
                                           Subset = col.samp.set), 
                                gp = gpar(col = "black"))
  
  # Heat map itself
  heat.p <- Heatmap(filt.res.wide, 
                     name = paste0("Coef (", i.lvl, ")"), 
                     rect_gp = gpar(col = "gray40", lwd = 0.1), 
                     column_names_rot = 45, col = col_fun , 
                     cluster_columns = FALSE,
                     top_annotation = top.anot)
  
  # Write results to list
  heat.name <- paste0("Summary of taxa identified as significant by ", 
                      paste(methods.names, collapse = ", "), 
                      " at ", i.lvl, " level.")
  
  heat.ls <- list("plot" = heat.p, 
                  "ncol" = ncol(filt.res.wide), 
                  "nrow" = nrow(filt.res.wide), 
                  "name" = heat.name)
  
  # Heat map 
  heat.sum.ls[[i.lvl]] <- heat.ls
  
}


################################################################################
# Visualization per method
################################################################################

# OTU tables for plots 
otu.tabs.ls <- list()

for(i.lvl in unique(df.heat.sig.all$Tax_level)) {
  
  otu.inst <- pss.ls[["all"]][[i.lvl]] %>% 
                otu_table() %>%
                as.matrix() %>%
                as.data.frame()
  
  otu.inst.prop <- apply(otu.inst, 2, function(x){x/sum(x)})
  
  otu.tabs.ls[[i.lvl]][["clr"]] <- qiime2R::make_clr(otu.inst)
  
  otu.tabs.ls[[i.lvl]][["rel"]] <- otu.inst.prop*100
  
  otu.tabs.ls[[i.lvl]][["log"]] <-  pss.ls[["all"]][[paste0("rare_", i.lvl)]] %>% 
                    transform_sample_counts(., function(x){log((x+1))}) %>%            
                                                              otu_table() %>%
                                                              as.matrix() %>%
                                                              as.data.frame()
}

method.vis.ls <- list()

for(i.lvl in unique(df.heat.sig.all$Tax_level)) {
  
  for(i.method in unique(df.heat.sig.all$Method)) {
    
     sig.tax <- df.heat.sig.all %>% 
                    filter(Method == i.method, 
                           Tax_level == i.lvl) 
     
     sig.tax.uniq <- sig.tax %>% 
                       group_by(Taxa) %>%
                       summarise(Coef = sum(abs(Coefficient))) %>% 
                       arrange(desc(Coef)) %>%
                       pull(Taxa)
     
     #--------------------------------------------------------------------------
     # Heat map with other methods
     #--------------------------------------------------------------------------
     df.heat.inst <- df.heat.sig.all %>% 
                         filter(Tax_level == i.lvl, 
                                Taxa %in% sig.tax.uniq) %>% 
                         select(Taxa, Coefficient, ColID) %>% 
                         pivot_wider(names_from = ColID, 
                                     values_from = Coefficient) %>% 
                         as.data.frame() %>% 
                         column_to_rownames(var = "Taxa") %>% 
                         relocate(order(colnames(.))) %>% 
                         .[sig.tax.uniq, ]
     
     df.heat.inst[is.na(df.heat.inst)] <- 0
     
     # Top annotation 
     meth.v <- gsub("\\|.*", "", colnames(df.heat.inst))

     cid.v <- gsub(".*set_", "", colnames(df.heat.inst)) %>% 
                 gsub("\\|.*", "", .)
     
     
     top.anot <- HeatmapAnnotation(Method = meth.v, 
                                   Subset = cid.v, 
                                   col = list(Method = col.method, 
                                              Subset = col.samp.set), 
                                   gp = gpar(col = "black"))
     
     
     # Heat map 
     heat.sig.p <- Heatmap(df.heat.inst, 
                           name = paste0("Coef (", i.lvl, ")"), 
                           rect_gp = gpar(col = "gray40", lwd = 0.1), 
                           column_names_rot = 45, col = col_fun , 
                           cluster_columns = FALSE, 
                           cluster_rows = FALSE,
                           top_annotation = top.anot)
     
     
     heat.name <- paste0("Taxa identified as significant by ", 
                         methods.names[i.method], " at ", i.lvl, 
                         " level, across other methods.")
     
     heat.ls <- list("plot" = heat.sig.p, 
                     "ncol" = ncol(df.heat.inst), 
                     "nrow" = nrow(df.heat.inst), 
                     "name" = heat.name)
     
     method.vis.ls[[i.lvl]][[i.method]][["det_heat"]] <- heat.ls
     
     
     #--------------------------------------------------------------------------
     # Normalization
     for(i.norm in c("clr", "rel", "log")) {
       
      abund.name <- paste0("Abundance(", i.norm, ")")
       
      otu.tabs.f <- otu.tabs.ls[[i.lvl]][[i.norm]] %>% 
                      .[sig.tax.uniq, ] %>% 
                      t() %>% 
                      as.data.frame()
      
      
      #-------------------------------------------------------------------------
      # Spaghetti plots
      #-------------------------------------------------------------------------
      spag.df <- bind_cols(otu.tabs.f, 
                           meta.ls[["all"]]) %>% 
                pivot_longer(cols = all_of(unique(sig.tax$Taxa)), 
                             names_to = "Taxa", 
                             values_to = abund.name) %>% 
                mutate(Taxa = factor(Taxa, 
                        levels = sig.tax.uniq))
      
      # Mean and median data frame for spaghetti plots
      spag.df.s <- spag.df %>% 
                    group_by(across(c("Taxa", 
                                      da.vars.ls$time_var, 
                                      da.vars.ls$group_var))) %>% 
                    summarise(Median = median(.data[[abund.name]], 
                                              na.rm = TRUE), 
                              Mean = mean(.data[[abund.name]], 
                                          na.rm = TRUE), 
                              SD = sd(.data[[abund.name]], 
                                      na.rm = TRUE)) %>% 
                    mutate(!!"Groups(Median)" := Group, 
                           !!"Groups(Mean)" := Group, 
                           SD_up = Mean + SD,
                           SD_down = Mean - SD)

      # Base Plot 
      base.p <- ggplot() + 

        facet_wrap("Taxa", scales = "free", ncol = 4) + 
        theme_bw() + 
        scale_x_continuous(breaks = unique(spag.df[[da.vars.ls$time_var]])) + 
        scale_fill_manual(values = aest.ls$color_gr) + 
        scale_color_manual(values = aest.ls$color_country)  + 
        guides(
          # color = guide_legend(override.aes = list(linewidth = 2.5, alpha = 1), 
          #                      order = 4), 
               linetype = guide_legend(order = 2), 
               shape = guide_legend(order = 3), 
               fill = guide_legend(order = 1, title = "Groups(lm)"))
      
      # Volatility and trend plots
      if(i.method %in% c("LD_tr", "LD_vol")) {
  
        spag.p <- base.p + 
                    geom_smooth(data = spag.df, 
                                aes(x = .data[[da.vars.ls$time_var]], 
                                    y = .data[[abund.name]], 
                                    group = .data[[da.vars.ls$group_var]], 
                                    fill = .data[[da.vars.ls$group_var]]), 
                                color = "gray40",
                                linewidth = 1, method = "lm") +
                    geom_point(data = spag.df.s, 
                               aes(x = .data[[da.vars.ls$time_var]], 
                                   y = Mean, 
                                   shape = .data[["Groups(Mean)"]]), 
                               size = 1.5) + 
                    geom_line(data = spag.df.s, 
                              aes(x = .data[[da.vars.ls$time_var]], 
                                  y = Mean, 
                                  group = Group, 
                                  linetype = .data[["Groups(Mean)"]]), 
                              linewidth = 0.75)
        
      } else { 
        
        spag.p <- base.p + 
                    geom_line(data = spag.df.s, 
                              aes(x = .data[[da.vars.ls$time_var]], 
                                  y = Mean, 
                                  group = Group), 
                              linewidth = 0.5, 
                              color = "gray50") +
                    geom_point(data = spag.df.s, 
                               aes(x = .data[[da.vars.ls$time_var]], 
                                   y = Mean, 
                                   shape = .data[["Groups(Mean)"]]), 
                               size = 2) + 
                    geom_line(data = spag.df.s, 
                              aes(x = .data[[da.vars.ls$time_var]], 
                                  y = Mean, 
                                  group = .data[[da.vars.ls$group_var]], 
                                  linetype = .data[["Groups(Mean)"]]), 
                              linewidth = 0.75)
      }
      
      
      spag.name <- paste0("Taxa identified as significant by ", 
                          methods.names[i.method], " at ", i.lvl, 
                          " level. Abundance ", i.norm, " normalized.")
      
      spag.ls <- list("plot" = spag.p, 
                      "nrow" = ceiling(length(unique(spag.df$Taxa))/4), 
                      "name" = spag.name)
      
      method.vis.ls[[i.lvl]][[i.method]][["spagg"]][[i.norm]] <- spag.ls
      
      
      #-------------------------------------------------------------------------
      # Heat map composition
      #-------------------------------------------------------------------------
      otu.tabs.f.heat <- otu.tabs.f %>% 
                               t()
      
      #-------------------------------------------------------------------------
      # All samples 
      meta.heat <- meta.ls[["all"]]
      
      # Heat map 
      top.anot <- HeatmapAnnotation(Group = meta.heat[[da.vars.ls$group_var]], 
                                    col = list(Group = aest.ls$color_gr))
      
      heat.p1 <- Heatmap(otu.tabs.f.heat, 
                    name = paste0("Abundance (", i.norm, ")"), 
                    cluster_column_slices = FALSE,
                    cluster_columns = TRUE, 
                    cluster_rows = FALSE,
                    column_split = meta.heat[["CIDGroup"]],
                    show_column_names = FALSE, 
                    top_annotation = top.anot, 
                    border = TRUE)
      
      
      heat.name1 <- paste0("Taxa identified as significant by ", 
                          methods.names[i.method], " at ", i.lvl, 
                          " level. Abundance ", i.norm, " normalized.", 
                          " All sequenced samples.")
      
      heat.ls1 <- list("plot" = heat.p1, 
                      "ncol" = ncol(otu.tabs.f.heat), 
                      "nrow" = nrow(otu.tabs.f.heat), 
                      "name" = heat.name1)
      
      #-------------------------------------------------------------------------
      # Heat map CIDs_4
      meta.heat <- meta.ls[["CIDs_4"]] %>% 
                      arrange(CID, Group, Country, Subject)
    
      otu.heat <- otu.tabs.f.heat[, rownames(meta.heat)]

      top.anot <- HeatmapAnnotation(Group = meta.heat[[da.vars.ls$group_var]], 
                                    col = list(Group = aest.ls$color_gr))
      
      heat.p2 <- Heatmap(otu.heat, 
                        name = paste0("Abundance (", i.norm, ")"), 
                        cluster_column_slices = FALSE,
                        cluster_columns = FALSE, 
                        column_split = meta.heat[["CIDGroup"]], 
                        column_labels = meta.heat[[da.vars.ls$subj_var]],
                        show_column_names = FALSE, 
                        top_annotation = top.anot, 
                        border = TRUE)
      
      heat.name2 <- paste0("Taxa identified as significant by ", 
                          methods.names[i.method], " at ", i.lvl, 
                          " level. Abundance ", i.norm, " normalized.", 
                          " Samples with all time points.")
      
      heat.ls2 <- list("plot" = heat.p2, 
                      "ncol" = ncol(otu.heat), 
                      "nrow" = nrow(otu.heat), 
                      "name" = heat.name2)
      
      heat.ls.comb <- list(all = heat.ls1, 
                           CIDs4 = heat.ls2)

      method.vis.ls[[i.lvl]][[i.method]][["heat_comp"]][[i.norm]] <- heat.ls.comb
      
      }
     
  }
    
}
  

save(list = c("qvalue.cutoff", 
              "heat.sum.ls",
              "method.vis.ls"),
     file = "out/supp/da_vis.Rdata")


