#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "ggsignif", "broom", "vegan", 
               "metagenomeSeq", "qiime2R", "picante", "phangorn", "FSA", 
               "knitr", "MicrobiomeStat", "ggvegan", "ggpmisc", "Maaslin2")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")

#-------------------------------------------------------------------------------
# Variables grid for loops 
#-------------------------------------------------------------------------------
var.grid <- expand.grid(beta.vars.ls$full_data_set, 
                        beta.vars.ls$used_ps, 
                        beta.vars.ls$Distances) %>% 
              mutate(across(everything(), as.character))

#-------------------------------------------------------------------------------
# Calculate distances 
#-------------------------------------------------------------------------------
comb.dist <- list()

for(i.grid in 1:nrow(var.grid)) {
  
   i.gr <- var.grid[i.grid, 1]
   i.lvl <- var.grid[i.grid, 2]
   i.dist <- var.grid[i.grid, 3]
    
   # Extract phyloseq 
   ps.inst <- pss.ls[[i.gr]][[i.lvl]]

   # Calculate distances
   if(i.dist == "jaccard") {

      comb.dist[[i.gr]][[i.lvl]][[i.dist]] <- distance(ps.inst,
                                                         i.dist,
                                                         binary = TRUE)
      } else {

      comb.dist[[i.gr]][[i.lvl]][[i.dist]] <- distance(ps.inst,
                                                         i.dist)

      }
}



################################################################################
# Full model
################################################################################

rda.comb.res.ls <- list()

rda.comb.plots.ls <- list()

pcoa.comb.plots.ls <- list()


for(i.grid in 1:nrow(var.grid)) {
  
  i.gr <- var.grid[i.grid, 1]
  i.lvl <- var.grid[i.grid, 2]
  i.dist <- var.grid[i.grid, 3]
  
  #-----------------------------------------------------------------------------
  # RDA
  #-----------------------------------------------------------------------------
  # Calculate RDA
  dist <- comb.dist[[i.gr]][[i.lvl]][[i.dist]]
  
  rda.form <- paste0("dist ~ ", beta.vars.ls$full_RDA_f) %>% 
    as.formula()
  
  rda.obj <- dbrda(rda.form, data = meta.ls[[i.gr]])
  
  a.res <- anova.cca(rda.obj,
                     permutations = how(nperm = beta.vars.ls$used_perm),
                     by = "terms")
  
  rda.comb.res.ls[[i.gr]][[i.lvl]][[i.dist]] <- a.res
  
  
  #-----------------------------------------------------------------------------
  # Extract data for plot
  #-----------------------------------------------------------------------------
  obj.sum <- summary(rda.obj)
  
  var.c <- round(obj.sum$cont$importance[2, 1:2]*100, 2)
  
  # Extract scores
  axis.score <- obj.sum$sites[, 1:2] %>%
                    as.data.frame() %>%
                    bind_cols(meta.ls[[i.gr]]) 
 
  # Extract centroids and vectors
  rda.p.auto <- autoplot(rda.obj) %>%
                    ggplot_build()
  
  centr.df <- rda.p.auto$data[[4]] %>%
                    mutate(Centroids = gsub(beta.vars.ls$test_var, "", label))
  
  vectors.df <- rda.p.auto$data[[3]] %>%
                    mutate(Group = gsub(beta.vars.ls$test_var, "", label))
  
  
  # Extract statistics
  stat.tab <- a.res %>% tidy() %>%
                    mutate(!!"R2(%)" := (round(SumOfSqs/sum(SumOfSqs)*100, 3))) %>%
                    filter(term != "Residual") %>%
                    select(c("term", "p.value", "R2(%)"))
  
  info.tab <- data.frame(Category = c("Subset", "Level", "Distance"),
                         Value = c(i.gr, i.lvl, i.dist))
  
  # Plot
  g.plot <- ggplot() +
                  geom_path(data = axis.score,
                            aes(x = dbRDA1, y = dbRDA2,
                                color = .data[[beta.vars.ls$p_color_var]],
                                group = .data[[beta.vars.ls$p_group_var]]),
                            alpha = 0.1) +
                  geom_point(data = axis.score,
                             aes(x = dbRDA1, y = dbRDA2,
                                 color = .data[[beta.vars.ls$p_color_var]],
                                 shape = .data[[beta.vars.ls$p_shape_var]]),
                             alpha = 0.5) +
                  geom_point(data = centr.df,
                             aes(x = x, y = y, fill = Centroids),
                             size = 2, shape = 23, color = "black") +
                  geom_segment(data = vectors.df,
                               aes(x = 0, y = 0, xend = x, yend = y),
                               arrow = arrow(length = unit(0.2, "cm"))) +
                  geom_text(data = vectors.df,
                            aes(x = x, y = y, label = Group), hjust = -0.1) +
                  scale_color_manual(values = aest.ls$color_gr) +
                  scale_fill_manual(values = aest.ls$color_gr) +
                  scale_shape_manual(values = aest.ls$shape_cid) +
                  coord_fixed() +
                  theme_bw() +
                  xlab(paste0(colnames(axis.score)[1], " [", var.c[1], "%]")) +
                  ylab(paste0(colnames(axis.score)[2], " [", var.c[2], "%]")) +
                  ggtitle(paste0("Formula: ", beta.vars.ls$full_RDA_f)) + 
                  guides(color = guide_legend(order=1),
                         fill = guide_legend(order=2),
                         shape = guide_legend(order=3))
                
                rda.comb.plots.ls[[i.gr]][[i.lvl]][[i.dist]] <- g.plot
                
                
  #-------------------------------------------------------------------------------
  # PCoA
  #-------------------------------------------------------------------------------
  beta <- betadisper(dist, group = meta.ls[[i.gr]][[beta.vars.ls$test_var]])
                
  axis.pcoa <- beta$vectors[, 1:2] %>%
                  as.data.frame() %>%
                  bind_cols(., meta.ls[[i.gr]])
                
  cent.pcoa <- beta$centroids[, 1:2] %>%
                  as.data.frame() %>%
                  setNames(c("center_x", "center_y")) %>%
                  mutate(!!beta.vars.ls$test_var := rownames(.))
                
  beta.plot.pcoa <- left_join(axis.pcoa, cent.pcoa, 
                              by = beta.vars.ls$test_var)
                
  # Make PCoA object
  a.pcoa <- ape::pcoa(dist)
                
  # Extract percentage
  var.c <- round((a.pcoa$values$Relative_eig[1:2]*100), 1)
                
                
  g.pcoa.p <- ggplot(beta.plot.pcoa) +
                  geom_segment(aes(x = center_x, y = center_y,
                                   xend = PCoA1, yend = PCoA2,
                                   color = .data[[beta.vars.ls$p_color_var]]),
                               alpha = 0.1) +
                  geom_point(aes(x = PCoA1, y = PCoA2,
                                 color = .data[[beta.vars.ls$p_color_var]],
                                 shape = .data[[beta.vars.ls$p_shape_var]]),
                             alpha = 0.5) +
                  geom_point(aes(x = center_x, y = center_y,
                                 fill = .data[[beta.vars.ls$p_color_var]]),
                             size = 2, shape = 23, color = "black") +
                  annotate(geom = 'table',
                           x=min(beta.plot.pcoa[, 1]),
                           y=max(beta.plot.pcoa[, 2])*1.5,
                           label=list(stat.tab), vjust = 0.85, hjust = 0.05) +
                  annotate(geom = 'table',
                           x=max(beta.plot.pcoa[, 1]),
                           y=max(beta.plot.pcoa[, 2])*1.5,
                           label=list(info.tab), vjust = 0.85, hjust = 0.9) +
                  scale_color_manual(values = aest.ls$color_gr) +
                  scale_fill_manual(values = aest.ls$color_gr) +
                  scale_shape_manual(values = aest.ls$shape_cid) +
                  theme_bw() +
                  coord_equal() +
                  xlab(label = paste0("PCoA1 [", var.c[1] ,"%]")) +
                  ylab(label = paste0("PCoA1 [", var.c[2] ,"%]")) +
                  ggtitle(paste0("Formula: ", beta.vars.ls$full_RDA_f))
                
  pcoa.comb.plots.ls[[i.gr]][[i.lvl]][[i.dist]] <- g.pcoa.p
  
}

save(list = c("rda.comb.res.ls", 
              "rda.comb.plots.ls", 
              "pcoa.comb.plots.ls"), 
     file = "out/supp/beta_res.Rdata")
