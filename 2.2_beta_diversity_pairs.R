#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "ggsignif", "broom", "vegan", 
               "metagenomeSeq", "FSA", 
               "knitr",  "ggvegan", "ggpmisc", "cowplot")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")


################################################################################
# Beta diversity per time point 
################################################################################

#-------------------------------------------------------------------------------
# Variables grid for loops 
#-------------------------------------------------------------------------------
var.grid <- expand.grid(beta.vars.ls$strata_var_levels, 
                        beta.vars.ls$used_ps, 
                        beta.vars.ls$Distances, 
                        beta.vars.ls$strata_data_set) %>% 
              mutate(across(everything(), as.character))


rda.pairs <- list()

pcoa.pair.plots <- list()

beta.stat.summary <- NULL


for (i.grid in 1:nrow(var.grid)) {
  
  # Variables instances 
  i.tp <- var.grid[i.grid, 1]
  i.lvl <- var.grid[i.grid, 2]
  i.dist <- var.grid[i.grid, 3]
  i.set <- var.grid[i.grid, 4] 
  
  #-----------------------------------------------------------------------------
  # Exact metadata for time point
  meta.inst <- meta.ls[[i.set]] %>%
                filter(.data[[beta.vars.ls$strata_var]] == i.tp)
  
  #-----------------------------------------------------------------------------
  # Calculate distances
  ps.inst <- pss.ls[[i.set]][[i.lvl]] %>% 
                      prune_samples(rownames(meta.inst), .) 
  
  if(i.dist == "jaccard") {
    
    dist.inst <- phyloseq::distance(ps.inst, i.dist, binary = TRUE)
    
  } else {
    
    dist.inst <- phyloseq::distance(ps.inst, i.dist)
    
  }
  
  
  #-------------------------------------------------------------------------
  # Calculate betadisper
  beta <- betadisper(dist.inst,
                     group = meta.inst[[beta.vars.ls$test_var]])
  
  beta.disp <-  anova(beta)
  
  
  #-------------------------------------------------------------------------
  # Calculate RDA
  rda.form <- paste0("dist.inst ~ ", beta.vars.ls$strata_RDA_f) %>% 
                as.formula()
  
  rda.obj <- dbrda(rda.form, data = meta.inst)
  
  a.res <- anova.cca(rda.obj,
                     permutations = how(nperm = beta.vars.ls$used_perm)) 
  
  stat.tab <- a.res %>% 
                tidy(.) %>% 
                mutate(!!"R2(%)" := (round(SumOfSqs/sum(SumOfSqs)*100, 
                                           3))) %>%
                filter(term != "Residual") %>%
                select(c("term", "p.value", "R2(%)")) %>%
                suppressWarnings()
  
  rda.pairs[[i.lvl]][[i.dist]][[i.tp]] <- a.res
  
  
  #-----------------------------------------------------------------------------
  # Calculate ADONIS2
  # No covariates
  adon.form <- paste0("dist.inst ~ ", beta.vars.ls$strata_adonis_f) %>% 
                  as.formula()
  
  adon.res <- adonis2(adon.form, data = meta.inst, 
                      permutations = how(nperm = beta.vars.ls$used_perm)) %>% 
                  tidy(.) %>% 
                  suppressWarnings() %>%
                  mutate(Method = "ADONIS2: No covariates", 
                         formula = Reduce(paste, deparse(adon.form)), 
                         !!"R2(%)" := R2*100, 
                         GroupDispPval = beta.disp$`Pr(>F)`[1])
  
  # With covariates  
  adon.form.cov <- paste0("dist.inst ~ ", beta.vars.ls$strata_adonis_f_cov) %>% 
                  as.formula()
  
  adon.res.cov <- adonis2(adon.form.cov, data = meta.inst, 
                          permutations = how(nperm = beta.vars.ls$used_perm)) %>% 
                      tidy(.) %>% suppressWarnings() %>%
                      mutate(Method = "ADONIS2: With covariates", 
                             formula = Reduce(paste, deparse(adon.form.cov)), 
                             !!"R2(%)" := R2*100)
  
  
  ##############################################################################
  # Plot ordinations
  ##############################################################################
  info.tab <- data.frame(Category = c("Level", "Distance", "Time"),
                         Value = c(i.lvl, i.dist, i.tp))
  
  
  axis.pcoa <- beta$vectors[, 1:2] %>%
                    as.data.frame() %>%
                    bind_cols(., meta.inst)
  
  cent.pcoa <- beta$centroids[, 1:2] %>%
                    as.data.frame() %>%
                    setNames(c("center_x", "center_y")) %>%
                    mutate(!!beta.vars.ls$test_var := rownames(.))
  
  beta.plot.pcoa <- left_join(axis.pcoa, cent.pcoa, 
                              by = beta.vars.ls$test_var)
  
  # Make PCoA object
  a.pcoa <- ape::pcoa(dist.inst)
  
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
                    theme_bw() +
                    coord_equal() +
                    xlab(label = paste0("PCoA1 [", var.c[1] ,"%]")) +
                    ylab(label = paste0("PCoA1 [", var.c[2] ,"%]")) +
                    ggtitle(paste0("Formula: ", as.character(rda.form)[3]))
  
  pcoa.pair.plots[[i.lvl]][[i.dist]][[i.tp]] <- g.pcoa.p
  
  
  #-------------------------------------------------------------------------
  # Summary table statistics
  beta.stat.summary <- stat.tab %>% 
                          mutate(Method = "dbRDA", 
                                 formula = Reduce(paste, deparse(rda.form)), 
                                 term = beta.vars.ls$test_var) %>% 
                          slice_head(n = 1) %>% 
                          bind_rows(., adon.res, adon.res.cov) %>%
                          filter(term == "Group") %>%
                          select(-all_of(c("term", "df", 
                                           "SumOfSqs", "statistic"))) %>%
                          mutate(!!"Time Point" := i.tp, 
                                 !!"Taxa Level" := i.lvl,
                                 !!"Distance" := i.dist) %>% 
                          bind_rows(beta.stat.summary, .)
  
}

# Arrange table 
beta.stat.summary <- beta.stat.summary %>% 
                        arrange(.data[["Taxa Level"]], Distance, 
                                Method, .data[["Time Point"]])

# Arrange plots
pcoa.pair.plots.arr <- list()

for (i.lvl in names(pcoa.pair.plots)) {
    
  for (i.dist in names(pcoa.pair.plots[[i.lvl]])) {
    
    p.ls.inst <- pcoa.pair.plots[[i.lvl]][[i.dist]]
    
    legd.inst <- get_legend(p.ls.inst[[1]])
    
    p.ls.inst.f <- lapply(p.ls.inst, 
                          function(x) {x + theme(legend.position="none")})
    
    p.grid.1 <- plot_grid(plotlist = p.ls.inst.f, ncol = 2)
    
    p.grid.2 <- plot_grid(p.grid.1, legd.inst, ncol = 2, rel_widths = c(1, .1))
    
    
    pcoa.pair.plots.arr[[i.lvl]][[i.dist]] <- p.grid.2
    
}}
  
pair.forms <- c(rda = Reduce(paste, deparse(rda.form)), 
                adon.no.cov = Reduce(paste, deparse(adon.form)),
                adon.cov = Reduce(paste, deparse(adon.form.cov)))

save(list = c("beta.stat.summary", 
              "pcoa.pair.plots", 
              "rda.pairs", 
              "pcoa.pair.plots.arr", 
              "pair.forms"), 
     file = "out/supp/beta_pairs_res.Rdata")
