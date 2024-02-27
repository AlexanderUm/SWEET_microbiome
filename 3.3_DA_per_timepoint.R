#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse",  "metagenomeSeq", "knitr", 
               "MicrobiomeStat", "Maaslin2")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")


#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------

tax.lvls <- c("Genus", "Family")

prev.cutoff.tp <- 0.5

i.ps <- "all"

test.var <- as.character(var.use["Group"])

adj.var.tp <- as.character(var.use["Country"])

ref.mas <- "Sugar"

norm.maas <- "TSS"


################################################################################
# Cross-sectional analysis per time point
#-------------------------------------------------------------------------------

da.per.tp <- list()

linda.df <- NULL  

maas.df <- NULL


# Subsets preparation

for(i.lvl in tax.lvls) {
    
    ps.inst <- pss.ls[[i.ps]][[i.lvl]]
    
    for (i.tp in unique(meta.ls[[i.ps]][[var.use[["CID"]]]]))  {
      
      meta.inst.tp <- meta.ls[[i.ps]] %>%
                      filter(.data[[var.use["CID"]]] == i.tp)
      
      ps.inst.tp <- prune_samples(rownames(meta.inst.tp), ps.inst)

      # LinDA cross-sectional method
      mstat <- mStat_convert_phyloseq_to_data_obj(ps.inst.tp)
      
      # Since LinDA convert everything to character I will replace metadata 
      # within the object with original metadata 
      mstat$meta.dat <- meta.inst.tp
      
      #-------------------------------------------------------------------------
      # Run cross-sectional LinDA
      LD.res <- generate_taxa_test_single(
                        data.obj = mstat,
                        group.var = test.var,
                        adj.vars = adj.var.tp,
                        feature.dat.type = "count",
                        feature.level = "original",
                        prev.filter = prev.da.cut.off,
                        abund.filter = 0)
      
      # Results to a data frame
      linda.df <- LD.res$original[[1]] %>% 
                          mutate(Tax_level = i.lvl, 
                                 Subset = i.tp, 
                                 Refference = "Sug",
                                 Method ="LD_TP", 
                                 Prevalence_CutOff = prev.da.cut.off) %>% 
                          bind_rows(linda.df, .)
      
      da.per.tp[["LD_TP"]][[i.lvl]][[i.tp]] <- LD.res$original[[1]]
      
      #-------------------------------------------------------------------------
      # Run cross-sectional MaasLin2
      
      # Output directory
      out.dir <- paste0("out/DA/maaslin/all/", i.lvl, "/", i.tp)
      
      dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
      
      # OTU table 
      otu.inst.tp <- ps.inst.tp %>%
                        otu_table() %>%
                        as.matrix() %>%
                        as.data.frame()
      
      maas.out <- Maaslin2(
                        input_data =  otu.inst.tp,
                        input_metadata = meta.inst.tp,
                        output = out.dir,
                        random_effects = adj.var.tp,
                        fixed_effects = test.var,
                        reference = ref.mas,
                        normalization = norm.maas,
                        transform = "LOG",
                        analysis_method = "LM",
                        min_abundance = 0,
                        min_prevalence = prev.da.cut.off,
                        min_variance = 0,
                        max_significance = 0.2,
                        correction = "BH",
                        standardize = TRUE,
                        cores = 4,
                        plot_heatmap = FALSE,
                        plot_scatter = FALSE)
      
      da.per.tp[["Maas_TP"]][[i.lvl]][[i.tp]] <- maas.out
      
      maas.df <- maas.out$results %>% 
                        mutate(Tax_level = i.lvl, 
                               Subset = i.tp, 
                               Method ="Maas_TP", 
                               Refference = "Sug", 
                               Prevalence_CutOff = prev.da.cut.off) %>% 
                        bind_rows(maas.df, .) 
      
      
}}
    

da.per.tp.df <- list(LinDA = linda.df, 
                     MaasLin = maas.df)  


#-------------------------------------------------------------------------------
# Write results 
#-------------------------------------------------------------------------------
save(list = c("prev.cutoff.tp", "da.per.tp.df", "da.per.tp", "adj.var.tp"), 
     file = "out/supp/da_res3.3_per_tp.Rdata")


