#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
qiime.path <- "data/qiime2/"

meta.path <- "data/meta_comb_to_check_mss(24)_MP_25012024_MMS.csv"

min.reads.tax <- 50

prev.da.cut.off <- 0.5

glom.lvls <- c("Genus", "Family")

#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
set.seed(2395692)

dir.create("out/supp/", 
           recursive = TRUE, 
           showWarnings = FALSE)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "metagenomeSeq", "qiime2R")

for (i in libs.list) {library(i, character.only = TRUE, )}

rm(list = c("i", "libs.list"))

source("R/phy_shorten_tax_names.R")


#-------------------------------------------------------------------------------
# Read data/ Create phyloseq
#-------------------------------------------------------------------------------
ps1 <- qza_to_phyloseq(features = paste0(qiime.path, "asv_table.qza"), 
                       tree = paste0(qiime.path, "tree/rooted-tree.qza"), 
                       taxonomy = paste0(qiime.path, "taxonomy_07.qza"))

# Samples metadata
ps1.meta <- read.csv(meta.path) %>% 
                mutate(across(everything(), \(x) trimWhiteSpace(x))) %>% 
                mutate(rownamecol = SeqID) %>% 
                column_to_rownames("rownamecol") %>% 
                select(-X.1) %>% 
                mutate(SeqID = rownames(.))


# Samples to keep
ps1.meta.f <- ps1.meta %>% 
                  filter(!is.na(SampleID)) %>% 
                  filter(!is.na(CID)) %>% 
                  filter(group %in% c(1, 2)) %>% 
                  group_by(across("UniqPartID")) %>% 
                  filter(Nreads == max(Nreads)) %>% 
                  ungroup() %>% 
                  mutate(Time = as.numeric(case_match(CID, 
                                                      "CID_1" ~ 0, 
                                                      "CID_2" ~ 2, 
                                                      "CID_3" ~ 6, 
                                                      "CID_4" ~ 12)), 
                         Group = factor(case_match(as.character(group), 
                                                      "1" ~ "Sugar", 
                                                      "2" ~ "S&SEs"), 
                                        levels = c("Sugar", "S&SEs")), 
                         Country = as.factor(country), 
                         Subject = as.factor(X_subject_id), 
                         RowNAMES = SeqID)  %>% 
                  mutate(GroupTime = interaction(Group, Time)) %>% 
                  column_to_rownames("RowNAMES")



# Trim out samples 
ps1 <-  prune_samples(rownames(ps1.meta.f), ps1) 

ps1.meta.f <- ps1.meta.f[sample_names(ps1), ]

sample_data(ps1) <- ps1.meta.f


#-------------------------------------------------------------------------------
# Filter out taxa
#-------------------------------------------------------------------------------
# filter taxa with with less than X reads in total   
ps1 <- prune_taxa(taxa_sums(ps1) >= min.reads.tax, ps1)

# Remove ASVs: 
# Kingdom: "d__Eukaryota", "Unassigned"
# Genus: "Mitochondria"
ps1 <- prune_taxa(!tax_table(ps1)[, "Genus"] %in% "Mitochondria", ps1)

ps1 <- prune_taxa(!tax_table(ps1)[, "Genus"] %in% "Chloroplast", ps1)

ps1 <- prune_taxa(!is.na(tax_table(ps1)[, "Phylum"])[, "Phylum"], ps1)

ps1 <- prune_taxa(tax_table(ps1)[, "Kingdom"] %in% c("d__Bacteria", "d__Archaea"), 
                  ps1)


#-------------------------------------------------------------------------------
# Glom to higher taxonomic level
#-------------------------------------------------------------------------------
pss.ls <- list(all = list(ASV = ps1))

for(i.lvl in glom.lvls)  {
  
  pss.ls$all[[i.lvl]] <- tax_glom(ps1, i.lvl)

}


#-------------------------------------------------------------------------------
# Adjust taxa names 
#-------------------------------------------------------------------------------
for(i.lvl in names(pss.ls$all)) {

  taxa_names(pss.ls$all[[i.lvl]]) <- phy_shorten_tax_names(pss.ls$all[[i.lvl]]) %>% 
                                                make.unique()
}


#-------------------------------------------------------------------------------
# Create filtered subsets
#-------------------------------------------------------------------------------
gr.ls = list(CIDs_4 = c("CID_1", "CID_2", "CID_3", "CID_4"), 
             CIDs_3 = c("CID_2", "CID_3", "CID_4"))

for(gr in names(gr.ls)) {
  
  samp.id <- ps1.meta.f %>% 
                filter(CID %in% gr.ls[[gr]]) %>% 
                group_by(Subject) %>% 
                filter(n() == length(gr.ls[[gr]])) %>% 
                pull(SeqID)
  
  for(i.lvl in names(pss.ls$all)) {
    
   pss.ls[[gr]][[i.lvl]] <- prune_samples(samp.id, pss.ls$all[[i.lvl]]) 
    
  }
  
}


#-------------------------------------------------------------------------------
# Add ps objects with normalized counts  
#-------------------------------------------------------------------------------
for(i.set in names(pss.ls)) {
  
  for(i.lvl in names(pss.ls[[i.set]])) {
    
    ps.inst <- pss.ls[[i.set]][[i.lvl]]
    
    # Rarefaction 
    ps.inst.rare <- rarefy_even_depth(ps.inst, rngseed = 30476048)
    
    pss.ls[[i.set]][[paste0("rare_", i.lvl)]] <- ps.inst.rare
    
    # CSS normalization
    ps.inst.css <- ps.inst
    
    otu_table(ps.inst.css) <- ps.inst.css %>% 
                                phyloseq_to_metagenomeSeq(.) %>% 
                                cumNorm(., p=cumNormStatFast(.)) %>% 
                                MRcounts(., norm=TRUE, log=TRUE) %>% 
                                as.data.frame() %>% 
                                otu_table(., taxa_are_rows = TRUE)
    
    pss.ls[[i.set]][[paste0("css_", i.lvl)]] <- ps.inst.css
    
  }
}

#-------------------------------------------------------------------------------
# Extract Metadata per samples subset
#-------------------------------------------------------------------------------
meta.ls <- list()

for(i.ps in names(pss.ls)) {
  
  meta.ls[[i.ps]] <- pss.ls[[i.ps]][[1]] %>% 
    sample_data() %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    mutate(Time = as.numeric(trimws(Time)), 
           Group = factor(Group, levels = c("Sugar", "S&SEs")), 
           Plate = as.factor(Plate), 
           CID = as.factor(CID), 
           Country = as.factor(Country), 
           Subject = as.factor(Subject), 
           GroupCID = as.factor(paste0(Group, ":", gsub("CID_", "c", CID))), 
           CIDGroup = as.factor(paste0(gsub("CID_", "c", CID), ":", Group)))
  
}


var.use <- c(Time = "Time", Subject = "Subject", 
             Group = "Group", Country = "Country", 
             Plate = "Plate", CID = "CID") 


################################################################################
# Variables for following analysis 
################################################################################

#-------------------------------------------------------------------------------
# Aesthetics for plots 
#-------------------------------------------------------------------------------
aest.ls <- list(color_gr = setNames(c("#377EB8", "red4"), 
                                    unique(ps1.meta.f[[var.use["Group"]]])), 
                shape_cid = setNames(c(15:18), 
                                     na.omit(unique(ps1.meta.f[[var.use["CID"]]]))), 
                color_country = setNames(brewer.pal(4, "Set3"), 
                                         unique(ps1.meta.f[[var.use["Country"]]])))


#-------------------------------------------------------------------------------
# Variables for alpha diversity
#-------------------------------------------------------------------------------
alpha.vars.ls <- list(alpha_ind = c("shannon", "simpson", 
                                    "observed_species", "chao1"), 
                      used_ps = c("rare_ASV", "rare_Genus"), 
                      data_set = c("all", "CIDs_3", "CIDs_4"),
                      time_var = "Time", 
                      subject_var = "Subject", 
                      group_var = "Group", 
                      adjust_var = "Country")

#-------------------------------------------------------------------------------
# Variables for beta diversity
#-------------------------------------------------------------------------------
beta.vars.ls <- list(Distances = c(#"unifrac", "wunifrac", 
                                   "jaccard", "bray"),
                    used_ps = c("css_ASV", "css_Genus"), 
                    used_perm = 19,
                    test_var = "Group", 
                    full_RDA_f = "Time*Group + Condition(Country)",
                    full_data_set = c("all", "CIDs_3", "CIDs_4"),
                    strata_var = "CID",
                    strata_var_levels = c("CID_1", "CID_2", "CID_3", "CID_4"),
                    strata_RDA_f = "Group + Condition(Country)", 
                    strata_adonis_f = "Group",
                    strata_adonis_f_cov = "Country + Group",
                    strata_data_set = "all", 
                    p_color_var = "Group", 
                    p_shape_var = "CID", 
                    p_group_var = "Subject")


save(list = c("pss.ls", "meta.ls", "aest.ls", 
              "alpha.vars.ls", "beta.vars.ls"), 
     file = "out/supp/data_bundel.Rdata")
