#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
qiime.path <- "data/qiime2/"

meta.path <- "data/meta_comb_to_check_mss(24)_MP_25012024_MMS.csv"

min.reads.tax <- 50

prev.da.cut.off <- 0.5

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

# Change ASV IDs to shortened taxonomic names 
taxa_names(ps1) <- make.unique(phy_shorten_tax_names(ps1))

# Genus level 
ps1.genus <- tax_glom(ps1, "Genus")

taxa_names(ps1.genus) <- make.unique(phy_shorten_tax_names(ps1.genus))

# Family level 
ps1.family <- tax_glom(ps1, "Family")

taxa_names(ps1.family) <- make.unique(phy_shorten_tax_names(ps1.family))


pss.ls <- list(all = list(ASV = ps1, 
                          Genus = ps1.genus, 
                          Family = ps1.family))



#-------------------------------------------------------------------------------
# Filter and normalized
#-------------------------------------------------------------------------------
for(i.gr in c(0, 3, 4)) {
  
  for(i.ps in names(pss.ls$all)) {
    
    if(i.gr == 3) {
      
      meta <- filter(ps1.meta.f, CID != "CID_1") 
      
    } else {meta <- ps1.meta.f}
    
    
    if(i.gr == 0) {
      
      ps.f <- pss.ls$all[[i.ps]] 
      
      name.add <- "all"
      
    } else {
      
      name.add <- paste0("CIDs_", i.gr)   
      
      meta <- meta %>% 
                group_by(Subject) %>% 
                filter(n() == i.gr) %>% 
                ungroup() %>% 
                droplevels() %>% 
                mutate(RowNames = SeqID) %>% 
                column_to_rownames("RowNames") 
      
      ps.f <- prune_samples(rownames(meta), pss.ls$all[[i.ps]]) 
      
    }
    
    # Normalize count
    pss.ls[[name.add]][[i.ps]] <- ps.f
    
    ps.f.rare <- rarefy_even_depth(ps.f, rngseed = 30476048)
    
    pss.ls[[name.add]][[paste0("rare_", i.ps)]] <- ps.f.rare
    
    otu_table(ps.f) <- ps.f %>% 
      phyloseq_to_metagenomeSeq(.) %>% 
      cumNorm(., p=cumNormStatFast(.)) %>% 
      MRcounts(., norm=TRUE, log=TRUE) %>% 
      as.data.frame() %>% 
      otu_table(., taxa_are_rows = TRUE)
    
    pss.ls[[name.add]][[paste0("css_", i.ps)]] <- ps.f
    
  }
  
}


# Metadata 
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

#-------------------------------------------------------------------------------
# Aestetics for plots 
#-------------------------------------------------------------------------------
aest.ls <- list(color_gr = setNames(c("#377EB8", "red4"), 
                                    unique(ps1.meta.f[[var.use["Group"]]])), 
                shape_cid = setNames(c(15:18), 
                                     na.omit(unique(ps1.meta.f[[var.use["CID"]]]))), 
                color_country = setNames(brewer.pal(4, "Set3"), 
                                         unique(ps1.meta.f[[var.use["Country"]]])))


save(list = c("pss.ls", "meta.ls", "aest.ls", "var.use", 
              "meta.path", "min.reads.tax", "prev.da.cut.off"), 
     file = "out/supp/data_bundel.Rdata")
