################################################################################
# Data preparation 
################################################################################

# Load parameters list 
load("PRM.Rdata")

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE, )}
rm(list = c("i"))

# Set seed 
set.seed(PRM$general$seed)

# Custom functions 
source("R/phy_shorten_tax_names.R")
source("R/phy_norm_count.R")

# Output directory 
DirOut <- PRM$data$out_dir

dir.create(DirOut, recursive = TRUE, showWarnings = FALSE)

# List for objects
DataComb <- list()

#-------------------------------------------------------------------------------
# Read data/ Create phyloseq
#-------------------------------------------------------------------------------
Ps <- qza_to_phyloseq(features = paste0(PRM$data$qiime_path, "asv_table.qza"), 
                       tree = paste0(PRM$data$qiime_path, "tree/rooted-tree.qza"), 
                       taxonomy = paste0(PRM$data$qiime_path, "taxonomy_07.qza"))

# Samples metadata
Meta <- read.csv(PRM$data$meta_path) %>% 
                  mutate(across(everything(), \(x) trimWhiteSpace(x))) %>% 
                  mutate(rownamecol = SeqID) %>% 
                  column_to_rownames("rownamecol") %>% 
                  select(-X.1) %>% 
                  mutate(SeqID = rownames(.))


# Samples to keep
MetaFilt <- Meta %>% 
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
                         RowNAMES = SeqID, 
                         Body_weight = as.numeric(Body.weight), 
                         Fasting_glucose = as.numeric(Fasting.glucose..mmol.L.), 
                         HbA1c = as.numeric(HbA1c))  %>% 
                  mutate(GroupTime = interaction(Group, Time)) %>% 
                  column_to_rownames("RowNAMES")


# Trim out samples 
Ps <-  prune_samples(rownames(MetaFilt), Ps) 

MetaFilt <- MetaFilt[sample_names(Ps), ]

sample_data(Ps) <- MetaFilt

DataComb[["all"]][["meta"]] <- MetaFilt


#-------------------------------------------------------------------------------
# Filter out taxa
#-------------------------------------------------------------------------------
#===============================================================================
# Revision 1 - 20 May 2025
PsAlpha <- Ps 

PsAlpha <- prune_taxa(!tax_table(PsAlpha)[, "Genus"] %in% "Mitochondria", PsAlpha)

PsAlpha <- prune_taxa(!tax_table(PsAlpha)[, "Genus"] %in% "Chloroplast", PsAlpha)

PsAlpha <- prune_taxa(!is.na(tax_table(PsAlpha)[, "Phylum"])[, "Phylum"], PsAlpha)

PsAlpha <- prune_taxa(tax_table(PsAlpha)[, "Kingdom"] %in% c("d__Bacteria", 
                                                             "d__Archaea"), 
                 PsAlpha)

# Add phyloseq to data list 
DataComb[["all"]][["ps"]][["ASV_alpha"]][["count"]] <- PsAlpha
#===============================================================================


# filter taxa with with less than X reads in total   
Ps <- prune_taxa(taxa_sums(Ps) >= PRM$data$min_reads_per_taxa, Ps)

# Remove ASVs: 
# Kingdom: "d__Eukaryota", "Unassigned"
# Genus: "Mitochondria"
Ps <- prune_taxa(!tax_table(Ps)[, "Genus"] %in% "Mitochondria", Ps)

Ps <- prune_taxa(!tax_table(Ps)[, "Genus"] %in% "Chloroplast", Ps)

Ps <- prune_taxa(!is.na(tax_table(Ps)[, "Phylum"])[, "Phylum"], Ps)

Ps <- prune_taxa(tax_table(Ps)[, "Kingdom"] %in% c("d__Bacteria", "d__Archaea"), 
                  Ps)

# Add phyloseq to data list 
DataComb[["all"]][["ps"]][["ASV"]][["count"]] <- Ps


#-------------------------------------------------------------------------------
# Glom to higher taxonomic level
#-------------------------------------------------------------------------------
for(i in PRM$data$glom_lvls)  {
  
  DataComb[["all"]][["ps"]][[i]][["count"]] <- tax_glom(Ps, i, 
                                                        NArm = PRM$data$glom_NArm)

}


#-------------------------------------------------------------------------------
# Adjust taxa names 
#-------------------------------------------------------------------------------
for(i in names(DataComb[["all"]][["ps"]])) {
  
  InstNames <- DataComb[["all"]][["ps"]][[i]][["count"]] %>% 
                  phy_shorten_tax_names(.) %>% 
                  make.unique()

  taxa_names(DataComb[["all"]][["ps"]][[i]][["count"]]) <- InstNames
  
}


#===============================================================================
# Revision 1
# Add PICRUSt2 predicted pathways into data set 
# Read in data, format and filter 
PathwayTab <- read_tsv(PRM$data$picrust_path, show_col_types = FALSE) %>% 
                as.data.frame() %>% 
                column_to_rownames(var = "pathway") %>% 
                select(all_of(sample_names(Ps)))

PathwayMockTax <- data.frame(Type = "Pathway", 
                             Pathway = rownames(PathwayTab), 
                             Rownames = rownames(PathwayTab)) %>% 
                    column_to_rownames(var = "Rownames") %>% 
                    as.matrix()

# Trim samples 
PsPathway <- phyloseq(otu_table = otu_table(PathwayTab, taxa_are_rows = TRUE), 
                      sample_data = sample_data(Ps), 
                      tax_table = tax_table(PathwayMockTax))

DataComb[["all"]][["ps"]][["pathways"]][["count"]] <- PsPathway
#===============================================================================

#-------------------------------------------------------------------------------
# Create filtered subsets
#-------------------------------------------------------------------------------
for(i in names(PRM$data$sample_subs)) { 
  
  InstMeta <- MetaFilt %>% 
                filter(CID %in% PRM$data$sample_subs[[i]]) %>% 
                filter(n() == length(PRM$data$sample_subs[[i]]), 
                       .by = Subject) %>% 
                droplevels()
  
  for(j in names(DataComb[["all"]][["ps"]])) { 
    
    DataComb[[i]][["ps"]][[j]][["count"]] <- 
                      prune_samples(rownames(InstMeta), 
                                    DataComb[["all"]][["ps"]][[j]][["count"]]) 
    
  }
  
  DataComb[[i]][["meta"]] <- 
                InstMeta[sample_names(DataComb[[i]][["ps"]][[j]][["count"]]), ]
  
}


#-------------------------------------------------------------------------------
# Normalize count   
#-------------------------------------------------------------------------------
NormGrid <- expand.grid("Set" = names(DataComb), 
                        "TaxLvl" = names(DataComb[[1]][["ps"]]), 
                        "Method" = PRM$data$count_norm, 
                        stringsAsFactors = FALSE)

for(i in 1:nrow(NormGrid)) {

  iSet <- NormGrid[i, "Set"]
  
  iTaxLvl <- NormGrid[i, "TaxLvl"]
  
  iMethod <- NormGrid[i, "Method"]
  
  DataComb[[iSet]][["ps"]][[iTaxLvl]][[iMethod]] <- 
        DataComb[[iSet]][["ps"]][[iTaxLvl]][["count"]] %>% 
        phy_norm_count(norm_type = iMethod, 
                       seed = PRM$general$seed, 
                       rare_depth = PRM$data$rare_depth)
                            
}


#-------------------------------------------------------------------------------
# Aesthetics for plots 
#-------------------------------------------------------------------------------
PlotsAttr <- list(color_gr = setNames(c("#377EB8", "red4"), 
                                    unique(MetaFilt$Group)), 
                  shape_cid = setNames(c(15:18), 
                                       na.omit(unique(MetaFilt$CID))), 
                  color_country = setNames(brewer.pal(4, "Set3"), 
                                           unique(MetaFilt$Country)))


#-------------------------------------------------------------------------------
# Write objects and clean environment 
#-------------------------------------------------------------------------------
save(list = c("DataComb", "PlotsAttr"), 
     file = paste0(DirOut, "/0_data.Rdata"))

rm(list = ls())
gc()
