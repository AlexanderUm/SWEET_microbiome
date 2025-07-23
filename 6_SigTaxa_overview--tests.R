################################################################################
# Investigate significant taxa 
################################################################################

load("PRM.Rdata")

load(paste0(PRM$data$out_dir, "/0_data.Rdata"))
load(paste0(PRM$data$out_dir, "/3_DA.Rdata"))
load(paste0(PRM$data$out_dir, "/4_Resp.Rdata"))
load(paste0(PRM$data$out_dir, "/4_RespRF.Rdata"))

source("R/RF_supp_functions.R")

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE)}

# Custom functions
source("R/phy_taxa_filter.R")

# Create directory 
DirOut <- PRM$overview$out_dir

dir.create(paste0(DirOut, "/plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(DirOut, "/tabs"), recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Overlap between significant taxa for DA and RF 
#-------------------------------------------------------------------------------
DaGrid <- expand.grid("da_lvl" = PRM$overview$da_taxa_lvl, 
                      "da_set" = PRM$overview$da_data_set, 
                      stringsAsFactors = FALSE)

RfColor <- RColorBrewer::brewer.pal(7, "Accent")

VennDataLs <- list()

OverTaxaDf <- data.frame(feature = NA)

for(i in 1:nrow(DaGrid)) {
  
  iDaLvl <- DaGrid[i, "da_lvl"]
  
  iDaSet <- DaGrid[i, "da_set"]
  
  DaSigTaxa <- DaRes_03$da_sig_tax[[iDaLvl]][[iDaSet]]$count
  
  VennDataLs[["DA"]] <- DaSigTaxa$Variable
  
  OverTaxaDf <- DaSigTaxa %>% 
                  mutate(count = 1) %>% 
                  select(Variable, count) %>% 
                  setNames(c("feature", paste0("DA_", iDaLvl, "_", iDaSet))) %>% 
                  full_join(OverTaxaDf, ., by = "feature")
  
  # Extract per RF model (Check that taxonomy levels are corresponding)
  for(j in 1:length(ResRFRes_04$RfPermuteSigRes)) {
    
    RfSigTaxa <- ResRFRes_04$RfPermuteSigRes[[j]]
    
    VennDataLs[[paste0("RF ", j)]] <- RfSigTaxa$feature
    
    OverTaxaDf <- RfSigTaxa %>% 
                    mutate(count = 1) %>% 
                    select(feature, count) %>% 
                    setNames(c("feature", 
                               paste0(paste(ResRFRes_04$ImpSetsInfo[[j]], collapse = "--"), 
                                 "(RF_", j, ")"))) %>% 
                    full_join(OverTaxaDf, ., by = "feature")
    
    # Save paired Venn diagram
    VennInst <- list(DaSigTaxa$Variable, RfSigTaxa$feature) %>% 
                      setNames(c("DA", paste0("RF ", j))) %>% 
                      ggvenn(fill_color = c("steelblue", RfColor[j]), 
                             set_name_size = 5) 
    
    ggsave(paste0(DirOut, "/plots/venn--DA_", 
                    paste0("RF", j), ".png"), 
           plot = VennInst, height = 5, width = 5)
    
  }
  
  # Add full taxonomy 
  Ps <- DataComb[[iDaSet]][["ps"]][[iDaLvl]][[1]]
  
  TaxTabInst <- Ps %>% 
                tax_table() %>% 
                as.matrix() %>% 
                as.data.frame()
  
  OverTaxaDf <- OverTaxaDf %>% 
                  mutate(Ntimes = rowSums(select(., -1), na.rm = TRUE))
  
  TaxTabInst[, colSums(!is.na(TaxTabInst)) > 0 ] %>% 
            rownames_to_column(var = "feature") %>% 
            right_join(OverTaxaDf, by = "feature") %>% 
            arrange(desc(Ntimes)) %>% 
    write.csv(paste0(DirOut, "/tabs/overlapping_taxa.csv"), na = "")
  
}

# All Random Forest 
VennDataLs[c(grep("RF", names(VennDataLs)))] %>% 
      ggvenn(., fill_color = RfColor[1:length(.)], 
             set_name_size = 5) %>% 
  ggsave(paste0(DirOut, 
                "/plots/venn--RF_combined.png"), 
         plot = ., height = 5, width = 5)


# All combined 
ggvenn(VennDataLs, 
       fill_color = c("steelblue", RfColor[1:length(VennDataLs)]), 
       set_name_size = 5) %>% 
  ggsave(paste0(DirOut, 
                "/plots/venn--All_combined.png"), 
         plot = ., height = 5, width = 5)
