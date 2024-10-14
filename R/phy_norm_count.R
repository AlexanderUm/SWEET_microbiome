#-------------------------------------------------------------------------------
# Normalize ASV count with CSS
#-------------------------------------------------------------------------------
phy_norm_count <- function(phyloseq, 
                           norm_type = c("CSS", "CSS_log2", 
                                         "TSS", "TSS_log2", 
                                         "CLR", "Rare"), 
                           seed = 42, 
                           rare_depth = min(sample_sums(phyloseq))) {
  
  require("phyloseq")
  require("tidyverse")
  require("metagenomeSeq")
  require("vegan")
  
  TaxaAreRows <- taxa_are_rows(phyloseq)
  
  if(TaxaAreRows) { 
    
    MARGIN = 2
  
  } else { 
    
    MARGIN = 1}
  
  #-----------------------------------------------------------------------------
  # CSS
  #-----------------------------------------------------------------------------
  if(norm_type == "CSS") {
    
    otu_table(phyloseq) <- phyloseq %>% 
                            phyloseq_to_metagenomeSeq(.) %>% 
                            cumNorm(., p=cumNormStat(.)) %>% 
                            MRcounts(., norm=TRUE, log=FALSE) %>% 
                            as.data.frame() %>% 
                            otu_table(., taxa_are_rows = TaxaAreRows)
  }
      
   
  #-----------------------------------------------------------------------------
  # CSS Log2
  #-----------------------------------------------------------------------------
  if(norm_type == "CSS_log2") {
                               
     otu_table(phyloseq) <- phyloseq %>% 
                                 phyloseq_to_metagenomeSeq(.) %>% 
                                 cumNorm(., p=cumNormStat(.)) %>% 
                                 MRcounts(., norm=TRUE, log=TRUE) %>% 
                                 as.data.frame() %>% 
                                 otu_table(., taxa_are_rows = TaxaAreRows)
  }
  
  
  #-----------------------------------------------------------------------------
  # TSS 
  #-----------------------------------------------------------------------------
  if(norm_type == "TSS") { 
    
    otu_table(phyloseq) <-  phyloseq %>% 
                              otu_table() %>% 
                              as.matrix() %>% 
                              decostand(., MARGIN, method = "total") %>% 
                              as.data.frame() %>% 
                              otu_table(., taxa_are_rows = TaxaAreRows)
  }
  
  
  #-----------------------------------------------------------------------------
  # TSS Log2
  #-----------------------------------------------------------------------------
  if(norm_type == "TSS_log2") { 
    
    otu.tab <- phyloseq %>% 
                  otu_table() %>% 
                  as.matrix()
    
    OtuTab <- decostand(otu.tab, MARGIN, method = "total") %>% 
                as.data.frame() 
    
    OtuTab <-  log((OtuTab + 1)) 
                
    
    otu_table(phyloseq) <- otu_table(OtuTab, 
                                     taxa_are_rows = TaxaAreRows)
    
  }
  
  
  #-----------------------------------------------------------------------------
  # CLR 
  #-----------------------------------------------------------------------------
  if(norm_type == "CLR") {  
    
      otu_table(phyloseq) <- phyloseq %>% 
                                otu_table() %>% 
                                as.matrix() %>% 
                                decostand(., 
                                          MARGIN, 
                                          method = "clr", 
                                          pseudocount = 1) %>% 
                                as.data.frame() %>% 
                                otu_table(otu.tab.trans, 
                                          taxa_are_rows = TaxaAreRows)
    
  }
  
  #-----------------------------------------------------------------------------
  # Rarefy even depth
  #-----------------------------------------------------------------------------
  if(norm_type == "Rare") {
    
    phyloseq <- rarefy_even_depth(phyloseq, 
                                  sample.size = rare_depth, 
                                  rngseed = seed)
    
  }
  
  
  return(phyloseq)
  
}
