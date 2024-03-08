#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "ggsignif", "broom", "vegan", 
               "metagenomeSeq", "qiime2R", "picante", "phangorn", "FSA", 
               "knitr", "MicrobiomeStat", "ggvegan", "ggpmisc", 
               "Maaslin2", "ComplexHeatmap", "circlize", "splinectomeR", )

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")


#-------------------------------------------------------------------------------
# Define responders and non-responders 
#-------------------------------------------------------------------------------
focus.var <- c("Body_weight", "HbA1c", "Fasting_glucose")

meta.var <- c("Subject", "CID", "Group", "Country")

sug.gr <- c("Sugar", "S&SEs")


resp.df.ls <- meta.ls$all %>% 
                group_by(Subject) %>% 
                reframe(across(focus.var, 
                               function(x){
                                  (x[CID == "CID_4"] - x[CID == "CID_2"])/
                                  (x[CID == "CID_1"] - x[CID == "CID_2"])})) %>% 
                mutate(across(focus.var, 
                              function(x){ifelse(x >= 0, "Resp", "NonResp")}, 
                              .names = "Resp_{.col}")) %>% 
                left_join(., meta.ls$all[, c("Subject", "Group")], 
                          by = "Subject", 
                          multiple = "first") %>% 
                split(., .[["Group"]])


#-------------------------------------------------------------------------------
# Use K-mean clustering to use all 3 variables
#-------------------------------------------------------------------------------
# Not numeric signs will be replaced with 0
library(factoextra)
library(cluster)

km.vars <- c("Body_weight", "Fasting_glucose", "HbA1c")

n.clust.ls <- list()
km.res.ls <- list()
kmean.df.ls <- list()

for(i.gr in names(resp.df.ls)) {
  
  kmean.df <- resp.df.ls[[i.gr]] %>% 
                  as.data.frame() %>% 
                  column_to_rownames("Subject") %>% 
                  select(all_of(km.vars)) %>% 
                  mutate(across(everything(), log)) %>% 
                  .[is.finite(rowSums(.)),]
  
  
  n.clust.ls[[i.gr]] <- fviz_nbclust(kmean.df, kmeans, method = "wss")
  
  km.res.ls[[i.gr]] <- kmeans(kmean.df, centers = 2, nstart = 2500)
  
  kmean.df.ls[[i.gr]] <- kmean.df
  
}

resp.df.ls$Sugar

fviz_cluster(km.res.ls$`S&SEs`, data = kmean.df.ls$`S&SEs`)


#-------------------------------------------------------------------------------
# Write out results 
#-------------------------------------------------------------------------------

resp.data.ls <- list(Resp_data = resp.df.ls)

save(list = c("resp.data.ls"), 
     file = "out/supp/5_Responders.Rdata")
