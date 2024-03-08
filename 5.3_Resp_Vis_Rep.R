#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(43957634)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "ggsignif", "broom", "vegan", 
               "metagenomeSeq", "qiime2R", "picante", "phangorn", "FSA", 
               "knitr", "MicrobiomeStat", "ggvegan", "ggpmisc", 
               "Maaslin2", "ComplexHeatmap", "circlize", "splinectomeR")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")
load("out/supp/5_Responders.Rdata")
load("out/supp/5.1_Resp_DA.Rdata")
load("out/supp/5.2_Resp_RF.Rdata")


#-------------------------------------------------------------------------------

resp.vis.var.ls <- list("Focus_var" = c("Body_weight", "HbA1c", 
                                        "Fasting_glucose"), 
                        "Meta_var" = c("Subject", "CID", "Group", "Country"), 
                        "Sub_Groups" = c("Sugar", "S&SEs"))

################################################################################
# General about responders 
################################################################################
resp.vis.gen <- list()

#-------------------------------------------------------------------------------
# Response variables overview 
#-------------------------------------------------------------------------------
meta.resp.long <- meta.ls$all %>% 
                    pivot_longer(all_of(resp.vis.var.ls[["Focus_var"]]))

resp.vis.gen[["Overview"]] <- ggplot(meta.resp.long) + 
                                    geom_line(aes(x = CID, 
                                                  y = value, 
                                                  group = Subject), 
                                              alpha = 0.75) + 
                                    facet_grid(name ~ Group, scales = "free") + 
                                    theme_bw()

#-------------------------------------------------------------------------------
# Distribution of calculated resp non resp
#-------------------------------------------------------------------------------
resp.df.long <- resp.df.ls %>% 
                  bind_rows() %>% 
                  pivot_longer(all_of(resp.vis.var.ls[["Focus_var"]])) 

resp.vis.gen[["Resp_dist"]] <- ggplot(resp.df.long, aes(x = Group, 
                                                        y = value)) + 
                                  geom_jitter(width = 0.2, 
                                              height = 0, 
                                              alpha = 0.25) +
                                  geom_violin(fill = NA) +
                                  facet_wrap(. ~ name, scales = "free") + 
                                  theme_bw()


#-------------------------------------------------------------------------------
# Summary tables
#-------------------------------------------------------------------------------
sum.tab.long <- resp.df.ls %>% 
                  bind_rows() %>% 
                  pivot_longer(starts_with("Resp_")) 

resp.vis.gen[["Resp_count"]]<- table(sum.tab.long[, c("name", 
                                                      "value", 
                                                      "Group")]) %>% 
                                  as.data.frame() %>% 
                                  pivot_wider(names_from = value, 
                                              values_from = Freq)


#-------------------------------------------------------------------------------
# Responders table
#-------------------------------------------------------------------------------
resp.vis.gen[["Resp_table"]] <- resp.df.ls


#-------------------------------------------------------------------------------
# Missing per variable
#-------------------------------------------------------------------------------
resp.df.long <- bind_rows(resp.df.ls)

missing.tab.ls <- list()

for(i in resp.vis.var.ls$Focus_var) {
  
  no.na.df <- resp.df.long %>% 
                  select(all_of(c(i, "Subject"))) %>% 
                  .[is.finite(.[[i]]),]
  
  missing.tab.ls[[i]] <- meta.ls$all %>% 
                          select(all_of(c(i, "Subject", "CID", 
                                          "Country", "Group"))) %>% 
                          filter(!.data[["Subject"]] %in% no.na.df[["Subject"]]) %>% 
                          pivot_wider(names_from = "CID", values_from = i)
  
}

resp.vis.gen[["Missing_table"]] <- missing.tab.ls

#-------------------------------------------------------------------------------
# Use K-mean clustering to use all 3 variables
#-------------------------------------------------------------------------------
# Not numeric signs will be replaced with 0
library(factoextra)
library(cluster)

km.vars.ls <- list("Two" = c("Body_weight", "Fasting_glucose"), 
                   "Three" = c("Body_weight", "Fasting_glucose", "HbA1c"))

km.res.ls <- list()

for(i.gr in names(resp.df.ls)) {
  for(i.vars in names(km.vars.ls)) {
    
     kmean.df <- resp.df.ls[[i.gr]] %>% 
                  as.data.frame() %>% 
                  column_to_rownames("Subject") %>% 
                  select(all_of(km.vars.ls[[i.vars]])) %>% 
                  mutate(across(everything(), log)) %>% 
                  .[is.finite(rowSums(.)), ]
     
     res.inst <- kmeans(kmean.df, centers = 2, nstart = 2500)
  
     clust.inst <- fviz_cluster(res.inst, data = kmean.df)
     
     km.res.ls[[i.gr]][[i.vars]] <- list("Data" = kmean.df, 
                                         "Results" = res.inst, 
                                         "Clust_plot" = clust.inst)
            
}}
  
resp.vis.gen[["K-mean"]] <- km.res.ls
  
################################################################################
# RF
################################################################################

#-------------------------------------------------------------------------------
# Classification accuracy
#-------------------------------------------------------------------------------

for(i in 1:length(resp.RF.res$RF_comb)) {
  
title.inst <- resp.RF.res$RF_comb[[i]][["Parameters"]] %>% 
                paste(., collapse = " | ") %>% 
                paste0("Parameters: ", .)
  
cons.out <- capture.output(resp.RF.res$RF_comb[[i]][["RF_res"]])
  
class.plot.df <- c(paste0("Responders", cons.out[11]), cons.out[12:14]) %>% 
                      str_squish() %>% 
                      str_split(pattern = " ", simplify = TRUE) %>% 
                      as.data.frame() %>% 
                      select(-c("V2", "V3")) %>% 
                      setNames(.[1, ]) %>% 
                      .[-1, ] %>% 
                      mutate(Responders = factor(Responders, 
                                                 levels = c("Overall", 
                                                            "Resp", 
                                                            "NonResp")), 
                             across(all_of(c("pct.correct", 
                                             "LCI_0.95", 
                                             "UCI_0.95")), as.numeric))


resp.RF.res$RF_comb[[i]][["Pred_plot"]] <- ggplot(class.plot.df) + 
                geom_hline(yintercept = 50, size = 1, 
                           color = "gray40") +
                geom_point(aes(x = Responders, 
                               y = pct.correct), size = 3) +
                geom_errorbar(aes(x = Responders, 
                                  ymin = LCI_0.95, 
                                  ymax = UCI_0.95), width=.1 ) + 
                theme_bw() + 
                ylab("Correctly Predicted (%)") +
                ggtitle(title.inst) + 
                theme(plot.title = element_text(size=8, 
                                                face="bold.italic"))

p.imp <- plotImportance(resp.RF.res$RF_comb[[i]][["RF_res"]], 
                        sig.only = TRUE, plot = FALSE)

resp.RF.res$RF_comb[[i]][["Pred_sig_var"]] <- p.imp
  
}

#-------------------------------------------------------------------------------
# Per time point 
#-------------------------------------------------------------------------------
RF.vis.tp <- list()

acc.long.tp <- NULL

for(i in 1:length(resp.RF.res$RF_per_tp)) {
  
  prm.inst <- resp.RF.res$RF_per_tp[[i]][["Parameters"]]
  
  title.inst <-  prm.inst %>% 
                  paste(., collapse = " | ") %>% 
                  paste0("Parameters: ", .)
  
  cons.out <- capture.output(resp.RF.res$RF_per_tp[[i]][["RF_res"]])
  
  acc.long.tp <- c(paste0("Responders", cons.out[11]), cons.out[12:14]) %>% 
                        str_squish() %>% 
                        str_split(pattern = " ", simplify = TRUE) %>% 
                        as.data.frame() %>% 
                        select(-c("V2", "V3")) %>% 
                        setNames(.[1, ]) %>% 
                        .[-1, ] %>% 
                        mutate(Responders = factor(Responders, 
                                                   levels = c("Overall", 
                                                              "Resp", 
                                                              "NonResp")), 
                               across(all_of(c("pct.correct", 
                                               "LCI_0.95", 
                                               "UCI_0.95")), as.numeric)) %>% 
                       cross_join(., data.frame(as.list(prm.inst))) %>%
                       bind_rows(acc.long.tp, .)
                       
}
  

save(list = c("resp.vis.gen"), file = "out/supp/5.3_Resp_Vis.Rdata")


