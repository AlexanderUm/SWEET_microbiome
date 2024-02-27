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
               "Maaslin2", "ComplexHeatmap", "circlize", "splinectomeR")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")


#-------------------------------------------------------------------------------
# Define responders and non-responders 
#-------------------------------------------------------------------------------
meta.resp <- meta.ls$all %>% 
                select(c("Body.weight", "HbA1c", "Fasting.glucose..mmol.L.", 
                         "Subject", "CID", "Group")) %>% 
                rename(Body_weight = Body.weight, 
                       Fasting_glucose = Fasting.glucose..mmol.L.) %>% 
                mutate(across(c("Body_weight", "HbA1c", "Fasting_glucose"), 
                              as.numeric))

meta.resp.long <- meta.resp %>% 
                  pivot_longer(all_of(c("Body_weight", 
                                        "HbA1c", 
                                        "Fasting_glucose")))

ggplot(meta.resp.long) + 
  geom_line(aes(x = CID, y = value, group = Subject)) + 
  facet_grid(name ~ Group, scales = "free") + 
  theme_bw()
