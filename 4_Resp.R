#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
load("PRM.Rdata")

load(paste0(PRM$data$out_dir, "/0_data.Rdata"))

set.seed(PRM$general$seed)

# Load libraries 
for (i in PRM$general$libs) {library(i, character.only = TRUE)}

# Custom function 
source("R/response_functions.R")

# Create directory 
DirOut <- PRM$resp$out_dir

dir.create(paste0(DirOut, "/plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(DirOut, "/tabs"), recursive = TRUE, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# Define responders and non-responders 
#-------------------------------------------------------------------------------
# Transform values log was tested but resulted in mostly identical outcomes. 
#-------------------------------------------------------------------------------
# Make a dataframe with indexes and response notations 

RespDf <- DataComb$all$meta %>% 
              reframe(Subject = unique(Subject), 
                      .by = Group)

for(i in names(PRM$resp$calc)) { 
  
  InstCol <- i
  
  InstIndexFun <- get(PRM$resp$calc[[i]]$index_fun)
  
  InstRespFun <- get(PRM$resp$calc[[i]]$resp_fun)
  
  RespDf <- DataComb$all$meta %>% 
                arrange(Subject, CID) %>% 
                reframe(across(all_of(InstCol), 
                                ~ InstIndexFun(.x, CID),
                               .names = "{.col}_index"), 
                        .by = c(Subject, Group)) %>% 
                mutate(across(all_of(paste0(InstCol, "_index")),
                               ~ InstRespFun(.x),
                              .names = "Resp_{.col}"), 
                       .by = Group) %>% 
                select(-Group) %>% 
                left_join(RespDf, ., by = "Subject")
}


#-------------------------------------------------------------------------------
# Plot responders and non-responders
#-------------------------------------------------------------------------------
RespColOrig <- names(PRM$resp$calc)

RespDfPlot <- RespDf %>% 
                  rename_with(~gsub("Resp_", "Response--", .x), 
                              starts_with("Resp")) %>% 
                  rename_with(~paste0("Index--", .x), 
                              all_of(paste0(RespColOrig, "_index"))) %>% 
                  pivot_longer(cols = -c(Subject, Group), 
                               names_to = c(".value", "Name"), 
                               names_sep = "--" ) %>% 
                  mutate(Name = gsub("_", " ", Name))


RespSumPlot <- ggplot(RespDfPlot,
                      aes(x = Group,
                          y = Index)) +
                  geom_jitter(aes(color = Response), 
                              width = 0.2,
                              height = 0,
                              alpha = 0.9, 
                              size = 1) +
                  geom_violin(fill = NA) +
                  facet_wrap(. ~ Name, scales = "free") +
                  theme_bw() + 
                  theme(axis.title.x = element_blank())

ggsave(filename = paste0(DirOut, "/Response_ind.svg"), 
       plot = RespSumPlot, width = 7, height = 3.5)

ggsave(filename = paste0(DirOut, "/Response_ind.png"), 
       plot = RespSumPlot, width = 7, height = 3.5)


#===============================================================================
# Formatting and writing out data about responders and non-responders 
# for publication: datafigures -> extanded figur 1A. 

bind_rows(
  RespDfPlot %>% 
    filter(!is.na(Index)) %>% 
    summarise(Mean = mean(Index), 
              SD = sd(Index), 
              .by = c(Group, Name)) %>% 
    mutate(Response = "Overall"),

  RespDfPlot %>% 
    filter(!is.na(Index)) %>% 
    summarise(Mean = mean(Index), 
              SD = sd(Index), 
              .by = c(Group, Name, Response))
) %>% 
  mutate(across(c(Mean, SD), function(x){sprintf("%.3f", round(x, 3))})) %>% 
  mutate(Experssion = paste0(Mean, "\u00B1", SD)) %>% 
  filter(!is.na(Response)) %>% 
  select(Name, Group, Response, Experssion) %>% 
  arrange(Name, Group, Response) %>% 
  pivot_wider(names_from = Response, values_from = Experssion) %>% 
  select(Name, Group, Overall, Resp, NonResp) %>% 
  write_tsv(paste0(DirOut, "/Response_ind_summary.tsv"))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Add Response data to original metadata 
#-------------------------------------------------------------------------------
MetaLsResp <- list()

RespDfFilt <- RespDf %>% 
                select(-Group)

for(i in names(DataComb)) { 
  
  MetaLsResp[[i]] <- DataComb[[i]]$meta %>% 
                          left_join(., RespDfFilt, by = "Subject") %>% 
                          mutate(RoWnAmEs = rownames(DataComb[[i]]$meta)) %>% 
                          mutate(across(starts_with("Resp"), as.factor)) %>% 
                          column_to_rownames(var = "RoWnAmEs")
}


#-------------------------------------------------------------------------------
# Write out results 
#-------------------------------------------------------------------------------
write.csv(RespDf, 
          file = paste0(DirOut, "/Response_ind.csv"))

RespData <- list("Resp_data" = RespDf, 
                 "Resp_meta_comb" = MetaLsResp)

save(list = c("RespData"), 
     file = paste0(PRM$data$out_dir, "/4_Resp.Rdata"))

rm(list = ls())
gc()
