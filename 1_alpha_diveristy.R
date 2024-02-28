#-------------------------------------------------------------------------------
# Load libraries and set the environment 
#-------------------------------------------------------------------------------
set.seed(349574)

# Load libraries 
libs.list <- c("phyloseq", "tidyverse", "picante", "MicrobiomeStat")

for (i in libs.list) {library(i, character.only = TRUE)}

rm(list = c("i", "libs.list"))

load("out/supp/data_bundel.Rdata")



#-------------------------------------------------------------------------------
# Alpha diversity
#-------------------------------------------------------------------------------
var.grid <- expand.grid(alpha.vars.ls$data_set, 
                        alpha.vars.ls$used_ps) %>% 
              mutate(across(everything(), as.character))


alpha.res.ls <- list()

alpha.plot.ls <- list()


for(i.grid in 1:nrow(var.grid)) {
  
  i.gr <- var.grid[i.grid, 1]
  i.lvl <- var.grid[i.grid, 2]
  
  # Extract data
  ps.inst <- pss.ls[[i.gr]][[i.lvl]] 
    
  otu.inst <- ps.inst %>% 
                  otu_table()

  # Calculate alpha diversity with mStat
  alpha.obj <- mStat_calculate_alpha_diversity(otu.inst, 
                                                 alpha.vars.ls$alpha_ind) 
    
  # Mstat object
  mstat.obj <- mStat_convert_phyloseq_to_data_obj(ps.inst)
    
  # Change to class 
  mstat.obj$meta.dat <- mstat.obj$meta.dat %>% 
              mutate(!!alpha.vars.ls$time_var := 
                       as.numeric(.data[[alpha.vars.ls$time_var]]), 
                     across(c(alpha.vars.ls$subject_var, 
                              alpha.vars.ls$group_var, 
                              alpha.vars.ls$adjust_var), as.factor))
    
    # Trends
  alpha.mstat <- generate_alpha_trend_test_long(
                                    data.obj = mstat.obj, 
                                    alpha.obj = alpha.obj,
                                    alpha.name = alpha.vars.ls$alpha_ind,
                                    time.var = alpha.vars.ls$time_var,
                                    subject.var = alpha.vars.ls$subject_var,
                                    group.var = alpha.vars.ls$group_var,
                                    adj.vars = alpha.vars.ls$adjust_var)
    
  alpha.res.ls[[i.gr]][[i.lvl]] <- alpha.mstat
    
  # Plot alpha 
  alpha.df <- as.data.frame(alpha.obj) %>% 
                      bind_cols(., mstat.obj$meta.dat) %>% 
                      pivot_longer(cols = alpha.vars.ls$alpha_ind, 
                                   names_to = "Index", 
                                   values_to = "Value")
    
  alpha.df.s <- alpha.df %>% 
                        group_by(across(c(alpha.vars.ls$time_var, 
                                          alpha.vars.ls$group_var, 
                                          "Index"))) %>% 
                        summarise(sum_val = mean(Value), 
                                  .groups = "keep")
    
    
  # Add statistics to plots 
  stat.text <- data.frame(Index=as.character(), 
                            Text=as.character())
    
  for(i.ind in names(alpha.mstat)) {
      
      alpha.test <- alpha.mstat[[i.ind]] 
      
      text <- paste0("Group(gr2):Time [P=", 
                     round(alpha.test[nrow(alpha.test), "P.Value"], 3), 
                     "; Est=", 
                     round(alpha.test[nrow(alpha.test), "Estimate"], 3), "]")
      
      stat.text <-  c("Index" = i.ind, "Text" = text) %>% 
        bind_rows(stat.text, .)
      
  } 
    
  sig.df <- alpha.df %>% 
                  group_by(across(c("Index"))) %>% 
                  summarise(sum_val = max(Value), 
                            .groups = "keep") %>% 
                  mutate(y_text = sum_val*1.1, 
                         y_fpoint = sum_val*1.2, 
                         !!alpha.vars.ls$time_var := 
                           sort(alpha.df[[alpha.vars.ls$time_var]])[1]) %>% 
                  left_join(., stat.text, by = "Index")
    
    
    
  alpha.plot <- ggplot(alpha.df) + 
                        geom_line(aes(x = .data[[alpha.vars.ls$time_var]], 
                                      y = Value, 
                                      group = .data[[alpha.vars.ls$subject_var]],
                                      color = .data[[alpha.vars.ls$group_var]]), 
                                  alpha = 0.25) + 
                        geom_line(data = alpha.df.s, 
                                  aes(x = .data[[alpha.vars.ls$time_var]], 
                                      y = sum_val, 
                                      group = .data[[alpha.vars.ls$group_var]],
                                      color = .data[[alpha.vars.ls$group_var]]), 
                                  linewidth = 1) +
                        geom_point(data = alpha.df.s, 
                                   aes(x = .data[[alpha.vars.ls$time_var]], 
                                       y = sum_val), 
                                   size = 1) +
                        geom_text(data = sig.df, 
                                  aes(label = Text,
                                      y = y_text, 
                                      x = Time), size =2.25, hjust = 0) +
                        facet_wrap(~Index, scales = "free_y") + 
                        theme_bw() + 
                        scale_color_manual(values = aest.ls$color_gr) + 
                        scale_x_continuous(breaks = c(0, 2, 6, 12)) + 
                        xlab("Time(Month)")
    
  alpha.plot.ls[[i.gr]][[i.lvl]] <- alpha.plot
    
}

save(list = c("alpha.res.ls", "alpha.plot.ls"), 
     file = "out/supp/alpha_res.Rdata")
