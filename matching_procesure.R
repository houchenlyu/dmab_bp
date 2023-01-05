###########################################################
###################  Matching procesure ###################
###########################################################

dmab # dataset of qualified denosumab users
control_full # dataset of qualified bisphosphoante users
# arrange by index data, so as to start with the earliest denosumab user
dmab <- dmab %>% arrange(index_date)
# save the caliper value
caliper_width <- 0.2*sd(result$log_ps)

### create datasets
# shared_control is dataset including all controls
# selected_control is the dataset including all selected controls
# used_list is a vector saving all the ID of the selected controls.
shared_control <- control_full 
selected_control <- NULL
used_list <- NULL

# covariates:
# subjid: ID;
# log_ps: log(ps)

### start to select ###
time1 <- Sys.time()
i <- 1
end <- nrow(dmab) # n = 4350
for (i in 1 : end) {
  # select the data for each denosumab user;
  single <- dmab %>% filter(subjid %in% dmab$subjid[i]) %>% mutate(log_ps_dmab = log_ps)
  # select the subjects in each exposure set for each dmab user;
  control <- shared_control %>% filter(matched_ID %in% dmab$subjid[i]) 
  # add log_ps_dmab to controls;
  control <- control %>% mutate(log_ps_dmab = single$log_ps_dmab)
  # create log_dif which measures the difference between log_ps_dmab and the log_ps of all possible controls in this exposure set;
  control <- control %>% mutate(log_dif = abs(log_ps_dmab - log_ps)) 
  # arrange the controls by value of log_dif (absolute value from small to large)
  control <- control %>% arrange(log_dif)
  # filter those within caliper with of 0.2SD of the PS on log scale;
  control <- control %>% filter(log_dif<= caliper_width)
  # select the best five controls if possible
  select <- control %>% slice(1:5)
  # save the used ID of selected controls
  used_list <- c(used_list, select$subjid)
  # reset single and control
  single <- NULL
  control <- NULL
  # save selected controls
  selected_control <- bind_rows(selected_control, select)
  # show the progress of the loop
  indicator <- round(as.numeric(100*i/end), digits = 8)
  print(paste0("seq=", i, ", ",indicator," %, ", Sys.time()-time1)) 
}
### end of selection ###
# this dataset “selected_control” includes best controls for all denosumab users
save.RDS(selected_control,"selected_control.RDS")
###########################################################
