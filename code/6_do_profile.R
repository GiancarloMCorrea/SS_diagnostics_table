source('code/0_set_params.R')

# -------------------------------------------------------------------------
# 5) Do profile (one dimensional):

# this_par = 'R0'
# par_vec = round(seq(from = 8, to = 11, length.out = n_cores), digits = 3)

this_par = 'NatM'
par_vec = round(seq(from = 0.2, to = 0.45, length.out = 10), digits = 3)
nCores = min(n_cores, length(par_vec))

for(i in seq_along(all_grids)) {
  
  dir_prof = file.path(grid_folder, all_grids[i], paste0("profile_", this_par))
  copy_SS_inputs(
    dir.old = file.path(grid_folder, all_grids[i]),
    dir.new = dir_prof,
    create.dir = TRUE,
    overwrite = TRUE,
    copy_par = TRUE,
    verbose = TRUE
  )
  
  # Modify starter
  starter <- SS_readstarter(file.path(dir_prof, "starter.ss"))
  starter[["ctlfile"]] <- "control_modified.ss"
  starter[["prior_like"]] <- 1
  starter[["init_values_src"]] = 0
  starter[["jitter_fraction"]] = 0
  starter[["parmtrace"]] = 0
  SS_writestarter(starter, dir = dir_prof, overwrite = TRUE)
  
  future::plan(future::multisession, workers = nCores)
  prof.table <- profile(
    dir = dir_prof,
    oldctlfile = "control.ss",
    newctlfile = "control_modified.ss",
    string = this_par, 
    profilevec = par_vec,
    exe = file.path(ss_folder, ss_exe),
    extras = '-nohess'
  )
  future::plan(future::sequential)
  
}


# -------------------------------------------------------------------------
# Make plots:

save_res = list()
for(i in seq_along(all_grids)) {
  
  # Get directory:
  dir_prof = file.path(grid_folder, all_grids[i], paste0("profile_", this_par))
  
  # Get outputs:
  profilemodels <- SSgetoutput(
    dirvec = dir_prof,
    keyvec = 1:length(par_vec), getcovar = FALSE
  )

  # Add original model:
  MLEmodel = SS_output(dir = file.path(grid_folder, all_grids[i]), verbose = FALSE, printstats = FALSE)
  profilemodels[["MLE"]] <- MLEmodel
  
  # Get summary
  profilesummary <- SSsummarize(profilemodels)
  
  # plot profile using summary created above
  results <- SSplotProfile(profilesummary,
                           profile.string = this_par, 
                           profile.label = this_par,
                           plotdir = file.path(output_folder, all_grids[i]),
                           print = TRUE) 
  save_res[[i]] = results %>% mutate(species = sel_sp, modname = all_grids[i])
  
  # Change name of png file created:
  file.rename(from = file.path(output_folder, all_grids[i], 'profile_plot_likelihood.png'),
              to = file.path(output_folder, all_grids[i], paste0('profile_plot_likelihood_', this_par, '.png')))
  
  # Make summary time series:
  # SSplotComparisons(profilesummary, legendlabels = paste(this_par, "=", par_vec))
  
  print(i)
}

saveRDS(dplyr::bind_rows(save_res), file = file.path(output_folder, paste0('profile_df_', this_par, '.rds')))
