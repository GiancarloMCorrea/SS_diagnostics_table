# Install ss3diags if needed:
# remotes::install_github("jabbamodel/ss3diags")
# remotes::install_github("r4ss/r4ss")
rm(list = ls())

# IMPORTANT: code only works on Windows

# Load libraries:
require(ss3diags)
require(r4ss)
require(dplyr)
require(doSNOW)
require(parallel)
require(funtimes)
require(future)
# You will need to set your WD first:
source('create_retro_files.r')
source('run_ASPM.R')

# Specify path where SS exe is located:
ss_folder = 'C:/Use/OneDrive - AZTI/Assessment_models'
# SS executable name:
ss_exe = 'ss.exe'

# Folder where diagnostics results will be saved:
output_folder = 'C:/Use/OneDrive - AZTI/Assessment_models/IOTC/2025/BET/diagnostics'
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# SS model grid folder
# All your models of interest are saved here, exclude other files or folders
grid_folder = 'C:/Use/OneDrive - AZTI/Assessment_models/IOTC/2025/UpdatedGridBET2025'

# Read all models in your grid folder:
all_grids = list.files(grid_folder)
# You could also just index the desired folders 
# e.g., all_grids = list.files(grid_folder)[c(2,3,4)]

# Number of peels for retrospective analysis:
n_peels = -5 # should be negative

# Species code:
sel_sp = 'BET'

# Maximum number of cores to use in analyses:
n_cores = 10


# -------------------------------------------------------------------------
# 1) Check if some models in the uncertainty grid need to be run:

models_to_run = NULL
for(i in seq_along(all_grids)) {
  # Look for Report.sso file. If not found, run model:
  is_rep = 'Report.sso' %in% list.files(file.path(grid_folder, all_grids[i]))
  if(!is_rep) {
    models_to_run = c(file.path(grid_folder, all_grids[i]), models_to_run)
  }
}

if(is.null(models_to_run)) {
  cat("All models are run", "\n")
} else {
  if(length(models_to_run) == 1) { # only one model, run it
    r4ss::run(dir = file.path(grid_folder, models_to_run), exe = file.path(ss_folder, ss_exe))
    cat(file.path(grid_folder, models_to_run), " done", "\n")
  } else { # run models in parallel
    nSim = length(models_to_run)
    nCores = min(n_cores, nSim)
    cl = makeCluster(nCores)
    registerDoSNOW(cl)
    # Run in parallel:
    foreach(ix = 1:nSim) %dopar% {
      require(r4ss)
      r4ss::run(dir = models_to_run[ix], exe = file.path(ss_folder, ss_exe))
    }
    # Stop cluster:
    stopCluster(cl)
  }
}

# -------------------------------------------------------------------------
# 2) Run retrospective analysis (in parallel):
# You can adapt this part and do it as you prefer.
# IMPORANT: This might take a while, run it overnight or over a weekend

# Create retro folders:
all_retro_dir = NULL
for(i in seq_along(all_grids)) {
  orig_files = list.files(path = file.path(grid_folder, all_grids[i]))
  dir.create(file.path(grid_folder, all_grids[i], "retrospectives", "retro0"), showWarnings = FALSE, recursive = TRUE)
  # Copy and paste all files for retro year == 0 to save some time:
  file.copy(paste (file.path(grid_folder, all_grids[i]), orig_files, sep = "/"), 
            paste (file.path(grid_folder, all_grids[i], "retrospectives", "retro0"), orig_files, sep = "/"), 
            overwrite = TRUE, recursive = FALSE, 
            copy.mode = TRUE)
  # Copy files for -1 to n_peels
  i_retro_dir = create_retro_files(dir = file.path(grid_folder, all_grids[i]), years = -1:n_peels, mod_ts = "quarter")
  all_retro_dir = c(all_retro_dir, i_retro_dir)
}

# Detect number of models:
nSim = length(all_retro_dir)

# Specify number of cores:
# If you do not how many available cores you have, then run:
# detectCores()
# WARNING: I recommend to use max 80% of your total number of cores
# NEVER use all of them
nCores = n_cores
cl = makeCluster(nCores)
registerDoSNOW(cl)

# Run in parallel:
foreach(ix = 1:nSim) %dopar% {
  require(r4ss)
  r4ss::run(dir = all_retro_dir[ix], extras = '-nohess', exe = file.path(ss_folder, ss_exe))
}

# Stop cluster:
stopCluster(cl)


# -------------------------------------------------------------------------
# 3) run jitter:

numjitter = 20
for(i in seq_along(all_grids)) {
  
  future::plan(future::multisession, workers = n_cores)
  jit.likes <- jitter(
    dir = file.path(grid_folder, all_grids[i]), Njitter = numjitter, # change number jitters
    jitter_fraction = 0.1, init_values_src = 1, exe = file.path(ss_folder, ss_exe),
    extras = '-nohess'
  )
  future::plan(future::sequential)

}


# -------------------------------------------------------------------------
# 4) Run ASPM per model:

dir.create(file.path(output_folder, 'ASPM'))
# Detect number of models:
nSim = length(all_grids)
nCores = min(n_cores, nSim)
cl = makeCluster(nCores)
registerDoSNOW(cl)

# Run in parallel:
foreach(ix = 1:nSim) %dopar% {
  run_ASPM(mod_dir = file.path(grid_folder, all_grids[ix]), 
           aspm_dir =  file.path(output_folder, 'ASPM', all_grids[ix]),
           ss_dir = file.path(ss_folder, ss_exe))
}

# Stop cluster:
stopCluster(cl)


# -------------------------------------------------------------------------
# 5) Do profile (one dimensional):

# this_par = 'NatM'
# par_vec = round(seq(from = 0.2, to = 0.45, length.out = n_cores), digits = 3)

this_par = 'R0'
par_vec = round(seq(from = 8, to = 11, length.out = n_cores), digits = 3)

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
  SS_writestarter(starter, dir = dir_prof, overwrite = TRUE)
  
  future::plan(future::multisession, workers = nCores)
  prof.table <- profile(
    dir = dir_prof,
    oldctlfile = "control.ss",
    newctlfile = "control_modified.ss",
    string = this_par, 
    profilevec = par_vec,
    exe = file.path(ss_folder, ss_exe)
  )
  future::plan(future::sequential)

}

# -------------------------------------------------------------------------
# 5) Calculate diagnostics:

# Diagnostics tables will be saved in these objects:
# Based on Carvalho et al. (2021, 10.1016/j.fishres.2021.105959)
save_t0 = list() # Likelihood and convergence
save_t1 = list() # Residual diagnostics CPUE
save_t2 = list() # Residual diagnostics mean length
save_t3 = list() # RMSE CPUE
save_t4 = list() # RMSE mean length
save_t5 = list() # MASE CPUE
save_t6 = list() # Mohn rho
save_t7 = list() # Trends in recdevs
save_t8 = list() # jitter
for(i in seq_along(all_grids)) {
  
  # Create folder to save images:
  dir.create(file.path(output_folder, all_grids[i]))
  
  # Read model:
  i_mod = r4ss::SS_output(dir = file.path(grid_folder, all_grids[i]), 
                          verbose = FALSE, printstats = FALSE)
  
  # Convergence diagnostics:
  i_test0 = data.frame(species = sel_sp, modname = all_grids[i])
  i_test0$max_grad = i_mod$maximum_gradient_component
  i_test0$tot_nll = i_mod$likelihoods_used$values[1]
  i_test0$n_pars = i_mod$N_estimated_parameters
  i_test0$n_pars_bounds = length(which(!(i_mod$parameters$Status %in% c(NA, 'OK', 'act'))))
  i_test0$hessian_inv = !is.null(i_mod$stdtable)
  
  # Run tests for CPUE (goodness-of-fit):
  i_test1 = ss3diags::SSrunstest(ss3rep = i_mod, quants = 'cpue')
  i_test1 = i_test1 %>% mutate(species = sel_sp, modname = all_grids[i])
  
  # Run tests for mean length (goodness-of-fit):
  i_test2 = ss3diags::SSrunstest(ss3rep = i_mod, quants = 'len')
  i_test2 = i_test2 %>% mutate(species = sel_sp, modname = all_grids[i])
  
  # Make run tests plots, which will be saved in the models' folder:
  ss3diags::SSplotRunstest(ss3rep = i_mod, print = TRUE, plotdir = file.path(output_folder, all_grids[i]), 
                           subplots = 'cpue', filenameprefix = 'cpue')
  ss3diags::SSplotRunstest(ss3rep = i_mod, print = TRUE, plotdir = file.path(output_folder, all_grids[i]), 
                           subplots = 'len', filenameprefix = 'len')
  
  # RMSE (goodness-of-fit), and save images in the models' folder:
  i_test3 = ss3diags::SSplotJABBAres(i_mod, subplots = 'cpue', plotdir = file.path(output_folder, all_grids[i]),
                                     print = TRUE, filenameprefix = 'cpue')
  i_test3 = i_test3 %>% mutate(species = sel_sp, modname = all_grids[i])
  
  i_test4 = ss3diags::SSplotJABBAres(i_mod, subplots = 'len', plotdir = file.path(output_folder, all_grids[i]),
                                     print = TRUE, filenameprefix = 'len')
  i_test4 = i_test4 %>% mutate(species = sel_sp, modname = all_grids[i])
  
  graphics.off() # to close plot windows open
  
  # Read retros:
  subdirs = list.files(path = file.path(grid_folder, all_grids[i], 'retrospectives'))
  subdirs = subdirs[c(length(subdirs), 1:(length(subdirs)-1))] # reorder
  retroModels = SSgetoutput(dirvec = file.path(grid_folder, all_grids[i], 'retrospectives', subdirs))
  yrvec = numeric(length(subdirs))
  yrlab = character(length(subdirs))
  for(k in seq_along(subdirs)) {
    yrvec[k] = retroModels[[k]]$Retro_year
    yrlab[k] = as.character(retroModels[[k]]$Retro_year)
  }
  yrlab[1] = 'Ref'
  retros = SSsummarize(retroModels)
  
  # MASE (prediction skill):
  i_test5 = ss3diags::SSmase(retros, quants = 'cpue')
  i_test5 = i_test5 %>% mutate(species = sel_sp, modname = all_grids[i])
  ss3diags::SSplotHCxval(retroSummary = retros, subplots = 'cpue', print = TRUE, 
                         plotdir = file.path(output_folder, all_grids[i]), filenameprefix = 'cpue',
                         endyrvec = yrvec, legendlabels = yrlab)
  
  # MASE (for mean length):
  # TODO: check if it works
  
  # Monhn rho, retrospective analysis (model consistency):
  i_test6 = rbind(ss3diags::SShcbias(retros, quants = 'SSB'),
                  ss3diags::SShcbias(retros, quants = 'F'))
  i_test6 = i_test6 %>% mutate(species = sel_sp, modname = all_grids[i])
  SSplotRetro(summaryoutput = retros, subplots = 'SSB', print = TRUE, 
              plotdir = file.path(output_folder, all_grids[i]), endyrvec = yrvec,
              legendlabels = yrlab)
  SSplotRetro(summaryoutput = retros, subplots = 'F', print = TRUE, 
              plotdir = file.path(output_folder, all_grids[i]), endyrvec = yrvec,
              legendlabels = yrlab)
  
  # Trends in recdevs (process error):
  # See: 10.1016/j.fishres.2022.106478
  rec_vec = i_mod$recruit %>% dplyr::filter(era == 'Main') %>% na.omit
  ltest = funtimes::notrend_test(x = rec_vec$dev)
  i_test7 = data.frame(species = sel_sp, modname = all_grids[i]) 
  i_test7$pval = ltest$p.value
  
  # Jitter results:
  profilemodels = SSgetoutput(dirvec = file.path(grid_folder, all_grids[i]), 
                              keyvec = 1:numjitter, getcovar = FALSE,
                              getcomp = FALSE, verbose = FALSE)
  profilesummary = SSsummarize(profilemodels)
  i_test8 = profilesummary[["likelihoods"]] %>% tidyr::pivot_longer(!Label, names_to = 'Iter', values_to = "Likelihood")
  i_test8 = i_test8 %>% mutate(species = sel_sp, modname = all_grids[i])
  
  # Save all test:
  save_t0[[i]] = i_test0
  save_t1[[i]] = i_test1
  save_t2[[i]] = i_test2
  save_t3[[i]] = i_test3
  save_t4[[i]] = i_test4
  save_t5[[i]] = i_test5
  save_t6[[i]] = i_test6
  save_t7[[i]] = i_test7
  save_t8[[i]] = i_test8
  
  print(i)
  
}

# Save outputs:
saveRDS(dplyr::bind_rows(save_t0), file = file.path(output_folder, 'convergence.rds'))
saveRDS(dplyr::bind_rows(save_t1), file = file.path(output_folder, 'run_test_index.rds'))
saveRDS(dplyr::bind_rows(save_t2), file = file.path(output_folder, 'run_test_len.rds'))
saveRDS(dplyr::bind_rows(save_t3), file = file.path(output_folder, 'rmse_index.rds'))
saveRDS(dplyr::bind_rows(save_t4), file = file.path(output_folder, 'rmse_len.rds'))
saveRDS(dplyr::bind_rows(save_t5), file = file.path(output_folder, 'mase_index.rds'))
saveRDS(dplyr::bind_rows(save_t6), file = file.path(output_folder, 'retrospective.rds'))
saveRDS(dplyr::bind_rows(save_t7), file = file.path(output_folder, 'recdev_trend.rds'))
saveRDS(dplyr::bind_rows(save_t8), file = file.path(output_folder, 'jitter.rds'))