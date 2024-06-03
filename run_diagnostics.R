# Install ss3diags if needed:
# remotes::install_github("jabbamodel/ss3diags")

# IMPORTANT: code only works on Windows

# Load libraries:
require(ss3diags)
require(r4ss)
require(dplyr)
require(doSNOW)
require(parallel)
require(funtimes)
# You will need to set your WD first:
source('create_retro_files.r')

# Specify path where SS exe is located:
ss_folder = 'C:/Use/OneDrive - AZTI/Codes'
# SS executable name:
ss_exe = 'ss.exe'

# Folder where diagnostics results will be saved:
output_folder = 'C:/Use/OneDrive - AZTI/Assessment_models/ICCAT/2024/YFT/diagnostics'
dir.create(output_folder)

# SS model grid folder
# All your models of interest are saved here, exclude other files or folders
grid_folder = 'C:/Use/OneDrive - AZTI/Assessment_models/ICCAT/2024/YFT'

# Read all models in your grid folder:
all_grids = list.files(grid_folder)
# You could also just index the desired folders 
# e.g., all_grids = list.files(grid_folder)[c(2,3,4)]

# Number of peels for retrospective analysis:
n_peels = -5 # should be negative

# Species code:
sel_sp = 'YFT'


# -------------------------------------------------------------------------
# 1) Check if some models in the uncertainty grid need to be run:

for(i in seq_along(all_grids)) {
  # Look for Report.sso file. If not found, run model:
  is_rep = 'Report.sso' %in% list.files(file.path(grid_folder, all_grids[i]))
  if(!is_rep) {
    r4ss::run(dir = file.path(grid_folder, all_grids[i]), exe = file.path(ss_folder, ss_exe))
    cat(file.path(grid_folder, all_grids[i]), " done", "\n")
  }
}

# -------------------------------------------------------------------------
# 2) Run retrospective analysis (in parallel):
# You can adapt this part and do it as you prefer.

# Create retro folders:
all_retro_dir = NULL
for(i in seq_along(all_grids)) {
  i_retro_dir = create_retro_files(dir = file.path(grid_folder, all_grids[i]), years = 0:n_peels)
  all_retro_dir = c(all_retro_dir, i_retro_dir)
}

# Detect number of models:
nSim = length(all_retro_dir)

# Specify number of cores:
# WARNING: I recommend to use max 80% of your total number of cores
# NEVER use all of them
nCores = 2
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
# 3) Calculate diagnostics:

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
for(i in seq_along(all_grids)) {

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
  ss3diags::SSplotRunstest(ss3rep = i_mod, print = TRUE, plotdir = file.path(grid_folder, all_grids[i]), 
                           subplots = 'cpue', filenameprefix = 'cpue')
  ss3diags::SSplotRunstest(ss3rep = i_mod, print = TRUE, plotdir = file.path(grid_folder, all_grids[i]), 
                           subplots = 'len', filenameprefix = 'len')
  
  # RMSE (goodness-of-fit), and save images in the models' folder:
  i_test3 = ss3diags::SSplotJABBAres(i_mod, subplots = 'cpue', plotdir = file.path(grid_folder, all_grids[i]),
                                     print = TRUE, filenameprefix = 'cpue')
  i_test3 = i_test3 %>% mutate(species = sel_sp, modname = all_grids[i])
  
  i_test4 = ss3diags::SSplotJABBAres(i_mod, subplots = 'len', plotdir = file.path(grid_folder, all_grids[i]),
                                     print = TRUE, filenameprefix = 'len')
  i_test4 = i_test4 %>% mutate(species = sel_sp, modname = all_grids[i])
  
  graphics.off() # to close plot windows open
  
  # Read retros:
  subdirs = list.files(path = file.path(grid_folder, all_grids[i], 'retrospectives'))
  subdirs = subdirs[c(length(subdirs), 1:(length(subdirs)-1))] # reorder
  retroModels = SSgetoutput(dirvec = file.path(grid_folder, all_grids[i], 'retrospectives', subdirs))
  retros = SSsummarize(retroModels)
  
  # MASE (prediction skill):
  i_test5 = ss3diags::SSmase(retros, quants = 'cpue')
  i_test5 = i_test5 %>% mutate(species = sel_sp, modname = all_grids[i])
  ss3diags::SSplotHCxval(retroSummary = retros, subplots = 'cpue', print = TRUE, 
                         plotdir = file.path(grid_folder, all_grids[i]), filenameprefix = 'cpue')
  
  # MASE (for mean length):
  # TODO: check if it works
  
  # Monhn rho, retrospective analysis (model consistency):
  i_test6 = rbind(ss3diags::SShcbias(retros, quants = 'SSB'),
                  ss3diags::SShcbias(retros, quants = 'F'))
  i_test6 = i_test6 %>% mutate(species = sel_sp, modname = all_grids[i])
  SSplotRetro(summaryoutput = retros, subplots = 'SSB', print = TRUE, plotdir = file.path(grid_folder, all_grids[i]))
  SSplotRetro(summaryoutput = retros, subplots = 'F', print = TRUE, plotdir = file.path(grid_folder, all_grids[i]))
  
  # Trends in recdevs (process error):
  # See: 10.1016/j.fishres.2022.106478
  ltest = funtimes::notrend_test(x = i_mod$recruit$dev[!is.na(i_mod$recruit$dev)])
  i_test7 = data.frame(species = sel_sp, modname = all_grids[i]) 
  i_test7$pval = ltest$p.value
  
  # Save all test:
  save_t0[[i]] = i_test0
  save_t1[[i]] = i_test1
  save_t2[[i]] = i_test2
  save_t3[[i]] = i_test3
  save_t4[[i]] = i_test4
  save_t5[[i]] = i_test5
  save_t6[[i]] = i_test6
  save_t7[[i]] = i_test7
  
  print(i)
  
}

# Merge outputs:
all_t0 = dplyr::bind_rows(save_t0)
all_t1 = dplyr::bind_rows(save_t1)
all_t2 = dplyr::bind_rows(save_t2)
all_t3 = dplyr::bind_rows(save_t3)
all_t4 = dplyr::bind_rows(save_t4)
all_t5 = dplyr::bind_rows(save_t5)
all_t6 = dplyr::bind_rows(save_t6)
all_t7 = dplyr::bind_rows(save_t7)

# Save outputs:
dir.create(path = file.path(output_folder, sel_sp))
save(all_t0, file = file.path(output_folder, sel_sp, 'test0.RData'))
save(all_t1, file = file.path(output_folder, sel_sp, 'test1.RData'))
save(all_t2, file = file.path(output_folder, sel_sp, 'test2.RData'))
save(all_t3, file = file.path(output_folder, sel_sp, 'test3.RData'))
save(all_t4, file = file.path(output_folder, sel_sp, 'test4.RData'))
save(all_t5, file = file.path(output_folder, sel_sp, 'test5.RData'))
save(all_t6, file = file.path(output_folder, sel_sp, 'test6.RData'))
save(all_t7, file = file.path(output_folder, sel_sp, 'test7.RData'))