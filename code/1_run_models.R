source('code/0_set_params.R')

# -------------------------------------------------------------------------
# 1) Check if some models in the uncertainty grid need to be run:

# Identify models to run:
models_to_run = NULL
for(i in seq_along(all_grids)) {
  # Look for Report.sso file. If not found, run model:
  is_rep = 'Report.sso' %in% list.files(file.path(grid_folder, all_grids[i]))
  if(!is_rep) {
    models_to_run = c(file.path(grid_folder, all_grids[i]), models_to_run)
  }
}

# Run models:
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
      r4ss::run(dir = models_to_run[ix], exe = file.path(ss_folder, ss_exe), extras = '-nohess')
    }
    # Stop cluster:
    stopCluster(cl)
  }
}
