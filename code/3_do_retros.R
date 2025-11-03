source('code/0_set_params.R')

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
  i_retro_dir = create_retro_files(dir = file.path(grid_folder, all_grids[i]), years = -1:n_peels, mod_ts = ts_model)
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

