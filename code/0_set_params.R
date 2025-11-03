# Install ss3diags if needed:
# remotes::install_github("nmfs-ost/ss3diags")
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
require(ggplot2)
require(future)
# You will need to set your WD first:
source('create_retro_files.r')
source('run_ASPM.R')

# Specify path where SS exe is located:
ss_folder = 'C:/Use/OneDrive - AZTI/Assessment_models'
# SS executable name:
ss_exe = 'ss.exe'

# Time step in your SS models:
ts_model = 'quarter' # quarter or year

# Folder where diagnostics results will be saved:
output_folder = 'C:/Use/OneDrive - AZTI/Assessment_models/IOTC/2025/BET/diagnostics_LS'
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