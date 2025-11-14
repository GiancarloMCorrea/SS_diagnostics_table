rm(list = ls())
require(r4ss)

# -------------------------------------------------------------------------

# Folder where your models are saved:
# Each subfolder is a different model configuration
grid_folder = 'C:/Use/OneDrive - AZTI/General - IOTC_YFT_2024/models/FinalGrid'

# Folder where the SS inputs will be saved:
# Subfolders with model configurations will be created
new_grid_folder = 'C:/Use/OneDrive - AZTI/Codes/YFT-Assessment-CPUE-update-2025/models_2024'
dir.create(new_grid_folder, showWarnings = FALSE, recursive = FALSE)

# Copy only input files:
all_models = list.files(grid_folder)
for(i in seq_along(all_models)) {
  copy_SS_inputs(dir.old = file.path(grid_folder, all_models[i]),
                 dir.new = file.path(new_grid_folder, all_models[i]))
}
