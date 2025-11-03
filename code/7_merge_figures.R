source('code/0_set_params.R')

all_models = list.dirs(path = output_folder, full.names = FALSE, recursive = FALSE)

# -------------------------------------------------------------------------

# Select figure name:
var_vector = c("retro_SSB", "retro_F", "retro_compare9_recruits", "retro_compare3_Bratio")

for(k in seq_along(var_vector)) {
  this_fig = var_vector[k]  
  png(filename = file.path(output_folder, paste0(this_fig, '.png')), width = 160, height = 120, 
      units = 'mm', res = 400)
  
  par(mfrow = c(2,3))
  
  for(i in seq_along(all_models)) {
    img1 = readPNG(file.path(output_folder, all_models[i], paste0(this_fig, '.png')))
    par(mar = c(0.2, 0.2, 0.65, 0.2))
    plot(NA, xlim = 0:1, ylim = 0:1, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    rasterImage(img1, 0, 0, 1, 1)
    title(all_models[i], cex.main = 1)
  }
  
  dev.off()
}


# -------------------------------------------------------------------------
# Likelihood profiles:

this_par = 'R0'  # R0 NatM

png(filename = file.path(output_folder, paste0('profile_plot_likelihood_', this_par, '.png')), 
    width = 160, height = 90, units = 'mm', res = 400)
par(mfrow = c(2,3))

for(i in seq_along(all_models)) {
  img1 = readPNG(file.path(output_folder, all_models[i], paste0('profile_plot_likelihood_', this_par, '.png')))
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot(NA, xlim = 0:1, ylim = 0:1, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  rasterImage(img1, 0, 0, 1, 1)
}

dev.off()

