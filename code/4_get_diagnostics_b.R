source('code/0_set_params.R')

# Diagnostics tables will be saved in these objects:
# Based on Carvalho et al. (2021, 10.1016/j.fishres.2021.105959)
save_t5 = list() # MASE CPUE
save_t6 = list() # Mohn rho
for(i in seq_along(all_grids)) {
  
  # Create folder to save images:
  dir.create(file.path(output_folder, all_grids[i]))
  
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
  rho_est = r4ss::SSmohnsrho(summaryoutput = retros, endyrvec = yrvec)
  i_test6 = data.frame(type = c('SSB', 'F', 'Rec', 'Bratio'), 
                       Rho = c(rho_est$AFSC_Hurtado_SSB, rho_est$AFSC_Hurtado_F,
                               rho_est$AFSC_Hurtado_Rec, rho_est$AFSC_Hurtado_Bratio))
  i_test6 = i_test6 %>% mutate(species = sel_sp, modname = all_grids[i])
  SSplotRetro(summaryoutput = retros, subplots = 'SSB', print = TRUE, 
              plotdir = file.path(output_folder, all_grids[i]), endyrvec = yrvec,
              legendlabels = yrlab)
  SSplotRetro(summaryoutput = retros, subplots = 'F', print = TRUE, 
              plotdir = file.path(output_folder, all_grids[i]), endyrvec = yrvec,
              legendlabels = yrlab)
  # Plot Bratio and R retro:
  SSplotComparisons(retros, endyrvec = yrvec, new = FALSE, legendlabels = yrlab, subplots = 3, 
                    plotdir = file.path(output_folder, all_grids[i]), print = TRUE, filenameprefix="retro_")
  SSplotComparisons(retros, endyrvec = yrvec, new = FALSE, legendlabels = yrlab, subplots = 9,
                    plotdir = file.path(output_folder, all_grids[i]), print = TRUE, filenameprefix="retro_")
  
  # Save all test:
  save_t5[[i]] = i_test5
  save_t6[[i]] = i_test6

  print(i)
  
}

# Save outputs:
saveRDS(dplyr::bind_rows(save_t5), file = file.path(output_folder, 'mase_index.rds'))
saveRDS(dplyr::bind_rows(save_t6), file = file.path(output_folder, 'retrospective.rds'))
