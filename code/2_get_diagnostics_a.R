source('code/0_set_params.R')

# Diagnostics tables will be saved in these objects:
# Based on Carvalho et al. (2021, 10.1016/j.fishres.2021.105959)
save_t0 = list() # Likelihood and convergence
save_t1 = list() # Residual diagnostics CPUE
save_t2 = list() # Residual diagnostics mean length
save_t3 = list() # RMSE CPUE
save_t4 = list() # RMSE mean length
save_t7 = list() # Trend in recdevs
save_t7b = list() # Rec devs for plotting
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
  
  # Trends in recdevs (process error):
  # See: 10.1016/j.fishres.2022.106478
  rec_vec = i_mod$recruit %>% dplyr::filter(era == 'Main') %>% na.omit
  ltest = funtimes::notrend_test(x = rec_vec$dev)
  i_test7 = data.frame(species = sel_sp, modname = all_grids[i]) 
  i_test7$pval = ltest$p.value
  
  # Save rec devs:
  i_test7b = data.frame(Yr = rec_vec$Yr, dev = rec_vec$dev,
                        species = sel_sp, modname = all_grids[i]) 
  
  # Save all test:
  save_t0[[i]] = i_test0
  save_t1[[i]] = i_test1
  save_t2[[i]] = i_test2
  save_t3[[i]] = i_test3
  save_t4[[i]] = i_test4
  save_t7[[i]] = i_test7
  save_t7b[[i]] = i_test7b
  
  print(i)
  
}

# Save outputs:
saveRDS(dplyr::bind_rows(save_t0), file = file.path(output_folder, 'convergence.rds'))
saveRDS(dplyr::bind_rows(save_t1), file = file.path(output_folder, 'run_test_index.rds'))
saveRDS(dplyr::bind_rows(save_t2), file = file.path(output_folder, 'run_test_len.rds'))
saveRDS(dplyr::bind_rows(save_t3), file = file.path(output_folder, 'rmse_index.rds'))
saveRDS(dplyr::bind_rows(save_t4), file = file.path(output_folder, 'rmse_len.rds'))
saveRDS(dplyr::bind_rows(save_t7), file = file.path(output_folder, 'recdev_trend.rds'))
saveRDS(dplyr::bind_rows(save_t7b), file = file.path(output_folder, 'recdev_trend_df.rds'))

# -------------------------------------------------------------------------


# Plot recdevs time series
save_t7b = readRDS(file.path(output_folder, 'recdev_trend_df.rds'))
plot_data = save_t7b

p1 = ggplot(data = plot_data, aes(x = Yr, y = dev)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() + ylab('Recruitment deviates') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'bottom') +
  facet_wrap(~ modname)
ggsave(filename = file.path(output_folder, 'rec_devs.png'), plot = p1, 
       width = 200, height = 220, units = 'mm', dpi = 400)