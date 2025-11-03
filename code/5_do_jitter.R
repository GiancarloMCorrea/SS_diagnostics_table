source('code/0_set_params.R')
numjitter = 20 # define number of iterations

# -------------------------------------------------------------------------
# run jitter:

# Run in parallel:
for(i in seq_along(all_grids)) {

  # Create a new folder in the model folder and paste all files in there 
  # This is done to avoid overwriting the input files in original folder
  dir.create(path = file.path(file.path(grid_folder, all_grids[i]), "jitter"), showWarnings = FALSE)
  file.copy(from = list.files(path = file.path(grid_folder, all_grids[i]), pattern="\\.", full.names = TRUE), 
            to = file.path(grid_folder, all_grids[i], "jitter"))
  
  future::plan(future::multisession, workers = n_cores)
  jit.likes <- jitter(
    dir = file.path(grid_folder, all_grids[i], "jitter"), Njitter = numjitter, # change number jitters
    jitter_fraction = 0.1, init_values_src = 1, exe = file.path(ss_folder, ss_exe),
    extras = '-nohess'
  )
  future::plan(future::sequential)
  
}


# -------------------------------------------------------------------------
# Extract outputs:

save_t8 = list() # jitter like
save_t9 = list() # jitter ts
for(i in seq_along(all_grids)) {
  
  # Jitter results:
  profilemodels = SSgetoutput(dirvec = file.path(grid_folder, all_grids[i], "jitter"), 
                              keyvec = c('', 1:numjitter), getcovar = FALSE,
                              getcomp = FALSE, verbose = FALSE)
  profilesummary = SSsummarize(profilemodels)
  # Likelihood:
  i_test8 = profilesummary[["likelihoods"]] %>% 
    tidyr::pivot_longer(!Label, names_to = 'Iter', values_to = "Likelihood") %>% 
    mutate(species = sel_sp, modname = all_grids[i])
  save_t8[[i]] = i_test8
  # Time series:
  i_test9 = profilesummary[["SpawnBio"]] %>% 
    tidyr::pivot_longer(!c(Yr, Label), names_to = 'Iter', values_to = "value") %>% 
    mutate(species = sel_sp, modname = all_grids[i])
  save_t9[[i]] = i_test9
  
  print(i)
  
}
saveRDS(dplyr::bind_rows(save_t8), file = file.path(output_folder, 'jitter-like.rds'))
saveRDS(dplyr::bind_rows(save_t9), file = file.path(output_folder, 'jitter-ts.rds'))


# -------------------------------------------------------------------------

# Plot jitter results:
save_t8 = readRDS(file.path(output_folder, 'jitter-like.rds'))
save_t9 = readRDS(file.path(output_folder, 'jitter-ts.rds'))

# jitter likelihood:
plot_data = save_t8 %>% filter(Label == 'TOTAL') %>% mutate(Type = if_else(Iter == 'replist', 'Base', 'Jitter'))

p1 = ggplot(data = plot_data, aes(x = Iter, y = Likelihood)) +
  geom_point(aes(color = Type)) +
  scale_color_manual(values = c('red', 'gray50')) +
  theme_bw() + ylab('Total likelihood') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.position = 'bottom') +
  facet_wrap(~ modname, scales = 'free_y')
ggsave(filename = file.path(output_folder, 'jitter-like.png'), plot = p1, 
       width = 180, height = 130, units = 'mm', dpi = 400)

# SSB time series:
plot_data = save_t9 %>% mutate(Type = if_else(Iter == 'replist', 'Base', 'Jitter'))

p1 = ggplot(data = plot_data, aes(x = Yr, y = value)) +
  geom_line(aes(group = Iter, color = Type), alpha = 0.5) +
  scale_color_manual(values = c('red', 'gray50')) +
  theme_bw() + ylab('SSB') + xlab("Time") +
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  facet_wrap(~ modname, scales = 'free_y')
ggsave(filename = file.path(output_folder, 'jitter-ssb.png'), plot = p1, 
       width = 180, height = 130, units = 'mm', dpi = 400)
