run_ASPM = function(mod_dir, aspm_dir, ss_dir) {
  
  require(r4ss)
  require(dplyr)
  
  # Copy SS input files:
  copy_SS_inputs(dir.old = mod_dir, dir.new = aspm_dir,  create.dir = TRUE,
                 use_ss_new = TRUE, copy_par = TRUE, verbose = FALSE)
  
  # Read input files:
  starter_file = SS_readstarter(file = file.path(aspm_dir, 'starter.ss'), verbose = FALSE)
  dat_file = SS_readdat(file = file.path(aspm_dir, starter_file$datfile), verbose = FALSE)
  ctl_file = SS_readctl(file = file.path(aspm_dir, starter_file$ctlfile), datlist = dat_file, verbose = FALSE)
  par_file = SS_readpar_3.30(parfile = file.path(aspm_dir, 'ss3.par'), 
                             datsource = dat_file, ctlsource = ctl_file, verbose = FALSE)
  
  # Turn off rec devs:
  par_file$recdev1[,2] = 0
  par_file$recdev_forecast[,2] = 0
  
  # Change starter:
  starter_file$init_values_src = 1
  
  # Fix all selex parameters
  if(!is.null(ctl_file$age_selex_parms)) ctl_file$age_selex_parms[,7] = -6
  if(!is.null(ctl_file$size_selex_parms)) ctl_file$size_selex_parms[,7] = -6
  
  # Set length data comp = 0 influence
  ctl_file$lambdas = rbind( as.data.frame(expand.grid(like_comp = 4:5, fleet = 1:dat_file$Nfleets, phase = 1, value = 0, sizefreq_method = 1) ),
                            ctl_file$lambdas %>% filter(!(like_comp %in% 4:5 )) )
  
  # Fix rec devs:
  ctl_file$recdev_phase = -3
  ctl_file$recdev_early_phase = -3
  
  # Turn off recruitment likelihood:
  ctl_file$lambdas = rbind( as.data.frame(expand.grid(like_comp = 10, fleet = 0, phase = 1, value = 0, sizefreq_method = 0) ),
                            ctl_file$lambdas %>% filter(!(like_comp %in% 10 )) )
  
  # Update number lambdas:
  ctl_file$N_lambdas = nrow(ctl_file$lambdas)
  
  # Write new input files
  SS_writepar_3.30(parlist = par_file, outfile = file.path(aspm_dir, 'ss3.par'), overwrite = TRUE)
  SS_writestarter(mylist = starter_file, dir = aspm_dir, overwrite = TRUE)
  SS_writectl(ctllist = ctl_file, outfile = file.path(aspm_dir, starter_file$ctlfile), overwrite = TRUE)
  
  # Run ASPM model:
  r4ss::run(dir = aspm_dir, extras = '-nohess', exe = ss_dir)
  
}