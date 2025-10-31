create_retro_files <- function(dir = getwd(), 
                               oldsubdir = "", newsubdir = "retrospectives",
                               subdirstart = "retro", years = 0:-5, overwrite = TRUE,
                               verbose = FALSE, mod_ts = 'year') {
  
  olddir <- file.path(dir, oldsubdir)
  newdir <- file.path(dir, newsubdir)
  if(mod_ts == 'year') tstep = 1
  if(mod_ts == 'quarter') tstep = 4
  
  # get model file names from olddir
  startfile <- dir(olddir)[tolower(dir(olddir)) == "starter.ss"]
  if (length(startfile) == 0) {
    stop("No starter.ss file found in ", olddir)
  }
  
  # read original starter (later written to each folder)
  startfile <- file.path(olddir, startfile)
  starter <- SS_readstarter(startfile, verbose = FALSE)
  subdirnames <- paste0(subdirstart, years)
  
  save_dirs = NULL
  # loop over retrospective years
  for (iyr in seq_along(years)) {
    newdir_iyr <- file.path(newdir, subdirnames[iyr])
    message("Creating retrospective files in ", newdir_iyr)
    
    # copy original input files to retro folder
    copy_SS_inputs(
      dir.old = olddir,
      dir.new = newdir_iyr,
      create.dir = TRUE,
      recursive = TRUE,
      overwrite = TRUE,
      copy_exe = FALSE,
      verbose = verbose
    )
    
    # change starter file to do retrospectives
    starter[["retro_yr"]] <- years[iyr]*tstep
    starter[["init_values_src"]] <- 0
    SS_writestarter(starter,
                    dir = newdir_iyr,
                    verbose = FALSE,
                    overwrite = TRUE
    )
    
    # delete covar file to avoid using file from previous model run
    # (not sure if this is necessary)
    if (file.exists("covar.sso")) {
      file.remove("covar.sso")
    }
    
    save_dirs = c(save_dirs, newdir_iyr)
    
  }
  
  return(save_dirs)
  
}