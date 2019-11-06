###setup, changing these two filepaths should be all that is required.
#set input directory that contains nd2 files
input_dir = "/slipstream/home/scope_user/data/TMA_Primary_0926_Runx2_input"
#set path to spectra file.  first column is reference spectra name ("Dapi" etc).  One column per wavelength. One row per reference.
#first column must be called "name".
spectra_csv = "ref_spec_4_no_bg.csv"
#tiffs will be digests in step x step bites, should be no need to change this.
step = 500

### Things likely to break scripts
# change to R environment on server
# anything that reorders files from ascending wavelength order.  

###Do not edit below here unless you really know what you're doing
ref_spec = data.table::fread(spectra_csv, header = TRUE)
ref_spec_dt = data.table::melt(ref_spec, value.name = "em", variable.name = "wl", id.vars = "name")

#this is trickier, path to a serialized R data object that must have wl, em, and name attributes.
#wl in spectra should match wavelengths in nd2 files.
saveRDS(ref_spec_dt, file = "ref_spec_to_use.Rds")
ref_spec_dt_rds = normalizePath("ref_spec_to_use.Rds")
# run this to see what it contains.
# readRDS(ref_spec_dt_rds) #file attribute is not used and can be removed
# ggplot(ref_spec_dt, aes(x = wl, y = em, group = name)) + geom_path() + facet_wrap(~name)

library(magrittr)
library(data.table)
library(ggplot2)
library(hsdar)
source("functions_plot_spectra.R")
source("functions_qsub.R")

stopifnot(file.exists(input_dir))
stopifnot(file.exists(ref_spec_dt_rds))

inputs = dir(input_dir, pattern = "nd2", full.names = TRUE)
nd2_file = inputs[1]
stopifnot(file.exists(nd2_file))

output_root_dir = file.path(dirname(input_dir), "unmixing3")
dir.create(output_root_dir, showWarnings = FALSE)

qsub_logging = FALSE
if(qsub_logging){
  QSUB_CMD = "qsub"  
}else{
  QSUB_CMD = "qsub -e /dev/null -o /dev/null"
}

nd2_file = inputs[1]
for(nd2_file in inputs){
  start_message = paste0("Starting (", which(nd2_file == inputs), "/", length(inputs), ") ", nd2_file)
  message(start_message)
  out_dir = file.path(output_root_dir, paste0(sub("_.+", "", basename(nd2_file)), "_unmix_v1"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  dim_str = system(paste("/slipstream/home/joeboyd/bin/showinf", nd2_file, "-nopix -no-sas | head -n 100 | awk '$0 ~ \"Width\" || $0 ~ \"Height\"'"), intern = TRUE)
  nd2_dim = as.numeric(regmatches(dim_str, regexpr(pattern = "[0-9]+", dim_str)))
  
  log_file = file.path(out_dir, "log.txt")
  write(paste0("input=", nd2_file), file = log_file)
  write(timestamp(suffix = " ---Start"), file = log_file, append = TRUE)
  
  #digest files - unmix independent
  file.copy("pipeline_digest_unmix.R", out_dir)
  dig_dir = file.path(out_dir, "digest")
  assembled_raw_dir = file.path(out_dir, "assembly_raw")
  max_dir = file.path(out_dir, "max_values")
  #unmix files
  unmix_dir = file.path(out_dir, "digest_umixed")
  
  assembled_unmix_dir = file.path(out_dir, "assembly_unmixed")
  
  
  # one_pass_unmix_dir = file.path(out_dir, "one_pass_umixed")
  
  message("Digesting large tiffs")
  write(timestamp(suffix = " ---Begin Digest"), file = log_file, append = TRUE)
  if(dir.exists(dig_dir)){
    message("directory exists! ", dig_dir, "\ndelete and rerun.")
  }else{
    dir.create(dig_dir)
    base_cmd = "/slipstream/home/joeboyd/bin/bfconvert -crop CROP_AREA IN_FILE OUT_PREFIX.Z%z_T%t_%wnm.tiff"
    
    xs = c(0, step)
    ys = c(0, step)
    
    wrap_x = TRUE
    wrap_y = TRUE
    cmds = character()
    while(wrap_x){
      if(max(xs) > nd2_dim[1]){
        xs[2] = nd2_dim[1]-1
        wrap_x = FALSE
      }
      # message(xs[1], "-", xs[2])
      while(wrap_y){
        if(max(ys) > nd2_dim[2]){
          ys[2] = nd2_dim[2]-1
          wrap_y = FALSE
        }
        # message("---", ys[1], "-", ys[2])
        crop_str = paste(min(xs), min(ys), diff(xs), diff(ys), sep = ",")
        cmd = base_cmd
        cmd = sub("CROP_AREA", crop_str, cmd)
        cmd = sub("IN_FILE", nd2_file, cmd)
        cmd = sub("OUT_PREFIX", file.path(dig_dir, paste0("raw.", gsub(",", "_", crop_str))), cmd)
        cmd_sub = paste0(QSUB_CMD, ' -v PATH=$PATH ~/../joeboyd/scripts/run_cmd.sh "', cmd, '"')
        cmds = c(cmds, cmd_sub)
        ys = ys + step
      }
      ys = c(0, step)
      wrap_y = TRUE
      xs = xs + step
    }
    qsub_and_wait(cmds)
  }
  write(timestamp(suffix = " ---End Digest"), file = log_file, append = TRUE)
  
  ### from digested tiffs to unmixed
  message("Finding global maximum")
  if(dir.exists(max_dir)){
    write(paste0("directory exists! ", max_dir, "\ndelete and rerun."), file = log_file, append = TRUE)
    message("directory exists! ", max_dir, "\ndelete and rerun.")
  }else{
    write(timestamp(suffix = " ---Begin Find Max"), file = log_file, append = TRUE)
    dir.create(max_dir)
    # roots = regBetween(dir(dig_dir, full.names = TRUE) %>% basename, reg_left = "_", reg_right = "\\.", reg_capture = "[0-9]{3}nm") %>% unique
    roots = regBetween(dir(dig_dir, full.names = TRUE) %>% basename, reg_left = "raw\\.", reg_right = "\\.", reg_capture = "[0-9_]+") %>% unique
    # roots = dir(dig_dir, full.names = TRUE) %>% basename %>% sub("\\..+", "", .) %>% unique
    
    spec_files_list = lapply(roots, function(f){
      spec_files = dir(dig_dir, pattern = paste0("^", f, "\\."), full.names = TRUE)
    })
    stopifnot(length(unique(lengths(spec_files_list))) == 1)
    file.copy("rscript_calcmax.R", file.path(max_dir, "rscript_calcmax.R"))
    
    cmds = character()
    for(i in seq_along(roots)){
      spec_files = dir(dig_dir, pattern = paste0("^raw\\.", roots[i], "\\."), full.names = TRUE)
      root = regBetween(spec_files, reg_left = "raw\\.", reg_capture = "[0-9_]+", reg_right = "\\.") %>% unique
      stopifnot(length(root) == 1)
      if(length(spec_files) < 1)stop("no spec_files")
      args = c(max_dir, root, spec_files)
      cmd = paste("Rscript rscript_calcmax.R", paste(args, collapse = " "))
      # cmd_sub = paste0('qsub -e /dev/null -o /dev/null -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
      cmd_sub = paste0(QSUB_CMD, ' -q slipstream_queue@cn-t630 -v PATH=$PATH ~/../joeboyd/scripts/run_cmd.sh "', cmd, '"')
      cmds = c(cmds, cmd_sub)
      
    }
    qsub_and_wait(cmds)
    write(timestamp(suffix = " ---End Find Max"), file = log_file, append = TRUE)
  }
  
  max_files = dir(max_dir, pattern = "_max.txt", full.names = TRUE)
  max_val = sapply(max_files, readLines) %>% as.numeric %>% max
  
  message("Unmixing")
  
  if(dir.exists(unmix_dir)){
    write(paste0("directory exists! ", unmix_dir, "\ndelete and rerun."))
    message("directory exists! ", unmix_dir, "\ndelete and rerun.")
  }else{
    
    write(timestamp(suffix = " ---Begin Unmix"), file = log_file, append = TRUE)
    dir.create(unmix_dir)
    file.copy(ref_spec_dt_rds, file.path(unmix_dir, "ummix_spectra_used.Rds"))
    # roots = dir(dig_dir, full.names = TRUE) %>% basename %>% sub("\\..+", "", .) %>% unique
    roots = regBetween(dir(dig_dir, full.names = TRUE) %>% basename, reg_left = "raw\\.", reg_right = "\\.", reg_capture = "[0-9_]+") %>% unique
    spec_files_list = lapply(roots, function(f){
      spec_files = dir(dig_dir, pattern = paste0("^raw\\.", f, "\\."), full.names = TRUE)
    })
    stopifnot(length(unique(lengths(spec_files_list))) == 1)
    cmds = character()
    for(i in seq_along(roots)){
      spec_files = dir(dig_dir, pattern = paste0("^raw\\.", roots[i], "\\."), full.names = TRUE)
      # root = sub("\\..+", "", basename(spec_files[1]))
      root = regBetween(spec_files, reg_left = "raw\\.", reg_capture = "[0-9_]+", reg_right = "\\.") %>% unique
      stopifnot(length(root) == 1)
      if(length(spec_files) < 1)stop("no spec_files")
      #NEED TO CALCULATE VAL MAX AND SUPPLY
      cmd = paste("Rscript rscript_unmix.R", paste(c(unmix_dir, root, as.character(max_val), ref_spec_dt_rds, spec_files), collapse = " "))
      cmd_sub = paste0(QSUB_CMD, ' -q slipstream_queue@cn-t630 -v PATH=$PATH ~/../joeboyd/scripts/run_cmd.sh "', cmd, '"')
      cmds = c(cmds, cmd_sub)
    }
    file.copy("rscript_unmix.R", file.path(unmix_dir, "rscript_unmix.R"))
    qsub_and_wait(cmds)
    write(timestamp(suffix = " ---End Unmix"), file = log_file, append = TRUE)
  }
  
  # if(F){
  message("Reassembling digests")
  cmds = character()
  if(dir.exists(assembled_unmix_dir)){
    write(paste0("directory exists! ", assembled_unmix_dir, "\ndelete and rerun."))
    message("directory exists! ", assembled_unmix_dir, "\ndelete and rerun.")
  }else{
    dir.create(assembled_unmix_dir)
    unmix_channels = dir(unmix_dir, full.names = TRUE, pattern = "tiff$") %>% basename %>% sub("\\..+", "", .) %>% unique
    chan = unmix_channels[2]
    for(chan in unmix_channels){
      if(chan != unmix_channels[length(unmix_channels)]){
        message(chan, ", ", appendLF = FALSE)
      }else{
        message(chan, appendLF = FALSE)
      }
      stitch_files = dir(unmix_dir, pattern = paste0(chan, "\\."), full.names = TRUE)
      out_f = file.path(assembled_unmix_dir, paste0("unmixed_", chan, ".tiff"))
      args = c(out_f, stitch_files)
      cmd = paste("Rscript rscript_restitch.R", paste(args, collapse = " "))
      cmd_sub = paste0(QSUB_CMD, ' -q slipstream_queue@cn-t630 -v PATH=$PATH ~/../joeboyd/scripts/run_cmd.sh "', cmd, '"')
      cmds = c(cmds, cmd_sub)
    }
  }
  message("Reassembling raw digests")
  if(dir.exists(assembled_raw_dir)){
    write(paste0("directory exists! ", assembled_raw_dir, "\ndelete and rerun."))
    message("directory exists! ", assembled_raw_dir, "\ndelete and rerun.")
  }else{
    dir.create(assembled_raw_dir)
    unmix_channels = dir(dig_dir, full.names = TRUE, pattern = "tiff$") %>% 
      basename %>% regmatches(., regexpr("Z.+nm", .)) %>% 
      sub("Z0_T0_", "", .) %>% unique
    
    for(chan in unmix_channels){
      if(chan != unmix_channels[length(unmix_channels)]){
        message(chan, ", ", appendLF = FALSE)
      }else{
        message(chan, appendLF = FALSE)
      }
      
      stitch_files = dir(dig_dir, pattern = paste0(chan, "\\."), full.names = TRUE)
      
      out_f = file.path(assembled_raw_dir, paste0(chan, ".tiff"))
      args = c(out_f, stitch_files)
      cmd = paste("Rscript rscript_restitch.R", paste(args, collapse = " "))
      cmd_sub = paste0(QSUB_CMD, ' -q slipstream_queue@cn-t630 -v PATH=$PATH ~/../joeboyd/scripts/run_cmd.sh "', cmd, '"')
      cmds = c(cmds, cmd_sub)
    }
  }
  
  
  if(length(cmds) == 0){
    
  }else{
    write(timestamp(suffix = " ---Begin Restitch"), file = log_file, append = TRUE)
    file.copy("rscript_restitch.R", file.path(unmix_dir, "rscript_restitch.R"))
    qsub_and_wait(cmds)
    write(timestamp(suffix = " ---End Restitch"), file = log_file, append = TRUE)
  }
  write(timestamp(suffix = " ---Stop"), file = log_file, append = TRUE)
}