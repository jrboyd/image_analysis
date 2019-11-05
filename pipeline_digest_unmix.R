library(magrittr)
library(data.table)
library(ggplot2)
library(hsdar)
source("functions_plot_spectra.R")
source("functions_qsub.R")
### Things likely to break scripts
# change to R environment on server
# anything that reorders files from ascending wavelength
# change in dimensions of images

### from nd2 to subsections of tiffs for each channel
# nd2_dim = c(5700, 5700)
step = 500

input_dir = "/slipstream/home/joeboyd/data/TMA_Primary_0926_Runx2_input"
stopifnot(file.exists(input_dir))
ref_spec_dt_rds = file.path(getwd(), "ref_spec_4_no_bg.Rds")
stopifnot(file.exists(ref_spec_dt_rds))


inputs = dir(input_dir, pattern = "nd2", full.names = TRUE)
nd2_file = inputs[1]
stopifnot(file.exists(nd2_file))
for(nd2_file in inputs){
    out_dir = file.path(dirname(input_dir), "unmixing", paste0(sub("_.+", "", basename(nd2_file)), "_unmix_v1"))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    dim_str = system(paste("/slipstream/home/joeboyd/bin/showinf", nd2_file, "-nopix -no-sas | head -n 100 | awk '$0 ~ \"Width\" || $0 ~ \"Height\"'"), intern = TRUE)
    nd2_dim = as.numeric(regmatches(dim_str, regexpr(pattern = "[0-9]+", dim_str)))
    
    log_file = file.path(out_dir, "log.txt")
    write(timestamp(suffix = " ---Start"), file = log_file)
    
    #digest files - unmix independent
    file.copy("pipeline_digest_unmix.R", out_dir)
    dig_dir = file.path(out_dir, "digest")
    assembled_raw_dir = file.path(out_dir, "assembly_raw")
    max_dir = file.path(out_dir, "max_values")
    #unmix files
    unmix_dir = file.path(out_dir, "unmixed_signal4")
    file.copy(ref_spec_dt_rds, file.path(unmix_dir, "ummix_spectra_used.Rds"))
    assembled_unmix_dir = file.path(out_dir, "assembly_signal4")
  
    
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
                cmd = sub("OUT_PREFIX", file.path(dig_dir, gsub(",", "_", crop_str)), cmd)
                cmd_sub = paste0('qsub -e /dev/null -o /dev/null -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
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
    
    
    
    #load tiffs
    
    
    
    
    
    message("Finding global maximum")

    if(dir.exists(max_dir)){
        write(paste0("directory exists! ", max_dir, "\ndelete and rerun."), file = log_file, append = TRUE)
        message("directory exists! ", max_dir, "\ndelete and rerun.")
    }else{
        write(timestamp(suffix = " ---Begin Find Max"), file = log_file, append = TRUE)
        dir.create(max_dir)
        roots = dir(dig_dir, full.names = TRUE) %>% basename %>% sub("\\..+", "", .) %>% unique
        
        spec_files_list = lapply(roots, function(f){
            spec_files = dir(dig_dir, pattern = paste0("^", f, "\\."), full.names = TRUE)
        })
        stopifnot(length(unique(lengths(spec_files_list))) == 1)
        file.copy("rscript_calcmax.R", file.path(max_dir, "rscript_calcmax.R"))
        
        cmds = character()
        for(i in seq_along(roots)){
            spec_files = dir(dig_dir, pattern = paste0("^", roots[i], "\\."), full.names = TRUE)
            root = sub("\\..+", "", basename(spec_files[1]))
            if(length(spec_files) < 1)stop("no spec_files")
            args = c(max_dir, root, spec_files)
            cmd = paste("Rscript rscript_calcmax.R", paste(args, collapse = " "))
            # cmd_sub = paste0('qsub -e /dev/null -o /dev/null -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
            cmd_sub = paste0('qsub -e /dev/null -o /dev/null -q slipstream_queue@cn-t630 -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
            cmds = c(cmds, cmd_sub)
            # cmd_sub = paste0('qsub -e /dev/null -o /dev/null -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
            
        }
        qsub_and_wait(
            
        )
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
        roots = dir(dig_dir, full.names = TRUE) %>% basename %>% sub("\\..+", "", .) %>% unique
        
        spec_files_list = lapply(roots, function(f){
            spec_files = dir(dig_dir, pattern = paste0("^", f, "\\."), full.names = TRUE)
        })
        stopifnot(length(unique(lengths(spec_files_list))) == 1)
        cmds = character()
        for(i in seq_along(roots)){
            spec_files = dir(dig_dir, pattern = paste0("^", roots[i], "\\."), full.names = TRUE)
            root = sub("\\..+", "", basename(spec_files[1]))
            if(length(spec_files) < 1)stop("no spec_files")
            #NEED TO CALCULATE VAL MAX AND SUPPLY
            cmd = paste("Rscript rscript_unmix.R", paste(c(unmix_dir, root, as.character(max_val), ref_spec_dt_rds, spec_files), collapse = " "))
            cmd_sub = paste0('qsub -e /dev/null -o /dev/null -q slipstream_queue@cn-t630 -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
            cmds = c(cmds, cmd_sub)
        }
        file.copy("rscript_unmix.R", file.path(unmix_dir, "rscript_unmix.R"))
        qsub_and_wait(cmds)
        write(timestamp(suffix = " ---End Unmix"), file = log_file, append = TRUE)
    }
    
    # if(F){
    message("Reassembling digests")
    cmds = character()
    if(F){#dir.exists(assembled_unmix_dir)){
        write(paste0("directory exists! ", assembled_unmix_dir, "\ndelete and rerun."))
        message("directory exists! ", assembled_unmix_dir, "\ndelete and rerun.")
    }else{
        dir.create(assembled_unmix_dir)
        unmix_channels = dir(unmix_dir, full.names = TRUE, pattern = "tiff$") %>% basename %>% sub("\\..+", "", .) %>% unique
        chan = unmix_channels[2]
        for(chan in unmix_channels){
            message(chan)
            stitch_files = dir(unmix_dir, pattern = paste0(chan, "\\."), full.names = TRUE)
            out_f = file.path(assembled_unmix_dir, paste0("unmixed_", chan, ".tiff"))
            args = c(out_f, stitch_files)
            cmd = paste("Rscript rscript_restitch.R", paste(args, collapse = " "))
            cmd_sub = paste0('qsub -e /dev/null -o /dev/null -q slipstream_queue@cn-t630 -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
            cmds = c(cmds, cmd_sub)
            
            # stitch_regions =  regmatches(basename(stitch_files), 
            #                              regexpr("(?<=\\.)[0-9].+[0-9](?=\\.)", 
            #                                      basename(stitch_files), 
            #                                      perl = TRUE))
            # reg_dt = data.table(id = stitch_regions)
            # reg_dt[, c("xmin", "ymin", "width", "height") := tstrsplit(id, "_")]
            # reg_dt$xmin = as.numeric(reg_dt$xmin)
            # reg_dt$ymin = as.numeric(reg_dt$ymin)
            # reg_dt$width = as.numeric(reg_dt$width)
            # reg_dt$height = as.numeric(reg_dt$height)
            # w = reg_dt[, max(xmin + width -1)]
            # h = reg_dt[, max(ymin + height -1)]
            # tiff_mat = matrix(0, ncol = w, nrow = h)
            # for(i in seq_along(stitch_files)){
            #     xs = seq(reg_dt$xmin[i], reg_dt$xmin[i] + reg_dt$width[i] - 1)
            #     ys = seq(reg_dt$ymin[i], reg_dt$ymin[i] + reg_dt$height[i] - 1)
            #     tiff_mat[ys, xs] = tiff::readTIFF(stitch_files[i])
            # }
            # # tiff_mat = t(tiff_mat)
            # # range(tiff_mat)
            # if(FALSE){
            #     plot(0, xlim = 0:1, ylim = 0:1, main="", xlab="", ylab="", axes=FALSE, mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
            #     rasterImage(tiff_mat, 0, 0, 1, 1)
            # }
            # tiff::writeTIFF(tiff_mat, file.path(assembled_unmix_dir, paste0("unmixed_", chan, ".tiff")))
        }
    }
    # }
    
    # if(F){
    message("Reassembling raw digests")
    if(F){#dir.exists(assembled_raw_dir)){
        write(paste0("directory exists! ", assembled_raw_dir, "\ndelete and rerun."))
        message("directory exists! ", assembled_raw_dir, "\ndelete and rerun.")
    }else{
        dir.create(assembled_raw_dir)
        unmix_channels = dir(dig_dir, full.names = TRUE, pattern = "tiff$") %>% 
            basename %>% regmatches(., regexpr("Z.+nm", .)) %>% 
            sub("Z0_T0_", "", .) %>% unique
        
        for(chan in unmix_channels){
            message(chan)
            stitch_files = dir(dig_dir, pattern = paste0(chan, "\\."), full.names = TRUE)
            
            out_f = file.path(assembled_raw_dir, paste0("unmixed_", chan, ".tiff"))
            args = c(out_f, stitch_files)
            cmd = paste("Rscript rscript_restitch.R", paste(args, collapse = " "))
            cmd_sub = paste0('qsub -e /dev/null -o /dev/null -q slipstream_queue@cn-t630 -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
            cmds = c(cmds, cmd_sub)
            # stitch_regions =  basename(stitch_files) %>%
            #     sub("\\..+", "", .)
            # reg_dt = data.table(id = stitch_regions)
            # reg_dt[, c("xmin", "ymin", "width", "height") := tstrsplit(id, "_")]
            # reg_dt$xmin = as.numeric(reg_dt$xmin)
            # reg_dt$ymin = as.numeric(reg_dt$ymin)
            # reg_dt$width = as.numeric(reg_dt$width)
            # reg_dt$height = as.numeric(reg_dt$height)
            # w = reg_dt[, max(xmin + width -1)]
            # h = reg_dt[, max(ymin + height -1)]
            # tiff_mat = matrix(0, ncol = w, nrow = h)
            # for(i in seq_along(stitch_files)){
            #     xs = seq(reg_dt$xmin[i], reg_dt$xmin[i] + reg_dt$width[i] - 1)
            #     ys = seq(reg_dt$ymin[i], reg_dt$ymin[i] + reg_dt$height[i] - 1)
            #     tiff3chan = tiff::readTIFF(stitch_files[i])
            #     
            #     tiff_mat[ys, xs] = apply(tiff3chan, 1:2, sum)
            # }
            # # tiff_mat = t(tiff_mat)
            # # range(tiff_mat)
            # if(FALSE){
            #     plot(0, xlim = 0:1, ylim = 0:1, main="", xlab="", ylab="", axes=FALSE, mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
            #     rasterImage(tiff_mat/max(tiff_mat), 0, 0, 1, 1)
            # }
            # tiff::writeTIFF(tiff_mat/max(tiff_mat), file.path(assembled_raw_dir, paste0("raw_", chan, ".tiff")))
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
# if(dir.exists(one_pass_unmix_dir)){
#     message("directory exists! ", one_pass_unmix_dir, "\ndelete and rerun.")
# }else{
#     dir.create(one_pass_unmix_dir)
#     spec_files = dir(assembled_raw_dir, full.names = TRUE)
#     root = "one_pass"
#     if(length(spec_files) < 1)stop("no spec_files")
#     # spec_files = c(unmix_dir, root, spec_files)
#     cmd = paste("Rscript rscript_unmix.R", paste(c(unmix_dir, root, spec_files), collapse = " "))
#     # cmd_sub = paste0('qsub -e /dev/null -o /dev/null -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
#     cmd_sub = paste0('qsub -q slipstream_queue@cn-t630 -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
#     # cmd_sub = paste0('qsub -e /dev/null -o /dev/null -v PATH=$PATH ~/scripts/run_cmd.sh "', cmd, '"')
#     system(cmd_sub)
# }
# }