#!/usr/bin/env Rscript
library(magrittr)
library(data.table)
library(hsdar)
setDTthreads(threads = 1)
args = commandArgs(trailingOnly=TRUE)
stitch_dir = args[1]
root = args[2]
spec_files = args[-1:-2]
# message(system("which R", intern = TRUE))
# message(paste(spec_files, collapse = "\n"))

#load spectra tiffs
rres = regexpr("(?<=_)[0-9]{3}(?=nm)", spec_files, perl = TRUE)
spec_wl =  regmatches(spec_files, rres) %>% as.numeric()
spec_dt = lapply(seq_along(spec_files), function(j){
    f = spec_files[j]
    dt = as.data.table(tiff::readTIFF(f, as.is = TRUE))
    setnames(dt, c("j", "i", "z", "val"))
    dt$wl = spec_wl[j]
    dt
}) %>% rbindlist()
spec_dt = spec_dt[, .(val = sum(val)), .(i, j, wl)]
max_val = max(spec_dt$val)

out_f = file.path(stitch_dir, paste0(root, "_max.txt"))
writeLines(as.character(max_val), out_f)
