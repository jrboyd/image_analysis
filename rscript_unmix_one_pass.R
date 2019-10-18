#!/usr/bin/env Rscript
library(magrittr)
library(data.table)
library(hsdar)
setDTthreads(threads = 1)
spec_files = commandArgs(trailingOnly=TRUE)
stitch_dir = spec_files[1]
root = spec_files[2]
spec_files = spec_files[-1:-2]
# message(system("which R", intern = TRUE))
# message(paste(spec_files, collapse = "\n"))

#load ref spectra
refspec_files = dir("~/data/spectra/", pattern = ".txt$", full.names = TRUE)
names(refspec_files) = refspec_files %>% basename %>% sub(" ", "_", .)  %>% sub("\\.txt", "", .)
refspec_files= refspec_files[names(refspec_files) != "Opal_570 09_20"]
refspec_dt = lapply(refspec_files, fread) %>% rbindlist(., use.names = TRUE, idcol = "file")
refspec_dt$res = NULL

refspec_dt[, n5 := round((wl - 5) / 10) * 10 + 5]
refspec_dt = unique(refspec_dt[, .(file, wl = n5, em)])

refspec_dt
#manually override RBC bg
rbc_override = readRDS("RBC_spectra_manual_override.Rds")
refspec_dt = refspec_dt[file != "Red_blood cell background"]
refspec_dt = rbind(refspec_dt, rbc_override[, .(file, wl, em = value)])
if(FALSE){
    ggplot(refspec_dt, aes(x = wl, y = em, color = file)) + geom_path() +
        facet_wrap(~file)
}

refspec_dt.wide = dcast(refspec_dt, file~wl, value.var = "em")
file2name = data.table(file = refspec_dt.wide$file)
file2name$name = c("background", "MANCR", "tp63", "Runx2", "RBC", "Dapi")
refspec_dt.wide = merge(file2name, refspec_dt.wide, by = "file")
refspec_dt.wide$file = NULL
refspec_mat = as.matrix(refspec_dt.wide[,-1])
rownames(refspec_mat) = refspec_dt.wide$name

#load spectra tiffs
rres = regexpr("(?<=_)[0-9]{3}(?=nm)", spec_files, perl = TRUE)
spec_wl =  regmatches(spec_files, rres) %>% as.numeric()

abindlist = function(mats){
    out = mats[[1]]
    for(i in seq(2, length(mats))){
        out = abind::abind(out, mats[[i]], along = 3)
    }
    out
}

spec_mat = lapply(seq_along(spec_files), function(j){
    f = spec_files[j]
    mat = tiff::readTIFF(f, as.is = TRUE)
    mat
}) %>% abindlist()
odim = dim(spec_mat)[1:2]
dimnames(spec_mat) = list(NULL, NULL, spec_wl)
spec_mat = apply(spec_mat, 3,  function(x)x)
# spec_dt = spec_dt[, .(val = sum(val)), .(i, j, wl)]
# ggplot(spec_dt[, .(val = sum(val)), .(i, j)], aes(x = i, y = j, fill = val)) +
#     geom_raster() +
#     scale_fill_viridis_c() +
#     scale_y_reverse() +
#     theme(panel.background = element_blank()) #+ facet_wrap(~wl)

# spec_dt = dcast(spec_dt, i+j~wl, value.var = "val")
# spec_mat = as.matrix(spec_dt[, -1:-2])
# rownames(spec_mat) = paste(spec_dt$i, spec_dt$j)
#ensure reference and target spectra wavelength compatibility
stopifnot(colnames(spec_mat) == colnames(refspec_mat))
wl = as.numeric(colnames(spec_mat))

#create Speclib objects
ref_spec = new("Speclib", spectra = refspec_mat, wavelength = wl)
ref_names = rownames(refspec_mat)
obs_spec = new("Speclib", spectra = spec_mat/max(spec_mat), wavelength = wl)

#run unmix
unmix_res <- unmix(obs_spec, ref_spec)
#create unmixed images
fractions = t(unmix_res[[1]])
head(fractions)
colnames(fractions) = ref_names
errors = unmix_res[[2]]

id_dt = as.data.table(expand.grid(seq_len(odim[1]), seq_len(odim[2])))
setnames(id_dt, c("j", "i"))
#ok to here i think
unmix_dt = cbind(as.data.table(fractions), errors = errors, id_dt)
# unmix_dt = merge(unmix_dt, spec_dt[, .(id = paste(i, j), i, j)], by = "id")

# unmix_dt$id = NULL


plot_stuff = FALSE
if(plot_stuff){
    unmix_dt.tall = melt(unmix_dt, id.vars = c("i", "j"))
    unmix_dt.tall[, value_floored := value - median(value), .(variable)]
    unmix_dt.tall[value_floored < 0, value_floored := 0]
    unmix_dt = dcast(unmix_dt.tall, i+j~variable, value.var = "value_floored")
    
    plot_rgb(unmix_dt, rv = "Runx2", gv = "MANCR", bv = "Dapi", norm1 = TRUE)
    plot_rgb(unmix_dt, rv = "tp63", gv = "MANCR", bv = "Dapi", norm1 = TRUE)
    
    plot_rgb(unmix_dt, rv = "Runx2", gv = "MANCR", bv = "Dapi", norm1 = TRUE)
    
    p_unmix = ggplot(unmix_dt.tall,
                     aes(x = i, y = j, fill = value_floored)) +
        geom_raster() +
        facet_wrap(~variable) +
        scale_fill_viridis_c() +
        scale_y_reverse() +
        theme(panel.background = element_blank(), strip.background = element_blank())
}

td = colnames(unmix_dt)
td = td[! td %in% c("i", "j")]

for(nam in td){
    print(nam)
    x = unmix_dt[, .(value = get(nam), i, j)]
    x = x[order(j)][order(i)]
    arr = array(x$value, dim = c(max(x$j), max(x$i), 1))
    out_f = paste0(nam, ".", root, ".tiff")
    tiff::writeTIFF(arr/max(arr), file.path(stitch_dir, out_f))
}

# unmix_list = split(unmix_dt.tall, unmix_dt.tall$variable)
# unmix_list = lapply(names(unmix_list), function(nam){
#     x = unmix_list[[nam]]
#     x = x[order(j)][order(i)]
#     arr = array(x$value, dim = c(max(x$j), max(x$i), 1))
#     out_f = paste0(nam, ".", root, ".tiff")
#     tiff::writeTIFF(arr, file.path(stitch_dir, out_f))
# })
