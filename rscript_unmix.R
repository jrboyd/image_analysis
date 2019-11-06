#!/usr/bin/env Rscript
library(magrittr)
library(data.table)
library(hsdar)
setDTthreads(threads = 1)
args = commandArgs(trailingOnly=TRUE)
stitch_dir = args[1]
root = args[2]
val_max = as.numeric(args[3])
refspec_dt = readRDS(args[4])
spec_files = args[-1:-4]
# message(system("which R", intern = TRUE))
# message(paste(spec_files, collapse = "\n"))

#load ref spectra
refspec_dt.wide = dcast(refspec_dt, name~wl, value.var = "em")

refspec_mat = as.matrix(refspec_dt.wide[,-1])
rownames(refspec_mat) = refspec_dt.wide$name
ggplot(refspec_dt, aes(x = wl, y = em)) + geom_path() + facet_wrap(~name)

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
# ggplot(spec_dt[, .(val = sum(val)), .(i, j)], aes(x = i, y = j, fill = val)) +
#     geom_raster() +
#     scale_fill_viridis_c() +
#     scale_y_reverse() +
#     theme(panel.background = element_blank()) #+ facet_wrap(~wl)

spec_dt = dcast(spec_dt, i+j~wl, value.var = "val")
spec_mat = as.matrix(spec_dt[, -1:-2])
rownames(spec_mat) = paste(spec_dt$i, spec_dt$j)
#ensure reference and target spectra wavelength compatibility
stopifnot(colnames(spec_mat) == colnames(refspec_mat))
wl = as.numeric(colnames(spec_mat))

#create Speclib objects
ref_spec = new("Speclib", spectra = refspec_mat, wavelength = wl)
ref_names = rownames(refspec_mat)
obs_spec = new("Speclib", spectra = spec_mat/val_max, wavelength = wl)

#run unmix
unmix_res <- unmix(obs_spec, ref_spec)
#create unmixed images
fractions = t(unmix_res[[1]])
head(fractions)
colnames(fractions) = ref_names
errors = unmix_res[[2]]
unmix_dt = cbind(as.data.table(fractions), errors = errors, id = rownames(spec_mat))
unmix_dt = merge(unmix_dt, spec_dt[, .(id = paste(i, j), i, j)], by = "id")

unmix_dt$id = NULL
unmix_dt.tall = melt(unmix_dt, id.vars = c("i", "j"))

plot_stuff = FALSE
if(plot_stuff){
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

unmix_list = split(unmix_dt.tall, unmix_dt.tall$variable)
unmix_list = lapply(names(unmix_list), function(nam){
    x = unmix_list[[nam]]
    x = x[order(j)][order(i)]
    arr = array(x$value, dim = c(max(x$j), max(x$i), 1))
    out_f = paste0(nam, ".", root, ".tiff")
    tiff::writeTIFF(arr, file.path(stitch_dir, out_f))
})
