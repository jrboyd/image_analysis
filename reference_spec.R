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

file2name = data.table(file = refspec_dt$file %>% unique %>% sort)
file2name$name = c("background", "MANCR", "tp63", "Runx2", "RBC", "Dapi")
refspec_dt = merge(refspec_dt, file2name, by = "file")
saveRDS(refspec_dt, file = "ref_spec_6_with_manual_RBC.Rds")
ggplot(refspec_dt, aes(x = wl, y = em)) + geom_path() + facet_wrap(~name)


