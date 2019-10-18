library(data.table)
library(ggplot2)
library(magrittr)
spec_files = dir("~/data/spectra/", pattern = ".txt$", full.names = TRUE)
names(spec_files) = spec_files %>% basename %>% sub(" ", "_", .)  %>% sub("\\.txt", "", .)
dt = lapply(spec_files, fread) %>% rbindlist(., use.names = TRUE, idcol = "file")
ggplot(dt, aes(x = wl, y = em, color = file)) + geom_path() +
    facet_wrap(~file)


ex_dt = fread("~/data/spectra/example.emn", col.names = c("wl", "em"))
ex_dt$file = "example"
dt2 = rbind(ex_dt, dt[, .(wl, em, file)])
ggplot(dt2[file %in% c("Opal_690 _09_20", "Red_blood cell background")], aes(x = wl, y = em, color = file)) + geom_path()
p_spectra = ggplot(dt2, aes(x = wl, y = em, color = file)) + geom_path() +
    facet_wrap(~file)
p_spectra
ggsave("~/data/spectra/spectra.png", p_spectra)

for(f in spec_files){
    txt = fread(f)
    out_f = sub("\\.txt", ".emn", f)
    fwrite(txt[, .(wl = round(wl), em)], out_f, sep = "\t", col.names = FALSE)
}
