library(data.table)
library(magrittr)
library(ggplot2)
source("functions_plot_spectra.R")
unmix_dir = "Region0018_unmix/assembly_unmixed/"
out_dir = "Region0018_unmix/assembly_raw/"
unmix_tiff = dir(unmix_dir, pattern = "tiff$", full.names = TRUE)
names(unmix_tiff) = regmatches(basename(unmix_tiff), regexpr("(?<=unmixed_).+(?=\\.tiff)" , basename(unmix_tiff), perl = TRUE))

names(unmix_tiff)
# raw_tiff = dir("")

ymin = 1150
ymax = 1350
xmin = 2150
xmax = 2450

ymin = 1286
ymax = 1299
xmin = 2384
xmax = 2398

rv = "Runx2"
bv = "Dapi"
gv = 'MANCR'

tiff_dt = pbapply::pblapply(unmix_tiff[c(rv, gv, bv)], function(f){
    tdat = tiff::readTIFF(f)
    dt = melt(tdat) %>% as.data.table()
    setnames(dt, c("j", "i", "value"))
    dt[i > xmin & i <= xmax & j > ymin & j <= ymax]
}) %>% rbindlist(., use.names = TRUE, idcol = "channel")


chan_dt = dcast(tiff_dt, i+j~channel, value.var = "value")

plot_rgb(channel_dt = chan_dt, rv = "Runx2", bv = "Dapi", gv = 'MANCR')


raw_tiff = dir(out_dir, pattern = "nm.tiff$", full.names = TRUE)
names(raw_tiff) = regmatches(raw_tiff, regexpr("(?<=/)[a-zA-Z0-9]+(?=nm)" , raw_tiff, perl = TRUE))
names(raw_tiff)
rawtiff_dt = pbapply::pblapply(raw_tiff, function(f){
    tdat = tiff::readTIFF(f)
    dt = melt(tdat) %>% as.data.table()
    setnames(dt, c("j", "i", "value"))
    dt = dt[i > xmin & i <= xmax & j > ymin & j <= ymax]
    dt$wl = basename(f) %>% sub("nm.tiff", "", .) %>% as.numeric()
    dt
}) %>% rbindlist(.)

basename(raw_tiff) %>% sub("nm.tiff", "", .) %>% as.numeric()
raw_tiff

ggplot(rawtiff_dt, aes(x = wl, y = value, group = paste(wl))) + geom_boxplot()

rawtiff_dt[, id := paste(j, i)]
clust_dt = seqsetvis::ssvSignalClustering(rawtiff_dt, column_ = "wl", facet_ = "", fill_ = "value")
agg_dt = clust_dt[, ]

rbc_dt = rawtiff_dt[, .(value = median(value)), .(wl)]
rbc_dt[, value := value / max(value)]
ggplot(rbc_dt, aes(x = wl, y = value)) + geom_path()
rbc_dt$file = "RBC"
saveRDS(rbc_dt, "RBC_spectra_manual_override.Rds")
