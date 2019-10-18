# install.packages("hsdar")
library(hsdar)
source("functions_plot_spectra.R")
?unmix
# 
# parameter <- data.frame(LAI = seq(0, 1, 0.01))
# spectral_data <- PROSAIL(parameterList = parameter)
# 
# avl <- USGS_get_available_files()
# grass_spectra <- USGS_retrieve_files(avl = avl, pattern = "grass-fescue")
# limestone <- USGS_retrieve_files(avl = avl, pattern = "limestone")
# 
# ## Integrate all spectra to Quickbird
# grass_spectra_qb <- spectralResampling(grass_spectra[1,], "Quickbird")
# limestone_qb <- spectralResampling(limestone, "Quickbird")
# spectral_data_qb <- spectralResampling(spectral_data, "Quickbird")
# 
# 
# em <- speclib(spectra = rbind(spectra(grass_spectra_qb), 
#                               spectra(limestone_qb))/100,
#               wavelength = wavelength(limestone_qb))
# 
# ## Unmix
# unmix_res <- unmix(spectral_data_qb, em)
# 
# unmix_res
# 
# plot(unmix_res$fractions[1,] ~ SI(spectral_data_qb)$LAI, type = "l",
#      xlab = "LAI", ylab = "Unmixed fraction of vegetation")
# 
# data(spectral_data)
# plot(wavelength(spectral_data), spectra(spectral_data)[1,], type="l")

### Start
library(data.table)
library(tiff)
library(ggplot2)
library(magrittr)
spec_files = dir("~/data/spectra/", pattern = ".txt$", full.names = TRUE)
names(spec_files) = spec_files %>% basename %>% sub(" ", "_", .)  %>% sub("\\.txt", "", .)
spec_files= spec_files[names(spec_files) != "Opal_570 09_20"]
spectra_dt = lapply(spec_files, fread) %>% rbindlist(., use.names = TRUE, idcol = "file")
spectra_dt$res = NULL


# spectra_dt[, wl := round(wl)]
spectra_dt[, n5 := round((wl - 5) / 10) * 10 + 5]
spectra_dt = unique(spectra_dt[, .(file, wl = n5, em)])
ggplot(spectra_dt, aes(x = wl, y = em, color = file)) + geom_path() +
    facet_wrap(~file)
spectra_dt.wide = dcast(spectra_dt, file~wl, value.var = "em")
file2name = data.table(file = spectra_dt.wide$file)
file2name$name = c("background", "MANCR", "tp63", "Runx2", "RBC", "Dapi")
spectra_dt.wide = merge(file2name, spectra_dt.wide, by = "file")
spectra_dt.wide$file = NULL
spectra_mat = as.matrix(spectra_dt.wide[,-1])
rownames(spectra_mat) = spectra_dt.wide$name
spectra_wavelengths = as.numeric(colnames(spectra_mat))

tiffs = dir("~/data/TMA_Primary_0926_Runx2/", pattern = "converted_crop.+1000.+Region0018", full.names = TRUE)
names(tiffs) = sapply(strsplit(basename(tiffs), "[_\\.]"), function(x)x[6])

options(mc.cores = 10)
t_dt = lapply(tiffs, function(f){
    dt = as.data.table(tiff::readTIFF(f, as.is = TRUE))
    setnames(dt, c("j", "i", "z", "val"))
    dt
})
t_dt = lapply(t_dt, function(x){#crop
    x[i < 1000 & j < 1000]
})
t_dt = rbindlist(t_dt, use.names = TRUE, idcol = "channel")
t_dt = t_dt[, .(val = sum(val)), .(channel, i, j, id = paste(i, j))]


ggplot(t_dt[channel %in% c("425nm", "535nm", "575nm", "615nm", "715nm")], aes(x = i, y = j, fill = val)) + geom_raster() + facet_wrap(~channel)

t_dt.wide = dcast(t_dt, id~channel, value.var = "val")
t_mat = as.matrix(t_dt.wide[,-1])
rownames(t_mat) = t_dt.wide$id
colnames(t_mat) = sub("nm", "", colnames(t_mat))
stopifnot(colnames(t_mat) == colnames(spectra_mat))


### unmix with hsdar
?`initialize,Speclib-method`
jrb_spec = fread("jrb_ref_spectra.csv", header = TRUE)
jrb_dt = melt(jrb_spec, id.vars = c("roi_name"))


jrb_dt[roi_name == "tp63+Runx2" & as.numeric(as.character(variable)) > 645, value := value / 5]
jrb_dt[roi_name == "tp63+Runx2", roi_name := "tp63"]
jrb_dt[roi_name == "mancr+tp63-2" & 
           as.numeric(as.character(variable)) >= 595 &
           as.numeric(as.character(variable)) < 645
       , 
       value := value / 3]
jrb_dt[roi_name == "mancr+tp63-2", roi_name := "mancr"]

jrb_dt[roi_name == "background_high" & 
           as.numeric(as.character(variable)) >= 555 
       , 
       value := value / 3]


ggplot(jrb_dt, aes(x = variable, y = value, group = roi_name)) + geom_path() +
    facet_wrap(~roi_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

jrb_spec = dcast(jrb_dt, roi_name~variable, value.var = "value" )
jrb_mat = as.matrix(jrb_spec[,-1])
jrb_mat = t(apply(jrb_mat, 1, function(x)x/max(x)))
rownames(jrb_mat) = jrb_spec$roi_name

# tp63 = jrb_mat["tp63+Runx2",, drop = F]
# tp63[, as.numeric(colnames(tp63)) > 645] = tp63[, as.numeric(colnames(tp63)) > 645] / 5
# plot(tp63)
# mancr = jrb_mat["mancr+tp63-2",, drop = F]

jrb_sel_mat = jrb_mat[c(1, 2, 4, 6, 9),]



ref_spec = new("Speclib", spectra = jrb_sel_mat, wavelength = spectra_wavelengths)
ref_names = rownames(jrb_sel_mat)
# ref_spec = new("Speclib", spectra = spectra_mat, wavelength = spectra_wavelengths)
# ref_names = rownames(spectra_mat)

obs_spec = new("Speclib", spectra = t_mat/max(t_mat), wavelength = spectra_wavelengths)

unmix_res <- unmix(obs_spec, ref_spec)
fractions = t(unmix_res[[1]])
head(fractions)
colnames(fractions) = ref_names
errors = unmix_res[[2]]

unmix_dt = cbind(as.data.table(fractions), errors = errors, id = rownames(t_mat))
head(unmix_dt)
unmix_dt = melt(unmix_dt, id.vars = "id")
unmix_dt[, c("i", "j") := tstrsplit(id, " ")]
unmix_dt$i = as.numeric(unmix_dt$i)
unmix_dt$j = as.numeric(unmix_dt$j)

unmix_dt[, hit := value > median(value) *1.2, .(variable)]

unmix_dt[, rnk := rank(-value), .(variable)]
unmix_dt = unmix_dt[order(rnk)]
ggplot(unmix_dt[rnk < 15000][rnk %% 10 == 1], 
       aes(x = rnk, y = value, color = hit)) + 
    geom_path() + facet_wrap(~variable) +
    coord_cartesian(xlim = c(0, 15000))

p_unmix_binary = ggplot(unmix_dt, 
                        aes(x = i, y = j, fill = hit)) + 
    geom_raster() + 
    facet_wrap(~variable) +
    # scale_fill_viridis_c() +
    scale_y_reverse()
ggsave("unmix_binary_jrb.pdf", p_unmix_binary, width = 6, height = 5)

unmix_dt[, value_floored := value - median(value), .(variable)]
unmix_dt[value_floored < 0, value_floored := 0]
p_unmix = ggplot(unmix_dt, 
                 aes(x = i, y = j, fill = value_floored)) + 
    geom_raster() + 
    facet_wrap(~variable) +
    scale_fill_viridis_c() +
    scale_y_reverse()
ggsave("unmix_jrb.pdf", p_unmix, width = 6, height = 5)

library(data.table)
library(scales)
unmix_dt[, value_norm := scales::rescale(value_floored, 0:1), .(variable)]
rgb_dt = data.table::dcast(unmix_dt, i+j~variable, value.var = "value_norm")

rv = "Runx2"
# gv = "575nm"
gv = "mancr"
bv = "nucleus"

rgb_dt[, rgbColor := rgb(rescale(get(rv), c(0, 1)), 
                         rescale(get(gv), c(0, 1)), 
                         rescale(get(bv), c(0, 1)))]
leg_col = c("red", "green", "blue")
names(leg_col) = c(rv, gv, bv)
legend_dt = data.table::data.table(label = c(rv, gv, bv))

p_rgb = ggplot() + 
    geom_point(data = legend_dt, aes(x = mean(rgb_dt$i), y = mean(rgb_dt$j), color = label)) +
    scale_color_manual(values = leg_col) +
    geom_raster(data = rgb_dt, aes(x = i, y = j, fill = rgbColor)) + 
    scale_fill_identity() +
    scale_y_reverse() +
    theme(panel.background = element_blank()) +
    labs(title = "unmixed RGB", x= "pixel", y = "pixel")

roi = c(190, 195, 746, 754)

if(F){
    full_dt = t_dt;
    channel_dt = rgb_dt;
    xmin = 191; xmax = 195; ymin = 748; ymax = 752; view_size = 20;
    rv = "Runx2"; gv = "mancr"; bv = "nucleus"
}

plot_roi(full_dt = t_dt, 
         channel_dt = rgb_dt, 
         xmin = 190, xmax = 194, ymin = 747, ymax = 751, view_size = 20, 
         rv = "Runx2", gv = "mancr", bv = "nucleus")


plot_roi(full_dt = t_dt, 
         channel_dt = rgb_dt, 
         xmin = 716, xmax = 719, ymin = 662, ymax = 665, view_size = 20, 
         rv = "Runx2", gv = "mancr", bv = "nucleus")


# plot_roi(full_dt = t_dt, 
#          channel_dt = rgb_dt, 
#          xmin = 192, xmax = 193, ymin = 748, ymax = 750, view_size = 20, 
#          rv = "Runx2", gv = "mancr", bv = "nucleus")
# p_rgb + annotate("rect", xmin = roi[1], xmax = roi[2], ymin = roi[3], ymax = roi[4], 
#                  fill = NA, color = "white")



# coord_cartesian(xlim = roi[1:2], ylim = roi[3:4])

t_dt.roi = t_dt[i > roi_exp[1] & i <= roi_exp[2] & j > roi_exp[3] & j <= roi_exp[4]]


brks = scales::rescale(0:2/2, range(t_dt.roi$val))

p_chan.roi = ggplot(t_dt.roi, aes(x = i, y = j, fill = val)) +
    geom_raster() +
    annotate("rect", xmin = roi[1], xmax = roi[2], ymin = roi[3], ymax = roi[4],
             fill = NA, color = "blue") +
    scale_x_continuous(breaks = roi_exp[1:2]) +
    scale_y_reverse(breaks = roi_exp[3:4]) +
    scale_fill_gradientn(colours = c("white", "black"), breaks = brks) +
    facet_wrap(~channel, ncol = 8) +
    theme(panel.background = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 6),
          axis.text = element_text(size = 6), 
          legend.position = "bottom", panel.spacing = unit(0, "npc")) +
    coord_fixed() +
    labs(x = "", y = "", fill = "raw value", title = "all channels")

cowplot::plot_grid(p_rgb.roi, p_chan.roi)

unmix_dt.wide[, cols := rgb(rescale(Runx2, c(0, 1)), 
                            rescale(tp63, c(0, 1)), 
                            rescale(nucleus, c(0, 1)))]
ggplot(unmix_dt.wide, aes(x = i, y = j, fill = cols)) + 
    geom_raster() + 
    scale_fill_identity() +
    scale_y_reverse()
