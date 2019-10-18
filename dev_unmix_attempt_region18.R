library(data.table)
library(tiff)
library(ggplot2)
library(magrittr)
spec_files = dir("~/data/spectra/", pattern = ".txt$", full.names = TRUE)
names(spec_files) = spec_files %>% basename %>% sub(" ", "_", .)  %>% sub("\\.txt", "", .)
spec_files= spec_files[names(spec_files) != "Opal_570 09_20"]
spectra_dt = lapply(spec_files, fread) %>% rbindlist(., use.names = TRUE, idcol = "file")
spectra_dt$res = NULL
ggplot(spectra_dt, aes(x = wl, y = em, color = file)) + geom_path() +
    facet_wrap(~file)
bg_dt = spectra_dt[file == file[1]]
bg_dt$file = "background"
bg_dt$em = 0

spectra_dt.bg = rbind(spectra_dt, bg_dt)


# ggplot(spectra_dt[file %in% c("Opal_690 _09_20", "Red_blood cell background")], aes(x = wl, y = em, color = file)) + geom_path()
p_spectra = ggplot(spectra_dt, aes(x = wl, y = em, color = file)) + geom_path() +
    facet_wrap(~file)
p_spectra

tiffs = dir("~/data/TMA_Primary_0926_Runx2/", pattern = "converted_crop.+1000.+Region0018", full.names = TRUE)
names(tiffs) = sapply(strsplit(basename(tiffs), "[_\\.]"), function(x)x[6])
f = tiffs[1]
t_dt = lapply(tiffs, function(f){
    dt = as.data.table(tiff::readTIFF(f, as.is = TRUE))
    setnames(dt, c("j", "i", "z", "val"))
    dt
})
# t_dt[[1]]

# t_dt = lapply(t_dt, function(x){#crop
#     x[i > 1500 & j > 1500]
# })
t_dt = rbindlist(t_dt, use.names = TRUE, idcol = "channel")
# t_dt[]

t_dt[, norm_val := val / max(c(val, 1)), .(z, channel)]
agg_dt = t_dt[, .(val = sum(val)), .(channel, i, j)]
agg_dt[, norm_val := val / max(c(val, 1)), .(channel)]

remove(t_dt)
#plot 3 spectra as RGB
# rv = "715nm"
# rv = "575nm"
# rv = "625nm"
rv = "575nm"
# gv = "575nm"
gv = "625nm"
bv = "435nm"

unique(agg_dt$channel)
p_spectra + 
    annotate("line", x = as.numeric(sub("nm", "", bv)), y = c(0, 1), color = "blue", size = 2) +
    annotate("line", x = as.numeric(sub("nm", "", gv)), y = c(0, 1), color = "green", size = 2) +
    annotate("line", x = as.numeric(sub("nm", "", rv)), y = c(0, 1), color = "red", size = 2)

rgb_dt = dcast(agg_dt[channel %in% c(rv, gv, bv)], formula = i+j~channel, value.var = "val")
# rgb_dt[[rv]] = rgb_dt[[rv]] / max(rgb_dt[[rv]])
# rgb_dt[[gv]] = rgb_dt[[gv]] / max(rgb_dt[[gv]])
# rgb_dt[[bv]] = rgb_dt[[bv]] / max(rgb_dt[[bv]])
quantile(rgb_dt[[rv]], 380:400/400)
rgb_dt[[rv]] = rgb_dt[[rv]] / quantile(rgb_dt[[rv]], .995)
rgb_dt[[gv]] = rgb_dt[[gv]] / quantile(rgb_dt[[gv]], .995)
rgb_dt[[bv]] = rgb_dt[[bv]] / quantile(rgb_dt[[bv]], .995)

set(rgb_dt, i = which(rgb_dt[[rv]] > 1), j = rv, value = 1)
set(rgb_dt, i = which(rgb_dt[[gv]] > 1), j = gv, value = 1)
set(rgb_dt, i = which(rgb_dt[[bv]] > 1), j = bv, value = 1)

rgb_dt[, chex := rgb(get(rv), get(gv), get(bv))]



legend_dt = data.table(label = c(rv, gv, bv))

cols = c("red", "green", "blue")
names(cols) = c(rv, gv, bv)
p_rgb = ggplot() + 
    geom_point(data = legend_dt, aes(x = mean(rgb_dt$i), y = mean(rgb_dt$j), color = label)) +
    geom_raster(data = rgb_dt, aes(x = i, y = j, fill = chex)) + 
    scale_color_manual(values = cols) +
    scale_y_reverse() +
    # facet_wrap(channel~.) +
    theme(panel.spacing = unit(0, "npc"), panel.background = element_blank(), axis.text = element_text(size = 6)) +
    scale_fill_identity() +
    labs(title = "raw values", x= "pixel", y = "pixel")
p_rgb

labs = list("background_high" = c(287, 295, 313, 325),
            "background_low" = c(480, 510, 842, 872),
            "nucleus" = c(483, 500, 484, 495),
            "nucleus+Runx2" = c(473, 480, 486, 497),
            "nucleus" = c(530, 545, 495, 510),
            "mancr+tp63-1" = c(857, 860, 661, 663),
            "Runx2" = c(955, 970, 545, 551),
            "mancr+tp63-2" = c(717, 720, 663, 666),
            "tp63+Runx2" = c(86, 89, 833, 836),
            "nucleus+mancr+Runx2" = c(74, 76, 27, 30))
lab_dt = lapply(labs, matrix, nrow = 1) %>% lapply(., as.data.frame) %>%
    rbindlist(., use.names = TRUE, idcol = "label") %>% setnames(., old = paste0("V", 1:4), c("xmin", "xmax", "ymin", "ymax"))

# xs = c(1684, 1690)-10
# ys = c(1980, 1987)-25
# xs = c(1800, 1850)-30
# ys = c(200, 270)-25
p_rgb.ann = p_rgb + geom_rect(data = lab_dt, 
                              mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                              fill= NA,
                              color = "white")
# p_rgb.ann
ggsave("ROI_global.png", p_rgb.ann)

pdf("ROI_details.pdf", width = 9, height = 6)

sel_dts = list()
for(i in seq_len(nrow(lab_dt))){
    
    d = lab_dt[i,]
    xlim = c(d$xmin, d$xmax)
    ylim = c(d$ymin, d$ymax)
    message(d$label)
    
    
    pad = max(diff(xlim), diff(ylim))/2
    
    sel_agg_dt_padded = agg_dt[i >= d$xmin-pad & i <= d$xmax+pad & j >= d$ymin-pad & j <= d$ymax+pad]
    p1 = ggplot() +
        geom_raster(data = sel_agg_dt_padded, aes(x = i, y = j, fill = val)) + 
        geom_rect(data = d, 
                  mapping = aes(xmin = xmin-.5, xmax = xmax-.5, ymin = ymin-.5, ymax = ymax-.5),
                  fill= NA,
                  color = "blue") +
        facet_wrap(~channel) +
        scale_fill_gradientn(colors = c("white", "black")) +
        scale_x_continuous(breaks = range(sel_agg_dt_padded$i)) +
        scale_y_reverse(breaks = range(sel_agg_dt_padded$j)) +
        theme(panel.background = element_blank(), 
              strip.background = element_blank(), 
              strip.text = element_text(size = 6),
              axis.text = element_text(size = 6), 
              legend.position = "bottom") +
        labs(x = "", y = "",  title = d$label) +
        coord_fixed()
    p1
    
    # sel__rgb_dt = rgb_dt[i >= d$xmin-pad & i <= d$xmax+pad & j >= d$ymin-pad & j <= d$ymax+pad]
    # p1 = ggplot() + 
    #     geom_point(data = legend_dt, aes(x = mean(sel__rgb_dt$i), y = mean(sel__rgb_dt$j), color = label)) +
    #     # geom_tile(data = sel__rgb_dt[`575nm` > .6 & `535nm` > .6 & `435nm` < .5],
    #     geom_tile(data = sel__rgb_dt, 
    #               aes(x = i, y = j, fill = chex)) + 
    #     geom_rect(data = d, 
    #               mapping = aes(xmin = xmin-.5, xmax = xmax-.5, ymin = ymin-.5, ymax = ymax-.5),
    #               fill= NA,
    #               color = "white") +
    #     scale_color_manual(values = cols) +
    #     scale_y_reverse() +
    #     # facet_wrap(channel~.) +
    #     theme(panel.spacing = unit(0, "npc"), panel.background = element_blank(), axis.text = element_text(size = 6)) +
    #     scale_fill_identity() +
    #     labs(title = "raw values", x= "pixel", y = "pixel") + 
    #     # coord_cartesian(xlim = xlim, ylim = ylim) +
    #     labs(title = d$label) +
    #     coord_fixed()
    # 
    
    pad = 0
    sel_agg_dt = agg_dt[i >= d$xmin-pad & i < d$xmax+pad & j >= d$ymin-pad & j < d$ymax+pad]
    p2 = ggplot(sel_agg_dt, 
                aes(x = channel, y = val)) +
        geom_boxplot() +
        theme(axis.text = element_text(angle = 45, hjust = 1, size = 6))
    sel_agg_dt[, id := paste(i, j)]
    if(length(unique(sel_agg_dt$id)) < 10){
        p3 = ggplot(sel_agg_dt, aes(x = wl, y = val, group = id)) + 
            geom_path() +
            stat_summary(geom = "path", mapping = aes(group = NA), color = "red", size = 1.5)
    }else{
        p3 = seqsetvis::ssvSignalHeatmap(sel_agg_dt, row_ = "id", column_ = "channel", fill_ = "val", facet_ = "") +
            scale_fill_viridis_c() +
            theme(axis.text.x = element_text(size = 6, hjust = 1, angle = 45))    
    }
    
    # cowplot::plot_grid(cowplot::plot_grid(p1, p2, ncol = 1),  p3) %>% show
    sel_dts[[d$label]] = sel_agg_dt
    cowplot::plot_grid(p1, cowplot::plot_grid(p2, p3, ncol = 1), 
                       rel_widths = c(1.3, 1), nrow = 1) %>% show
    
    # p1 + geom_tile(data = sel_agg_dt[channel == "425nm"], mapping = aes(x = i, y = j), fill = "red", alpha = .7)
    
}
dev.off()

roi_dt = rbindlist(sel_dts, use.names = TRUE, idcol = "roi_name")
unique(roi_dt[, .(id, roi_name)])$roi_name %>% table

# tk = roi_dt[, .(maxVal = max(norm_val)), .(id, roi_name)][maxVal > median(maxVal), .(id, roi_name)]
# roi_dt = merge(roi_dt, tk)
roi_dt[, refVal := sum(val), .(roi_name, id)]
roi_dt[, rnk := rank(-refVal), .(roi_name, wl)]
roi_dt[roi_name == "mancr+tp63-1"]
roi_dt[roi_name == "background_low"]$rnk %>% table
roi_dt = roi_dt[rnk < 10]

median_dt = roi_dt[, .(val = mean(val)), .(wl, roi_name)]
ggplot(roi_dt, aes(x = wl, y = val, group = id)) + 
    geom_path() +
    geom_path(data = median_dt, aes(group = NA), color = "red", size = 1.5) +
    facet_wrap(~roi_name, scales = "free_y")
    # stat_summary(geom = "path", mapping = aes(group = NA), color = "red", size = 1.5, fun.y = median) +
    # labs(title = roi_name) 

# split(median_dt, by = median_dt$roi_name)
fwrite(dcast(median_dt, roi_name~wl, value.var = "val"), "jrb_ref_spectra.csv")

for(i in seq_along(sel_dts)){
    roi_name = names(sel_dts)[i]
    dt = sel_dts[[i]]
    hist(dt[, max(norm_val), .(id)]$V1)
    
    

}


# p_rgb.ann + scale_y_reverse(limits = ylim) + scale_x_continuous(limits = xlim)
# (xlim = xs, ylim = ys)
#sel here
sel_dt = agg_dt[i >= min(xs) & i <= max(xs) & j >= min(ys) & j <= max(ys)]
sel_dt[, wl := as.numeric(sub("nm", "", channel))]
ggplot(sel_dt, aes(x = wl, y = val, group = wl)) + geom_boxplot()



agg_dt[, wl := as.numeric(sub("nm", "", channel))]


ggplot(agg_dt[i >= min(xs) & i <= max(xs) & 
                  j >= min(ys) & j <= max(ys)], 
       aes(x = i, y = j, fill = val)) + 
    geom_raster() +
    facet_wrap(~channel)

ggplot(agg_dt[i >= min(xs) & i <= max(xs) & 
                  j >= min(ys) & j <= max(ys)],
       aes(x = wl, y = val, group = paste(i, j))) +geom_path() +
    facet_wrap(~i+j)

#create kmeans centroids
spec_ref_wl = spectra_dt.bg$wl %>% unique
obs_wl = agg_dt$wl %>% unique
nearest_ref = sapply(obs_wl, function(x){
    o = order(abs(x - spec_ref_wl))
    spec_ref_wl[o][1]
})
spec_wide = dcast(spectra_dt.bg, file~wl, value.var = "em")
spec_wide = spec_wide[, c("file", as.character(nearest_ref)), with = FALSE]
spec_mat = as.matrix(spec_wide[, -1])
rownames(spec_mat) = spec_wide$file

agg_dt[, px_val := val / max(val), .(i, j)]
agg_dt[, id := paste(i, j)]
ggplot(agg_dt[id %in% sample(unique(id), 5)], aes(x = wl, y = px_val)) + geom_path() + facet_wrap(~id)

obs_dt = dcast(agg_dt, id~wl, value.var = "px_val")
obs_mat = as.matrix(obs_dt[,-1])
rownames(obs_mat) = obs_dt$id

dist(obs_mat)
spec_mat

tmp = t(apply(obs_mat, 1, function(x){
    sum(abs(x - spec_mat))
}))
dim(tmp)

exikm = kmeans(obs_mat, centers = spec_mat)
km_dt = data.table(cluster_id = km$cluster, id = names(km$cluster))
obs_dt = merge(obs_dt, km_dt, by = "id")
tall_dt = melt(obs_dt, id.vars = c("id", "cluster_id"))
tall_dt.mean = tall_dt[, .(value = mean(value)), .(cluster_id, variable)]
ggplot(tall_dt.mean, aes(x = variable, y = value, group = cluster_id)) + facet_grid(cluster_id~.) + geom_path()

s1 = spec_mat[,1]
s2 = spec_mat[,2]
s3 = spec_mat[,3]
group = rownames(spec_mat)
lm.s3 = lm(s1 + s2 + s3 ~ group - 1)
summary(lm.s3)
## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
lm.D90 <- lm(weight ~ group - 1) # omitting intercept

anova(lm.D9)
summary(lm.D90)

opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.D9, las = 1)      # Residuals, Fitted, ...
par(opar)

spec_wide
spec_mat
obs_mat
cor_mat = cor(t((obs_mat)), t(spec_mat))
cor_dt = as.data.table(cor_mat, keep.rownames = TRUE)
setnames(cor_dt, "rn", "id")
cor_dt$background = NULL
cor_dt[, c("i", "j") := tstrsplit(id, " ")]
cor_dt$i = as.numeric(cor_dt$i)
cor_dt$j = as.numeric(cor_dt$j)
cor_dt = melt(cor_dt, id.vars = c("id", "i", 'j'))
pg = cowplot::plot_grid(p_rgb, 
                        ggplot(cor_dt, aes(x = i, y = j, fill = value)) + 
                            facet_wrap(~variable) +
                            geom_raster() + 
                            scale_y_reverse() +
                            scale_fill_viridis_c(limits = c(.5, NA))
)
ggsave("tmp.png", pg, width = 8, height = 4, units = "in")
# p_rgb
# ggplot(cor_dt, aes(x = i, y = j, fill = value)) + 
#     facet_wrap(~variable) +
#     geom_raster() + 
#     scale_y_reverse() +
#     scale_fill_viridis_c()

