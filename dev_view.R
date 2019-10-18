library(data.table)
library(tiff)
library(ggplot2)
# source("functions.R")
# source("function_my_processFISH.R")
tiffs = dir("~/data/TMA_Primary_0926_Runx2/", pattern = "converted_crop.+Region0001", full.names = TRUE)
names(tiffs) = sapply(strsplit(basename(tiffs), "[_\\.]"), function(x)x[7])
f = tiffs[1]
t_dt = lapply(tiffs, function(f){
    dt = as.data.table(tiff::readTIFF(f, as.is = TRUE))
    setnames(dt, c("j", "i", "z", "val"))
    dt
})
t_dt[[1]]

t_dt = lapply(t_dt, function(x){#crop
    x[i > 1500 & j > 1500]
})

# p_dt = t_dt[[1]]
t_dt = rbindlist(t_dt, use.names = TRUE, idcol = "channel")

t_dt[, norm_val := val / max(c(val, 1)), .(z, channel)]
# p1 = ggplot(t_dt, aes(x = i, y = j, fill = norm_val)) + 
#     geom_raster() + 
#     facet_grid(channel~z) +
#     theme(panel.spacing = unit(0, "npc"), panel.background = element_blank(), axis.text = element_text(size = 6)) +
#     scale_fill_gradientn(colors = c("white", "black"))
# ggsave("tmcolsp.pdf", p1, height = 20, width = 4)


agg_dt = t_dt[, .(val = sum(val)), .(channel, i, j)]
# p2 = ggplot(agg_dt, aes(x = i, y = j, fill = val)) + 
#     geom_raster() + 
#     facet_wrap(channel~.) +
#     theme(panel.spacing = unit(0, "npc"), panel.background = element_blank(), axis.text = element_text(size = 6)) +
#     scale_fill_gradientn(colors = c("white", "black")) +
#     labs(title = "raw values")
# p2

agg_dt[, norm_val := val / max(c(val, 1)), .(channel)]
# p2n = ggplot(agg_dt, aes(x = i, y = j, fill = norm_val)) + 
#     geom_raster() + 
#     facet_wrap(channel~.) +
#     theme(panel.spacing = unit(0, "npc"), panel.background = element_blank(), axis.text = element_text(size = 6)) +
#     scale_fill_gradientn(colors = c("white", "black")) +
#     labs(title = "raw values")
# p2n

#plot 3 spectra as RGB
rv = "715nm"
gv = "575nm"
bv = "455nm"
rgb_dt = dcast(agg_dt[channel %in% c(rv, gv, bv)], formula = i+j~channel, value.var = "val")
rgb_dt[[rv]] = rgb_dt[[rv]] / max(rgb_dt[[rv]])
rgb_dt[[gv]] = rgb_dt[[gv]] / max(rgb_dt[[gv]])
rgb_dt[[bv]] = rgb_dt[[bv]] / max(rgb_dt[[bv]])
rgb_dt[, chex := rgb(get(rv), get(gv), get(bv))]

lab_dt = data.table(label = c(rv, gv, bv))

cols = c("red", "green", "blue")
names(cols) = c(rv, gv, bv)
p_rgb = ggplot() + 
    geom_point(data = lab_dt, aes(x = mean(rgb_dt$i), y = mean(rgb_dt$j), color = label)) +
    geom_raster(data = rgb_dt, aes(x = i, y = j, fill = chex)) + 
    scale_color_manual(values = cols) +
    scale_y_reverse() +
    # facet_wrap(channel~.) +
    theme(panel.spacing = unit(0, "npc"), panel.background = element_blank(), axis.text = element_text(size = 6)) +
    scale_fill_identity() +
    labs(title = "raw values", x= "pixel", y = "pixel")
p_rgb

# cluster across spectra
agg_dt[, row := paste(i, j)]
agg_dt$facet = "main"
wide_dt = dcast(agg_dt, i+j~channel, value.var = "norm_val")
mat = as.matrix(wide_dt[, -1:-2])
rownames(mat) = paste(wide_dt$i, wide_dt$j)
library(seqsetvis)      
clust_dt = ssvSignalClustering(agg_dt, 
                               column_ = "channel", row_ = "row", 
                               fill_ = "norm_val", facet_ = "", 
                               max_rows = 15000, nclust = 12)
ssvSignalHeatmap(clust_dt,
                 column_ = "channel", row_ = "row", 
                 fill_ = "norm_val", facet_ = "")
aclust_dt = clust_dt[, .(y = mean(norm_val)), .(channel, cluster_id)]
aclust_dt[, .N, .(channel, cluster_id)]
ggplot(aclust_dt, aes(x = channel, y = y, group = cluster_id)) + 
    geom_path() +
    facet_wrap(cluster_id~.)

dim(mat)
set.seed(0)
pmat = mat[sample(seq_len(nrow(mat)), 2e4),]
library(BiocFileCache)
bfc = BiocFileCache()
perp = 1500
tsne_res = bfcif(bfc, paste0("dev_unmix_tsne_v1", digest::digest(list(pmat, perp))), function(){
    Rtsne::Rtsne(pmat, num_threads = 30, perplexity = perp)    
})


plot(tsne_res$Y)
tsne_dt = as.data.table(tsne_res$Y)
setnames(tsne_dt, c("tx", "ty"))

wide_dt[, id := paste(i, j)]
tsne_dt$id = rownames(pmat)
tsne_dt = merge(tsne_dt, wide_dt, by = "id")
tsne_dt[, tx := scales::rescale(tx, c(0, 1))]
tsne_dt[, ty := scales::rescale(ty, c(0, 1))]

tsne_plot = melt(tsne_dt, id.vars = c("id", "tx", "ty", "i", "j"))

#green spots
xs = c(.9, 1)
ys = c(.61, .8)
#dapi w runx
xs = c(.4, .5)
ys = c(0, .2)

p_all = ggplot(tsne_plot, 
               aes(x = tx, y = ty, color = value)) +
    facet_wrap(~variable) + 
    geom_point(size = .2) + 
    scale_color_viridis_c()+
    annotate("rect", 
             xmin = min(xs), xmax = max(xs), 
             ymin = min(ys), ymax = max(ys), 
             fill = NA, color = "red")   
# p_all + coord_cartesian(xlim = c(.75, 1), ylim = c(.48, 1))
library(seqtsne)
bins = 40
tsne_plot[, bx := bin_values(tx, bins)]
tsne_plot[, by := bin_values(ty, bins)]


ggplot(tsne_plot[, .N, .(bx, by)], aes(x = bx, y = by, fill = N)) + geom_tile() + 
    scale_fill_viridis_c(limits = c(0, 3000))

p = seqtsne::plot_binned_aggregates(
    tsne_plot, 
    xbins = bins, ybins = bins, 
    val = "value", facet_ = "variable") +
    annotate("rect", 
             xmin = min(xs), xmax = max(xs), 
             ymin = min(ys), ymax = max(ys), 
             fill = NA, color = "red")   
p
# p_all 
# p + annotate("rect", xmin = min(xs), xmax = max(xs), ymin = min(ys), ymax = max(ys), fill = NA, color = "red")    

sel_tsne = unique(tsne_plot[tx > min(xs) & tx < max(xs) & ty > min(ys) & ty < max(ys)][, .(i, j)])

cowplot::plot_grid(
    p_rgb + geom_point(data = sel_tsne, aes(x = i, y = j), shape = 1, color = "white", size = 3),
    p_rgb
)



p_rgb2 = ggplot() + 
    geom_point(data = lab_dt, aes(x = mean(rgb_dt$i), y = mean(rgb_dt$j), color = label)) +
    geom_tile(data = merge(rgb_dt, unique(tsne_plot[, .(i, j)])), aes(x = i, y = j, fill = chex)) + 
    scale_color_manual(values = cols) +
    # facet_wrap(channel~.) +
    theme(panel.spacing = unit(0, "npc"), panel.background = element_blank(), axis.text = element_text(size = 6)) +
    scale_fill_identity() +
    labs(title = "raw values used in tsne", x= "pixel", y = "pixel")

rgb_tsne_dt = merge(rgb_dt, tsne_dt[, .(i, j, tx, ty)])
ggplot(rgb_tsne_dt, aes(x = tx, y = ty, color = chex)) + geom_point() + scale_color_identity() +
    theme(panel.background = element_rect(fill = "gray20"), panel.grid = element_blank()) +
    labs(title = "rgb colorized tsne")
 