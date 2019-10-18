#!/usr/bin/env Rscript
library(magrittr)
library(data.table)
library(hsdar)
setDTthreads(threads = 1)
args = commandArgs(trailingOnly=TRUE)
out_tiff = args[1]
stitch_files = args[-1]

stitch_regions =  regmatches(basename(stitch_files), 
                             regexpr("(?<=\\.)[0-9].+[0-9](?=\\.)", 
                                     basename(stitch_files), 
                                     perl = TRUE))
reg_dt = data.table(id = stitch_regions)
reg_dt[, c("xmin", "ymin", "width", "height") := tstrsplit(id, "_")]
reg_dt$xmin = as.numeric(reg_dt$xmin)
reg_dt$ymin = as.numeric(reg_dt$ymin)
reg_dt$width = as.numeric(reg_dt$width)
reg_dt$height = as.numeric(reg_dt$height)
w = reg_dt[, max(xmin + width -1)]
h = reg_dt[, max(ymin + height -1)]
tiff_mat = matrix(0, ncol = w, nrow = h)
for(i in seq_along(stitch_files)){
    xs = seq(reg_dt$xmin[i], reg_dt$xmin[i] + reg_dt$width[i] - 1)
    ys = seq(reg_dt$ymin[i], reg_dt$ymin[i] + reg_dt$height[i] - 1)
    tiff_mat[ys, xs] = tiff::readTIFF(stitch_files[i])
}
# tiff_mat = t(tiff_mat)
# range(tiff_mat)
if(FALSE){
    plot(0, xlim = 0:1, ylim = 0:1, main="", xlab="", ylab="", axes=FALSE, mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
    rasterImage(tiff_mat, 0, 0, 1, 1)
}
tiff::writeTIFF(tiff_mat, out_tiff)
