setwd("~/R/PG_images/")

library(FISHalyseR)
library(MaxContrastProjection)
source("functions.R")

# f = system.file( "extdata", "SampleFISHgray.jpg", package="FISHalyseR")
root = file.path("~/knime-workspace/NCI_imaging/", 
                 "GLP-181017-Prachi-FISH-60X_20181017_121249/AssayPlate_PerkinElmer_CellCarrier-96 Ultra/")

files = dir(root, pattern = "Assay.+tif", full.names = TRUE)
mat = matrix(unlist(strsplit(files, "_")), ncol = 5, byrow = T)

#A01 is DAPI
#A02 is ?
#A03 is ?
#F are views?
#G are bio reps?



# nucleus
analyze_params(name = "flatten_nucleus", key = "G08.+F001L01A01Z[0-9].+.tif", 
               min_size = 1000, max_size = 25000, w = 150)
analyze_params(name = "flatten_probeA", key = "G08.+F001L01A02Z[0-9].+.tif", 
               min_size = 1, max_size = 100, w = 20)
analyze_params(name = "flatten_probeB", key = "G08.+F001L01A03Z[0-9].+.tif", 
               min_size = 1, max_size = 100, w = 20)

pdf("flatten_nucleus.pdf")
zimg_files = get_files("G08.+F001L01A01Z[0-9].+.tif", files)
max_contrast = contrastProjection(imageStack = combine_imgs(zimg_files), w_x = 150, w_y = 150,
                                  smoothing = 5, brushShape = "box", interpolation = 5)

# cimg = combine_imgs(zimg_files)
# mimg = apply(cimg, c(1,2), max)
# plot_imgN(list(rescale(max_contrast), rescale(mimg), rgbImage(blue = rescale(max_contrast), red = rescale(mimg))))

chan_names = c("nucleus", "probeA", "probeB")
# chans = as.list(get_files("G08_T0001F001L01A0[1-3]Z09C0[1-3].tif", files))
# names(chans) = chan_names
# # ccols = list(rgb(1,0,0), rgb(0,1,0), rgb(0,0,1))
ccols = list(c(1,0,0), c(0,1,0), c(0,0,1))
names(ccols) = chan_names

# setwd("~/R/PG_images/")
# processFISH(chans$nucleus,
#             channelSignals = chans[2:3], channelColours = ccols[2:3], writedir = "PG_1_", 
#             sizeNucleus = c(1500, 20000), sizeProbe = c(1, 100))

# setwd("~/R/PG_images/")
# processFISH(normalizePath("flatten_nucleus.tiff"),
#             channelSignals = as.list(normalizePath(c("flatten_probeA.tiff", "flatten_probeB.tiff"))), channelColours = ccols[2:3], writedir = "PG_1_",
#             sizeNucleus = c(1500, 20000), sizeProbe = c(1, 100))
