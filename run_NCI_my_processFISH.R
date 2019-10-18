library(FISHalyseR)
library(MaxContrastProjection)
library(BiocFileCache)
source("functions.R")
source("function_my_processFISH.R")

#FISH image processing
# 1) flatten z-stack images via max contrast
# 2) binarize, fill, size-select, and segment nuclei
# 3) outlier trim and binzarize FISH probe spots

#Measurement
# 4) within each nucleus, measure mean (default was minimum - many zeroes) distances between all spots
# 5) filter for nuclei containing exactly 2 spots for each probe

#Analysis
# 6) compare distrubutions of A-B distances per well/field

# f = system.file( "extdata", "SampleFISHgray.jpg", package="FISHalyseR")
root = file.path("~/knime-workspace/NCI_imaging/", 
                 "GLP-181017-Prachi-FISH-60X_20181017_121249/AssayPlate_PerkinElmer_CellCarrier-96 Ultra/")

files = dir(root, pattern = "Assay.+tif", full.names = TRUE)
mat = matrix(unlist(strsplit(basename(files), "_")), ncol = 5, byrow = T)
mat2 = matrix(unlist(strsplit(mat[,5], "[A-Z\\.]")), ncol = 8, byrow = T)
all_grps = apply(mat2[, 2:7], 2, unique)
names(all_grps) = c("T", "F", "L", "A", "Z", "C")

bfc = BiocFileCache("~/.cache_FISH")

g_todo = unique(mat[,4])
f_todo = all_grps$F

# g_todo = g_todo[2]
# f_todo = f_todo[5:6]

for(GROUP in g_todo){
    for(FIELD in f_todo){
        message(GROUP, " ", FIELD)
        base_key = paste0(GROUP, "_T0001F", FIELD, "L01")
        stopifnot(length(get_files(base_key, files)) == 42) #hardcode expected results size
        odir = paste0("results_pipeline_NCI_maxProjection/results_", GROUP, "_L", FIELD)
        dir.create(odir, showWarnings = FALSE, recursive = TRUE)
        flat_fish = my_flattenFISH(writedir = odir, base_key = base_key, bfc = bfc)
        # undebug(my_processFISH)
        # debug()
        my_processFISH(writedir = odir, combinedImg = flat_fish[[1]], channelSignals = flat_fish[2:3], bfc = bfc)
    }
}
