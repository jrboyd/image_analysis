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
res_dir = "results_pipeline_HLB_v7"
root = file.path("~/R/PG_images/HLB_FISH/alltiff/")

files = dir(root, pattern = ".+tif", full.names = TRUE)
tmp = strsplit(basename(files), "[_\\.]")
lens = sapply(tmp, length)
m1 = matrix(unlist(tmp[lens == 9]), ncol = 9, byrow = TRUE)
m2 = matrix(unlist(tmp[lens == 10]), ncol = 10, byrow = TRUE)

files = c(files[lens == 9], files[lens == 10])

m1 = cbind(m1[,1:4], "1", m1[, 5:9])
mat = rbind(m1, m2)

mat[,2] = sub("2p4-488", "2P4-488", mat[,2])
mat[,4] = sub("NPAT-648", "NPAT-647", mat[,4])
mat = mat[, 1:8]
all_df = as.data.frame(cbind(mat, basename(files)))
colnames(all_df) = c("cell", "regA", "regB", "protein", "session", "field", "chan", "z", "file")
all_df$file = as.character(file.path(root, all_df$file))
# mat2 = matrix(unlist(strsplit(mat[,5], "[A-Z\\.]")), ncol = 8, byrow = T)
all_grps = apply(all_df[, -ncol(all_df)], 2, unique)
all_grps

bfc = BiocFileCache("~/.cache_FISH")

# g_todo = unique(mat[,4])
# f_todo = all_grps$F
# 
# g_todo = g_todo[2]
# f_todo = f_todo[5:6]

for(cl in all_grps$cell){
    message(cl)
    for(a in all_grps$regA){
        message("  ", a)
        for(b in all_grps$regB){
            message("    ", b)
            sub_df = subset(all_df, cell == cl & regA == a & regB == b)
            grps = paste(sub_df$session, sub_df$field, sub_df$protein)
            # print(sapply(unique(grps), function(g)sum(grps == g))  )
            i = 1
            for(g in unique(grps)){
                message("      ", g)
                k = g == grps
                grp_df = sub_df[k,]
                nuc_f = subset(grp_df, chan == "c0")$file
                prb1_f = subset(grp_df, chan == "c1")$file
                prb2_f = subset(grp_df, chan == "c2")$file
                
                
                odir = file.path(res_dir, paste(cl, a, b, i, sep = "_"))
                
                if(length(dir(file.path(odir, "csv/"), ".+csv") > 0)) next
                dir.create(odir, showWarnings = FALSE, recursive = TRUE)
                
                flat_fish = list(
                    
                    make_flat(file.path(odir, "flatten_nucleus"), 
                              nuc_f, 
                              bfc, 150),
                    
                    make_flat(file.path(odir, paste0("flatten_probeA_", a)), 
                              prb1_f, 
                              bfc, 20),
                    
                    make_flat(file.path(odir, paste0("flatten_probeB_", b)), 
                              prb2_f, 
                              bfc, 20)
                )
                
                my_processFISH(writedir = odir, 
                               combinedImg = flat_fish[[1]], 
                               channelSignals = flat_fish[2:3], 
                               bfc = bfc,
                               top_points = 30, sizeProbe = c(4, 100))
                i = i + 1
            }
        }
    }
}

