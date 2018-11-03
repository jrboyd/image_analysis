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

pdf("view_A01_all_Z_F1.pdf")
plot_imgN(get_files("G08.+F001.+A01", files), nrow = 4)
plot_imgN(get_files("G09.+F001.+A01", files), nrow = 4)
dev.off()
pdf("view_A02_all_Z_F1.pdf")
plot_imgN(get_files("G08.+F001.+A02", files), nrow = 4)
plot_imgN(get_files("G09.+F001.+A02", files), nrow = 4)
dev.off()
pdf("view_A03_all_Z_F1.pdf")
plot_imgN(get_files("G08.+F001.+A03", files), nrow = 4)
plot_imgN(get_files("G09.+F001.+A03", files), nrow = 4)
dev.off()

pdf("view_A01_all_Z_F2.pdf")
plot_imgN(get_files("G08.+F002.+A01", files), nrow = 4)
plot_imgN(get_files("G09.+F002.+A01", files), nrow = 4)
dev.off()
pdf("view_A02_all_Z_F2.pdf")
plot_imgN(get_files("G08.+F002.+A02", files), nrow = 4)
plot_imgN(get_files("G09.+F002.+A02", files), nrow = 4)
dev.off()
pdf("view_A03_all_Z_F2.pdf")
plot_imgN(get_files("G08.+F002.+A03", files), nrow = 4)
plot_imgN(get_files("G09.+F002.+A03", files), nrow = 4)
dev.off()

pdf("view_F.pdf")
plot_imgN(get_files("G08.+[1-3]L01A01Z09C01.tif", files), nrow = 2)
plot_imgN(get_files("G09.+[1-3]L01A01Z09C01.tif", files), nrow = 2)
dev.off()

pdf("view_A_and_C.pdf")
plot_imgN(get_files("G08_T0001F005L01A0[1-3]Z10", files), nrow = 2)
dev.off()

