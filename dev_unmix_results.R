unmix_tiff = "~/data/TMA_Primary_0926_Runx2/nmfUNMIX_dapi_mancr_runx2_rbc.tif"

dat= tiff::readTIFF(unmix_tiff, as.is = TRUE)
dt = as.data.table(tiff::readTIFF(unmix_tiff))
setnames(dt, c("j", "i", "z", "val"))
dt
