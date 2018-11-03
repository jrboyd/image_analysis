load_CellMask = function(bgCorrMethod, imgG, gaussigma){
    illu <- NULL
    CellMask <- NULL
    if (is.list(bgCorrMethod)) {
        if (bgCorrMethod[[1]] > 0) {
            if (bgCorrMethod[[1]] == 1) {
                message("Multiple Gaussian blurring background subtraction...")
                if (is.numeric(bgCorrMethod[[2]])) {
                    illu = gblur(imgG, sigma = gaussigma)
                    for (i in 1:bgCorrMethod[[2]]) {
                        illu = gblur(illu, sigma = gaussigma)
                    }
                    CellMask <- imgG - illu + mean(illu)
                    CellMask[CellMask < 0] <- 0
                }else {
                    stop("Background subtraction arguments are not correct")
                }
            }
            else if (bgCorrMethod[[1]] == 2) {
                message("Using user specified illumination corrcetion image...")
                if (file.exists(bgCorrMethod[[2]])) {
                    illu = readImage(bgCorrMethod[[2]])
                    illu = gblur(illu, sigma = gaussigma)
                    CellMask = imageSubtract(imgG, illu)
                }else {
                    stop("Background subtraction arguments are not correct")
                }
            }else if (bgCorrMethod[[1]] == 3) {
                print("Computing illumination correction image from stack")
                illu <- computeIlluminationCorrection(bgCorrMethod[[2]],
                                                      bgCorrMethod[[3]], bgCorrMethod[[4]])
                CellMask = imageSubtract(imgG, illu)
            }else {
                stop("argument should 1 or 2 or 3")
            }
        }else {
            message("No illumination correction method selected")
            illu <- imageSubtract(imgG, imgG)
            CellMask <- imgG
        }
    } else {
        stop("bgCorrMethod should be type list")
    }
    list(CellMask, illu)
}

my_setUpDirectory = function (writedir, imageName, force_overwrite = FALSE) 
{
    maindir = paste(writedir)
    if(dir.exists(maindir) && !force_overwrite) stop("maindir exists")
    suppressWarnings({
        dir.create(file.path(maindir), recursive = TRUE)
        print(paste("Setting work directory to ", file.path(maindir), 
                    sep = ""))
        # setwd(file.path(maindir))
        dir.create(file.path(maindir, "csv"))
        dir.create(file.path(maindir, "cells"))
        dir.create(file.path(maindir, "RData"))
    })
    return(maindir)
}

processCellMask = function(fill_img, sizeNucleus){
    FISHalyseR:::processCombined(fill_img, sizeNucleus)
}

applyPerim = function(CellMask){
    FISHalyseR:::bwPerim(CellMask)
}