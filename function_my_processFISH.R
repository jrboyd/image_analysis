# function (combinedImg, writedir, bgCorrMethod = list(1, 100), 
#           channelSignals = NULL, channelColours = NULL, sizeNucleus = c(5, 
#                                                                         15000), sizeProbe = c(5, 100), gaussigma = 20, outputImageFormat = ".png") 
# {
library(BiocFileCache)
setwd("~/R/PG_images/")
source("functions.R")
source("function_distances.R")
source("helper_my_processFISH.R")
library(BiocFileCache)
library(digest)
library(FISHalyseR)

my_flattenFISH = function(writedir, base_key, bfc = BiocFileCache("~/.cache_FISH"), force_overwrite = FALSE){
    nuc_f = make_flat(file.path(writedir, "flatten_nucleus"), 
                      get_files(paste0(base_key, "A01Z[0-9].+.tif"), 
                                files = files), bfc, 150)
    pA_f = make_flat(file.path(writedir, "flatten_probeA"), 
                     get_files(paste0(base_key, "A02Z[0-9].+.tif"), 
                               files = files), bfc, 20)
    pB_f = make_flat(file.path(writedir, "flatten_probeB"), 
                     get_files(paste0(base_key, "A03Z[0-9].+.tif"), 
                               files = files), bfc, 20)
    return(list(nuc_f, pA_f, pB_f))
}

my_processFISH = function(writedir, combinedImg, channelSignals, bfc = BiocFileCache("~/.cache_FISH")){
    bgCorrMethod = list(1, 100)
    ccols = list(c(1,0,0), c(0,1,0), c(0,0,1))
    channelColours = ccols[2:3]
    names(channelColours) = c("probeA", "probeB")
    
    sizeNucleus = c(1500, 50000)
    sizeProbe = c(4, 100)
    gaussigma = 20
    outputImageFormat = ".png"
    
    
    bfcif = peakrefine::bfcif
    
    starttime <- Sys.time()
    FISHalyseR:::checkArguments(writedir, channelSignals, channelColours, 
                                sizeNucleus, sizeProbe)
    imageName <- basename(combinedImg)
    tryCatch({
        CombinedChannel = readImage(combinedImg)
    }, error = function(cond) {
        stop(paste("Composite image not found: ", combinedImg))
    })
    if (length(dim(CombinedChannel)) != 2) {
        imgG <- channel(CombinedChannel, "gray")
    }else {
        imgG <- CombinedChannel
    }
    dimension = dim(imgG)
    message("setting up directory structure")
    
    rootdir <- my_setUpDirectory(writedir, imageName, force_overwrite = TRUE)
    
    
    
    my_writeArguments(file.path(rootdir, "arguments.txt"),
        bgCorrMethod, channelColours, channelSignals, 
                                sizeNucleus, sizeProbe)
    
    key_load = digest(list(CombinedChannel, bgCorrMethod, gaussigma))
    c.load_CellMask = function(){
        load_CellMask(bgCorrMethod, imgG, gaussigma)
    }
    
    res = bfcif(bfc, paste0("raw_CellMask_", key_load), c.load_CellMask, force_overwrite = FALSE)
    CellMask = res[[1]]
    illu = res[[2]]
    
    # max_contrast = readImage(combinedImg)
    if(max(CellMask) > 1) CellMask = rescale(CellMask)
    bimg = binarize(CellMask)
    pimg = analyseParticles(bimg, max(sizeNucleus), min(sizeNucleus), 0)
    
    
    fill_img = 1-1*analyseParticles(1-1*pimg, Inf, 5000, 0)
    
    seg_cols = EBImage::colorLabels(EBImage::bwlabel(fill_img))
    
    png(file.path(rootdir, "nucleus_treatment.png"), width = 5, height = 4, units = "in", res = 150)
    plot_imgN(list(bimg, pimg, fill_img, seg_cols))
    dev.off()
    
    message("Extracting mask from the image....")
    
    # CellMask = FISHalyseR:::processCombined(CellMask, sizeNucleus)
    
    CellMask = fill_img#bfcif(bfc, paste0("process_CellMask_", key_nucleus, "_", digest(fill_img)), processCellMask)
    
    dim(CellMask) = c(dimension[1], dimension[2])
    message("Dividing nuclei")
    gradImage = bfcif(bfc, paste0("gradImage_", digest(CellMask)), function()applyPerim(CellMask))
    CellMask[gradImage > .03] <- 0
    CellMask = analyseParticles(CellMask, max(sizeNucleus), min(sizeNucleus), 0)
    seg_cols2 = EBImage::colorLabels(EBImage::bwlabel(CellMask))
    
    png(file.path(rootdir, "nucleus_segmentation.png"), width = 5, height = 4, units = "in", res = 150)
    plot_imgN(list(fill_img, seg_cols, CellMask, seg_cols2))
    dev.off()
    
    writeImage(CellMask, 
               file.path(rootdir, 
                         paste(imageName, "_Cellmask", outputImageFormat, 
                               sep = "")))
    ImCom <- rgbImage(red = 0 * gradImage, green = gradImage, 
                      blue = gradImage)
    message("Processing channel probes")
    channelImg <- list()
    img = readImage(channelSignals[[1]])
    if (FISHalyseR:::isColorImage(img) == 1 && length(channelSignals) == 1) {
        message("Extracting Red, Green, Blue from the combined cell image")
        channelImg[[1]] <- img[, , 1]
        channelImg[[2]] <- img[, , 2]
        channelImg[[3]] <- img[, , 3]
        channelColours = list(r = c(255, 0, 0), g = c(0, 255, 
                                                      0), b = c(0, 0, 255))
    } else {
        message("loading channels .....")
        for (i in 1:length(channelSignals)) {
            img = readImage(channelSignals[[i]])
            channelImg[[i]] <- img
        }
    }
    
    bwCM = bwlabel(CellMask)
    
    n <- names(channelColours)
    for (im in 1:length(channelImg)) {
        message(paste("loading channel : ", n[im], sep = ""))
        imgG <- FISHalyseR:::imageSubtract(channelImg[[im]], illu)
        
        ## method to trim outliers
        # tmp = imgG
        # i_max = quantile(tmp, .9999)
        # message(i_max)
        # imgG[imgG > i_max] = i_max
        # imgG = analyseParticles(binarize(imgG), sizeProbe[2], sizeProbe[1], isMask = 0)
        # browser()
        ## method to filter to nuclei
        tmp = imgG
        tmp[CellMask == 0] = 0
        # plot_imgN(list(imgG, tmp))
        imgG = tmp
        
        ## method to filter nuclei individually
        tmp = matrix(0, nrow = nrow(imgG), ncol = ncol(imgG))
        for(i in as.numeric(names(table(bwCM)[-1]))){
            message("in cell ", i)
            ri = range(which(apply(bwCM, 1, function(x) any(x == i))))
            ci = range(which(apply(bwCM, 2, function(x) any(x == i))))
            imgSub = imgG[seq(min(ri), max(ri)), seq(min(ci), max(ci))]
            
            imgSub_f = imgSub
            imgSub_f[imgSub_f == 0] = mean(imgSub_f[imgSub_f > 0])
            # plot_imgN(list(imgSub, imgSub_f))
            imgSub = imgSub_f
            imgSub_p = imgSub
            # t = calculateThreshold(imgSub_p)
            t = quantile(imgSub_p, .995)
            imgSub_p[imgSub_p<t] <- 0
            imgSub_p[imgSub_p>=t] <- 1
            imgSub_pp = analyseParticles(imgSub_p, sizeProbe[2], 5, isMask = 0)
            # imgSub_p = (binarize(imgSub))
            # plot_imgN(list(imgSub, imgSub_p, imgSub_pp))
            imgG[seq(min(ri), max(ri)), seq(min(ci), max(ci))] = imgSub_pp
        }
        plot_imgN(imgG)
        
        channelImg[[im]] <- imgG
        writeImage(imgG, file.path(rootdir, 
                                   paste(imageName, "_", names(channelColours)[im], 
                                         outputImageFormat, sep = "")))
    }
    # ImageOverlay <- FISHalyseR:::combineImage(channelImg, ImCom, channelColours, 
    #                              dimension[1], dimension[2])
    
    ImageOverlay = rgbImage(channelImg[[1]], channelImg[[2]], CellMask)
    FullOverlay = rgbImage(readImage(channelSignals[[1]]), readImage(channelSignals[[2]]), readImage(combinedImg))
    
    
    
    
    png(file.path(rootdir, "overlays.png"), width = 3, height = 6, units = "in", res = 150)
    plot_imgN(list(FullOverlay, ImageOverlay))
    dev.off()
    writeImage(ImageOverlay, file.path(rootdir, paste(imageName, "_Overlay.png", 
                                                      sep = "")))
    writeImage(ImageOverlay, file.path(rootdir, paste(imageName, ".png", sep = "")))
    message("extracting features from the cell mask image")
    LabelCellmask <- bwlabel(CellMask)
    tab = table(LabelCellmask)
    keep = as.numeric(names(tab)[tab > 50])
    
    LabelCellmask[!LabelCellmask %in% keep] = as.numeric(names(tab)[tab == max(tab)])
    LabelCellmask = bwlabel(LabelCellmask)
    message("saving label image with cell id")
    FISHalyseR:::plotLabelMatrix(LabelCellmask, 
                                 file.path(rootdir, 
                                           paste0(imageName, "_Label.jpg")))
    message("computing features .... ")
    FeaturesCellmask <- computeFeatures(LabelCellmask, CellMask, 
                                        methods.noref = c("computeFeatures.moment", "computeFeatures.shape"), 
                                        properties = FALSE)
    cellIndex = dim(FeaturesCellmask)[2] + 1
    cellInfo <- list()
    for (i in 1:length(channelImg)) {
        message(paste("Processing channel ", n[[i]], sep = ""))
        cellStr <- FISHalyseR:::GetStain(LabelCellmask, channelImg[[i]], cellIndex)
        cellInfo[[n[[i]]]] <- cellStr
        save(cellStr, file = file.path(rootdir, 
                                       "RData", 
                                       paste0(imageName, 
                                              "_", 
                                              n[i], 
                                              ".RData")))
    }
    maxCellProbe <- matrix(data = NA, 
                           nrow = max(LabelCellmask), 
                           ncol = length(n))
    for (iCell in 1:max(LabelCellmask)) {
        for (i in 1:length(channelColours)) {
            maxCellProbe[iCell, i] <- length(which(cellInfo[[n[i]]]$FCells[, 
                                                                           cellIndex] == iCell))
        }
    }
    message("Computing distances between probes")
    dMatChannels <- list()
    maxProbeDist <- list()
    datanamelist <- list()
    probeAreaNames <- list()
    
    for (i in 1:length(channelColours)) {
        for (j in 1:length(channelColours)) {
            if (i == j) {
                
                message(paste("computing distance between ", 
                              n[i], " and ", n[j], sep = ""))
                same1 <- cellInfo[[n[i]]]$FCells
                same2 <- cellInfo[[n[j]]]$FCells
                cname <- paste(n[i], n[j], sep = "_")
                distMatSame <- my_GetDistances(same1, same2, cellInfo[[n[i]]], 
                                               cellInfo[[n[j]]], isOneColour = 1)
                dMatChannels[[cname]] <- distMatSame
                # if(i == 2) browser()
                maxProbeDist[[cname]] <- my_findProbeMaxLength(distMatSame, 1)
                if (maxProbeDist[[cname]][1] != 0) {
                    ind <- 1
                    for (x in maxProbeDist[[cname]][1]:1) {
                        for (y in 1:x) {
                            datanamelist[length(datanamelist) + 1] <- paste(n[i], 
                                                                            ind, " ", n[j], y + ind, sep = "")
                        }
                        ind <- ind + 1
                    }
                    mprobesame <- maxProbeDist[[cname]][1] + 1
                    for (x in 1:mprobesame) {
                        probeAreaNames[length(probeAreaNames) + 1] <- paste("A", 
                                                                            n[i], x, sep = "")
                    }
                }
                else {
                    probeAreaNames[length(probeAreaNames) + 1] <- "not found"
                    datanamelist[length(datanamelist) + 1] <- "not found"
                }
            }
            else {
                if (j <= i) {
                    next
                }
                else {
                    message(paste("computing distance between ", 
                                  n[i], " and ", n[j], sep = ""))
                    diff1 <- cellInfo[[n[i]]]$FCells
                    diff2 <- cellInfo[[n[j]]]$FCells
                    cname <- paste(n[i], n[j], sep = "_")
                    distMatDiff <- my_GetDistances(diff1, diff1, cellInfo[[n[i]]], 
                                                   cellInfo[[n[j]]], isOneColour = 0)
                    dMatChannels[[cname]] <- distMatDiff
                    # maxProbeDist[[cname]] <- FISHalyseR:::findProbeMaxLength(distMatDiff, 
                                                                             # 0)
                    maxProbeDist[[cname]] <- my_findProbeMaxLength(distMatDiff, 
                                                                             0)
                    
                    if (maxProbeDist[[cname]][1] != 0 && maxProbeDist[[cname]][2] != 
                        0) {
                        for (x in 1:maxProbeDist[[cname]][1]) {
                            for (y in 1:maxProbeDist[[cname]][2]) {
                                datanamelist[length(datanamelist) + 1] <- paste(n[i], 
                                                                                x, " ", n[j], y, sep = "")
                            }
                        }
                    }
                    else {
                        datanamelist[length(datanamelist) + 1] <- "not found"
                    }
                }
            }
        }
    }
    
    # browser()
    
    Analysis.data <- data.frame()
    message("saving data to csv file ............ ")
    for (iCell in 1:max(LabelCellmask)) {
        message(paste("processing cell id: ", iCell))
        # Nucleus <- FISHalyseR:::GetNucleus(FullOverlay, iCell, FeaturesCellmask)
        Nucleus <- my_GetNucleus(FullOverlay, iCell, FeaturesCellmask)
        writeImage(Nucleus, 
                   file.path(rootdir, "cells", 
                             paste0("CellID", iCell, ".png")))
        # Nucleus <- FISHalyseR:::GetNucleus(ImageOverlay, iCell, FeaturesCellmask)
        Nucleus <- my_GetNucleus(ImageOverlay, iCell, FeaturesCellmask)
        writeImage(Nucleus, 
                   file.path(rootdir, "cells", 
                             paste0("CellID", iCell, "_a.png")))
        probeDist <- list()
        probeArea <- list()
        for (i in 1:length(channelColours)) {
            for (j in 1:length(channelColours)) {
                if (i == j) {
                    cname <- paste(n[i], n[j], sep = "_")
                    probeDist[[cname]] <- FISHalyseR:::CreateOutputDistanceVector(dMatChannels[[cname]], 
                                                                                  iCell, maxProbeDist[[cname]], isOneColour = 1)
                    mprobesame <- maxProbeDist[[cname]][1] + 1
                    probeArea[[cname]] <- my_CreateOutputAreaVector(cellInfo[[n[[i]]]]$FCells, 
                                                                              iCell, maxProbeDist[[cname]][1])
                } else {
                    if (j <= i) {
                        next
                    } else {
                        cname <- paste(n[i], n[j], sep = "_")
                        probeDist[[cname]] <- FISHalyseR:::CreateOutputDistanceVector(dMatChannels[[cname]], 
                                                                                      iCell, maxProbeDist[[cname]], isOneColour = 0)
                    }
                }
            }
        }
        Analysis.data <- rbind(Analysis.data, data.frame(paste(imageName, 
                                                               ".png", sep = ""), iCell, FeaturesCellmask[iCell, 
                                                                                                          "x.0.m.eccentricity"], rbind(maxCellProbe[iCell, 
                                                                                                                                                    ]), rbind(unlist(probeDist, use.names = FALSE)), 
                                                         FeaturesCellmask[iCell, "x.0.m.cx"], FeaturesCellmask[iCell, 
                                                                                                               "x.0.m.cy"], FeaturesCellmask[iCell, "x.0.s.area"], 
                                                         FeaturesCellmask[iCell, "x.0.s.perimeter"], FeaturesCellmask[iCell, 
                                                                                                                      "x.0.s.radius.mean"], rbind(unlist(probeArea, 
                                                                                                                                                         use.names = FALSE))))
    }
    dataColNames <- list()
    dataColNames[1] <- "filename"
    dataColNames[2] <- "nucleus ID"
    dataColNames[3] <- "eccentricity"
    for (i in 1:length(channelColours)) {
        dataColNames[length(dataColNames) + 1] <- paste("num of ", 
                                                        n[i], " probes", sep = "")
    }
    dataColNames <- c(dataColNames, datanamelist)
    dataColNames[length(dataColNames) + 1] <- "X center of mass"
    dataColNames[length(dataColNames) + 1] <- "Y center of mass"
    dataColNames[length(dataColNames) + 1] <- "area of nucleus"
    dataColNames[length(dataColNames) + 1] <- "perimeter of nucleus"
    dataColNames[length(dataColNames) + 1] <- "radius of nucleus"
    dataColNames <- c(dataColNames, probeAreaNames)
    colnames(Analysis.data) <- dataColNames
    write.table(Analysis.data, file.path(rootdir, "csv", 
                                         paste0(imageName, "_data.csv")), 
                sep = ",", row.names = FALSE, col.names = TRUE)
    save(Analysis.data, file = file.path(rootdir, "RData", 
                                         paste0(imageName, "_final.RData")))
    message("Done processing fish data ...")
    endtime <- Sys.time()
    timetaken <- endtime - starttime
    print(timetaken)
    
    invisible(rootdir)
}

my_writeArguments = function (fileConn, bgCorrMethod, channelColours, channelSignals, sizeNucleus, 
                              sizeProbe) 
{
    string1 = "Argument list"
    string1 = c(string1, "---- Background subtraction ----")
    if (is.list(bgCorrMethod)) {
        if (bgCorrMethod[[1]] == 1) {
            string1 <- c(string1, "Multiple Gaussian blurring")
        }
        else if (bgCorrMethod[[1]] == 2) {
            string1 <- c(string1, "Subtract user-specified illumination image")
        }
        else if (bgCorrMethod[[1]] == 3) {
            string1 <- c(string1, "Multidimensional Illumination Correction")
        }
    }
    if (is.list(channelColours) && is.list(channelSignals)) {
        string1 = c(string1, "Probes information")
        for (i in 1:length(channelColours)) {
            chval <- paste(channelColours[[i]], collapse = ", ")
            string1 = c(string1, paste(names(channelColours)[i], 
                                       chval, channelSignals[[i]], sep = " - "))
        }
    }
    string1 = c(string1, "---- Analyse Cells ----")
    string1 = c(string1, paste("Maximum area", sizeNucleus[1], 
                               sep = " - "), paste("Minimum area", sizeNucleus[2], sep = " - "))
    string1 = c(string1, "---- Analyse Probes ----")
    string1 = c(string1, paste("Maximum area", sizeProbe[1], 
                               sep = " - "), paste("Minimum area", sizeProbe[2], sep = " - "))
    
    writeLines(string1, fileConn)
}

my_CreateOutputAreaVector = function (Channel, iCell, MaxNElements) 
{
    if (MaxNElements != 0) {
        celIndex = ncol(Channel)
        nProbes <- length(Channel[which(Channel[, celIndex] == 
                                            iCell), "x.0.s.area"])
        aVector <- matrix(data = NA, ncol = MaxNElements + 1, 
                          nrow = 1)
        if (nProbes > 0) {
            aVector[1, 1:nProbes] <- Channel[which(Channel[, 
                                                           celIndex] == iCell), "x.0.s.area"]
        }
    }
    else {
        aVector = matrix(data = c(-999), nrow = 1, ncol = 1)
    }
    return(aVector)
}

my_findProbeMaxLength = function(Channel, isOneColour) 
{
    # browser()
    if(nrow(Channel) == 0) return(c(1, 1))
    cells <- unique(Channel[, 2])
    maxProbes <- c(0, 0)
    if (unique(!is.na(cells))) {
        for (iCell in cells) {
            if (isTRUE(dim(Channel[which(Channel[, 2, drop = FALSE] == iCell), 
                                   ])[1] > 1)) {
                SubChannel <- Channel[which(Channel[, 2, drop = FALSE] == iCell), 
                                      ]
            }
            else {
                SubChannel <- 
                    Channel[which(
                        Channel[, 2] == iCell), ,drop = FALSE]
            }
            if (isOneColour == 1) {
                lmax <- length(unique(SubChannel[, 3]))
                vlen <- dim(SubChannel)[1]
                if (maxProbes[1] < lmax) {
                    maxProbes[1] <- lmax
                }
                if (maxProbes[2] < vlen) {
                    maxProbes[2] <- vlen
                }
            }
            else {
                col1Probe <- length(unique(SubChannel[, 3, drop = FALSE]))
                r <- dim(SubChannel)[1]
                col2Probe <- r/col1Probe
                if (maxProbes[1] < col1Probe) {
                    maxProbes[1] <- col1Probe
                }
                if (maxProbes[2] < col2Probe) {
                    maxProbes[2] <- col2Probe
                }
            }
        }
    }
    return(maxProbes)
}
