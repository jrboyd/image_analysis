require("BiocFileCache")
bfc = BiocFileCache("~/.cache_FISH")
bfcif = peakrefine::bfcif
library(digest)

plot_imgN = function(imgs, nrow = NULL, ncol = NULL, pad_size = 50, dpi = 150, apply_rescale = TRUE, apply_binarize = FALSE){
    stopifnot(length(pad_size) == 1)
    stopifnot(length(dpi) == 1)
    stopifnot(length(apply_rescale) == 1)
    stopifnot(length(apply_binarize) == 1)
    stopifnot(is.numeric(pad_size))
    stopifnot(is.numeric(dpi) || is.na(dpi))
    stopifnot(is.logical(apply_rescale))
    stopifnot(is.logical(apply_binarize))
    
    if(is.character(imgs))
        imgs = as.list(imgs)
    if(!is.list(imgs)) imgs = list(imgs)
    stopifnot(is.list(imgs))
    if(is.null(nrow) && is.null(ncol)){
        nrow = ceiling(length(imgs)^.5)
        ncol = ceiling(length(imgs) / nrow)
    }
    stopifnot(nrow * ncol >= length(imgs))
    
    tmp_f = tempfile()
    if(!is.na(dpi)){
        png(tmp_f, width = dev.size()[1], height = dev.size()[2], units = "in", res = dpi)
    }
    par(mar = rep(0, 4))
    par(bg = rgb(.2,.5,.7,1))
    plot(0:1, 0:1, 
         axes = FALSE, 
         xaxs = "i", yaxs = "i", 
         xlab = "", ylab = "",
         type = "n")
    xpad = 1 / pad_size / ncol
    ypad = 1 / pad_size / nrow
    xw = 1 / ncol
    yw = 1 / nrow
    
    for(i in seq_along(imgs)){
        img = imgs[[i]]
        if(is.character(img)){
            img = readImage(img)
        }
        if(apply_binarize){
            img = binarize(img)
        }
        if(apply_rescale){
            img = rescale(img)
        }
        nc = (i - 1) %% ncol
        nr = floor((i - 1) / ncol)
        message(nc, " ", nr)
        rasterImage(img, 
                    xleft = nc * xw + xpad, 
                    ybottom = 1 - (nr * yw - ypad + yw), 
                    xright = nc * xw - xpad + xw, 
                    ytop = 1 - (nr * yw + ypad), 
                    interpolate = FALSE)
    }
    if(!is.na(dpi)){
        dev.off()   
        fimg = png::readPNG(tmp_f)
        par(mar = rep(0, 4))
        par(bg = rgb(.2,.5,.7,1))
        plot(0:1, 0:1, 
             axes = FALSE, 
             xaxs = "i", yaxs = "i", 
             xlab = "", ylab = "",
             type = "n")
        rasterImage(fimg, xleft = 0, ybottom = 0, xright = 1, ytop = 1, interpolate = FALSE)
    }
}

binarize = function(img){
    # t = calculateMaxEntropy(img)
    # if(t == 0)
    t = calculateThreshold(img)
    img[img<t] <- 0
    img[img>=t] <- 1
    img
}

view_params = function(key){
    mat2 = strsplit(mat[grepl(key, mat[,5]),5], "[A-Z\\.]")
    mat2 = matrix(unlist(mat2), ncol = 8, byrow = T)
    uniqs = lapply(seq_len(ncol(mat2)), function(i){
        unique(mat2[,i])
    })
    names(uniqs) = c("empty", "T", "F", "L", "A", "Z", "C", "ext")
    uniqs
}

get_files = function(key, files){
    k = grepl(key, files)
    files[k]
}

key_plot = function(key){
    # browser()
    k = mat[,5] == key    
    tf = files[k]
    img = readImage(tf[1])
    img2 = readImage(tf[2])
    raw = img
    raw2 = img2
    img = binarize(img)
    img2 = binarize(img2)
    anaImg <- analyseParticles(img, 20000, 2000,0)
    anaImg2 <- analyseParticles(img2, 20000, 2000,0)
    plot_imgN(list(raw, img, anaImg, raw2, img2, anaImg2))
}

single_chan_load = function(img_f){
    stopifnot(length(img_f) == 1)
    stopifnot(file.exists(img_f))
    img = imageData(readImage(img_f))
    if(length(dim(img)) == 3){
        sums = sapply(seq_len(dim(img)[3]), function(z)sum(img[,,z]))
        k = which(sums > 0)[1]
        img = img[,,k]
    }
    stopifnot(is.matrix(img))
    img
}

combine_imgs = function(img_files){
    if(is.list(img_files)) img_files = unlist(img_files)
    img = single_chan_load(img_files[1])
    zimg = array(0, dim = list(nrow(img), ncol(img), length(img_files)))
    zimg[, , 1] = img
    for(i in seq_along(img_files)[-1]){
        zimg[, , i] = single_chan_load(img_files[i])
    }
    zimg[is.na(zimg)] = 0
    zimg
}

rescale = function(img){
    if(max(img) == 0 && min(img) == 0) return(img)
    img = img - min(img)
    img = img / max(img)
    img
}

analyze_params = function(name, key, min_size, max_size, w){
    pdf(paste0(name, ".pdf"))
    zimg_files = get_files(key, files)
    
    img_combined = bfcic(
        bfc, 
        digest(list(zimg_files, "combining", combine_imgs)), 
        function(){
            combine_imgs(zimg_files)
        }) 
    brush = "box"
    smoothing = 5
    interpolation = 5
    max_contrast = bfcif(
        bfc, 
        digest(list(zimg_files, w, smoothing, brush, interpolation, contrastProjection)), 
        function(){
            contrastProjection(imageStack = img_combined, w_x = w, w_y = w,
                               smoothing = smoothing, brushShape = brush, interpolation = interpolation)
        })
    cimg = combine_imgs(zimg_files)
    mimg = apply(cimg, c(1,2), max)
    plot_imgN(list(rescale(max_contrast), rescale(mimg)))
    dim(max_contrast)
    tiff::writeTIFF(rescale(max_contrast), paste0(name, ".tiff"))
    plot_imgN(zimg_files)
    bimg = binarize(max_contrast)
    mc = (max_contrast * 1/max(max_contrast))
    mc[bimg == 1] = 0
    plot_imgN(rgbImage(green = bimg, blue = rescale(mc)))
    pimg = analyseParticles(bimg, max_size, min_size, 0)
    tiff::writeTIFF(pimg, paste0(name, "_binary.tiff"))
    plot_imgN(rgbImage(green = pimg, red = bimg - pimg))
    dev.off()
}

make_flat = function(name, zimg_files, bfc, w){
    message("flattening ", name, "...")
    # max_contrast = bfcif(bfc, digest(list(zimg_files, w, 5, "box", 5)), function(){
    #  contrastProjection(imageStack = combine_imgs(zimg_files), w_x = w, w_y = w,
    #                        smoothing = 5, brushShape = "box", interpolation = 5)   
    # })
    FUN_contrast = function(){
        img = combine_imgs(zimg_files)
        max_img = apply(img, MARGIN = c(1,2), max)
        # cont_img = bfcif(bfc, digest(list(zimg_files, w, 5, "box", 5)), function(){
        #      contrastProjection(imageStack = combine_imgs(zimg_files), w_x = w, w_y = w,
        #                            smoothing = 5, brushShape = "box", interpolation = 5)
        #     })
        # plot_imgN(list(max_img, cont_img), nrow = 1, ncol = 2)
        max_img
    }
    max_contrast = bfcif(bfc, 
                         digest(list(zimg_files, "PG_maximum", FUN_contrast)), FUN_contrast, force_overwrite = FALSE)
    tiff::writeTIFF(rescale(max_contrast), paste0(name, ".tiff"))
    return(paste0(name, ".tiff"))
}

my_GetNucleus = function (Image, iCell, Features) 
{
    # browser()
    # MaxRadius <- Features[iCell, "x.0.s.radius.max"]
    # Width <- (MaxRadius + 5) * 2
    # Height <- (MaxRadius + 5) * 2
    # OriginX <- Features[iCell, 1] - Height/2
    # OriginY <- Features[iCell, 2] - Width/2
    # if (OriginX < 0) {
    #     OriginX <- 0
    # }
    # if (OriginY < 0) {
    #     OriginY <- 0
    # }
    # if (OriginX + Width > dim(Image)[1]) {
    #     Width <- round(Width - ((OriginX + Width) - dim(Image)[2]))
    # }
    # if (OriginY + Height > dim(Image)[2]) {
    #     Height <- round(Height - ((OriginY + Height) - dim(Image)[1]))
    # }
    # Nucleus <- Image[round(OriginX):round(OriginX + Width), 
    #                  round(OriginY):round(OriginY + Height), ]
    
    # browser()
    MaxRadius <- Features[iCell, "x.0.s.radius.max"]
    Width <- (MaxRadius + 5) * 2
    Height <- (MaxRadius + 5) * 2
    OriginX <- Features[iCell, 1] - Height/2
    OriginY <- Features[iCell, 2] - Width/2
    if (OriginX < 0) {
        OriginX <- 0
    }
    if (OriginY < 0) {
        OriginY <- 0
    }
    if (OriginX + Width > dim(Image)[1]) {
        Width <- round(Width - ((OriginX + Width) - dim(Image)[1]))
    }
    if (OriginY + Height > dim(Image)[2]) {
        Height <- round(Height - ((OriginY + Height) - dim(Image)[2]))
    }
    xs = floor(OriginX:(OriginX + Width))
    ys = floor(OriginY:(OriginY + Height))
    Nucleus <- Image[xs, ys, ]
    
    return(Nucleus)
}
