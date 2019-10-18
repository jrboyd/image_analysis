my_CalcDistAll_centroid = function (chan1_raw, chan2_raw, chan1, chan2, s1, s2, CellId, isOneColour = 0) 
{
    
    Distance <- list()
    LengthS1 <- length(s1[["PList"]])
    LengthS2 <- length(s2[["PList"]])
    for (l1 in seq_len(LengthS1)) {
        if (s1[["PList"]][[l1]][3] == CellId) {
            for (l2 in seq_len(LengthS2)) {
                if (s2[["PList"]][[l2]][3] == CellId) {
                    if (isOneColour == 1) {
                        if (l2 <= l1) {
                            next
                        }
                    }
                    if (isTRUE(all.equal(unlist(s1[["PList"]][[l1]][1]), 
                                         unlist(s2[["PList"]][[l2]][1])))) {
                        print("Skipping probe ... comparison versus same probe")
                        next
                    }
                    else {
                        
                        LengthP1 <- length(unlist(s1[["PList"]][[l1]][[1]]))/2
                        LengthP2 <- length(unlist(s2[["PList"]][[l2]][[1]]))/2
                        mat_p1 <- unlist(s1[["PList"]][[l1]][[1]][[1]])
                        mat_p2 <- unlist(s2[["PList"]][[l2]][[1]][[1]])
                        
                        
                        mat_p1 = cbind(mat_p1, value = chan1_raw[mat_p1])
                        mat_p2 = cbind(mat_p2, value = chan1_raw[mat_p2])
                        
                        x1 = weighted.mean(mat_p1[,1], mat_p1[,3])
                        y1 = weighted.mean(mat_p1[,2], mat_p1[,3])
                        x2 = weighted.mean(mat_p2[,1], mat_p2[,3])
                        y2 = weighted.mean(mat_p2[,2], mat_p2[,3])
                        
                        d <- my_CalcDist(x1, y1, x2, y2)
                        
                        ang = asin(x = (y2 - y1) / d) / pi * 180
                        
                        paste(x2 - x1, y2 - y1, ang)
                        
                        Distance[length(Distance) + 1] <- list(list(Dist = d,
                                                                    ProbeId = s1[["PList"]][[l1]][2],
                                                                    Angle = ang))
                    }
                }
            }
        }
    }
    
    return(Distance)
}

#added dist metric function
my_CalcDistAll = function (s1, s2, CellId, isOneColour = 0, dist_metric_FUN = mean) 
{
    stop("don't run this")
    Distance <- list()
    Distance_centroid = list()
    LengthS1 <- length(s1[["PList"]])
    LengthS2 <- length(s2[["PList"]])
    for (l1 in seq_len(LengthS1)) {
        if (s1[["PList"]][[l1]][3] == CellId) {
            for (l2 in seq_len(LengthS2)) {
                if (s2[["PList"]][[l2]][3] == CellId) {
                    if (isOneColour == 1) {
                        if (l2 <= l1) {
                            next
                        }
                    }
                    if (isTRUE(all.equal(unlist(s1[["PList"]][[l1]][1]), 
                                         unlist(s2[["PList"]][[l2]][1])))) {
                        print("Skipping probe ... comparison versus same probe")
                        next
                    }
                    else {
                        LengthP1 <- length(unlist(s1[["PList"]][[l1]][[1]]))/2
                        LengthP2 <- length(unlist(s2[["PList"]][[l2]][[1]]))/2
                        ListP1 <- unlist(s1[["PList"]][[l1]][[1]][[1]])
                        ListP2 <- unlist(s2[["PList"]][[l2]][[1]][[1]])
                        tDistance <- list()
                        
                        dist_cent = my_CalcDist(
                            median(ListP1[,1]),
                            median(ListP1[,2]),
                            median(ListP2[,1]),
                            median(ListP2[,2])
                        )
                        
                        Distance_centroid[length(Distance_centroid) + 1] = 
                            list(list(Dist = dist_cent, ProbeId = s1[["PList"]][[l1]][2]))
                        
                        for (p1 in 1:LengthP1) {
                            for (p2 in 1:LengthP2) {
                                x1 <- ListP1[p1, 1]
                                x2 <- ListP2[p2, 1]
                                y1 <- ListP1[p1, 2]
                                y2 <- ListP2[p2, 2]
                                d <- my_CalcDist(x1, y1, x2, y2)
                                tDistance[length(tDistance) + 1] <- d
                            }
                        }
                        Distance[length(Distance) + 1] <- list(list(Dist = dist_metric_FUN(unlist(tDistance)), 
                                                                    ProbeId = s1[["PList"]][[l1]][2]))
                    }
                }
            }
        }
    }
    
    # # print(Distance)
    # print(sapply(Distance, function(x)x$Dist))
    # # print(Distance_centroid)
    # print(sapply(Distance_centroid, function(x)x$Dist))
    return(Distance)
}

my_GetDistances = function (Channel1, Channel2, Channel1PList, Channel2PList, isOneColour, Img1_raw, Img2_raw, Img1, Img2) 
{
    celIndex = ncol(Channel1)
    Channel1 <- na.exclude(Channel1)
    Channel2 <- na.exclude(Channel2)
    RowCounter <- 1
    if (length(Channel1) > 0) {
        IdOfCells <- unique(na.exclude(Channel1[, celIndex], 
                                       na.rm = TRUE))
        if (isOneColour == 1) {
            dMat <- matrix(data = NA, nrow = 30 * length(IdOfCells), 
                           ncol = 4)
        }
        else {
            dMat <- matrix(data = NA, nrow = 160 * length(IdOfCells), 
                           ncol = 4)
        }
        if (max(IdOfCells) >= 0) {
            for (iCell in IdOfCells) {
                IdChannel1 <- which(na.exclude(Channel1[, celIndex] == 
                                                   iCell))
                IdChannel2 <- which(na.exclude(Channel2[, celIndex] == 
                                                   iCell))
                MaxDots1 <- length(IdChannel1)
                MaxDots2 <- length(IdChannel2)
                # d <- my_CalcDistAll(Channel1PList, Channel2PList, 
                #                  iCell, isOneColour)
                d <- my_CalcDistAll_centroid(
                    Img1_raw, Img2_raw, Img1, Img2,
                    Channel1PList, Channel2PList, 
                    iCell, isOneColour)
                if (length(d) >= 1) {
                    for (l in 1:length(d)) {
                        dMat[RowCounter, 1] <- as.numeric(d[[l]][["Dist"]])
                        dMat[RowCounter, 2] <- iCell
                        dMat[RowCounter, 3] <- as.numeric(d[[l]][["ProbeId"]])
                        dMat[RowCounter, 4] <- as.numeric(d[[l]][["Angle"]])
                        RowCounter <- RowCounter + 1
                    }
                }
            }
            dMat <- unique(na.exclude(dMat))
        }
    }
    else {
        message("No probe found!")
        if (isOneColour == 1) {
            dMat <- matrix(NA, 3, 4)
        }
        else {
            dMat <- matrix(NA, 16, 4)
        }
    }
    return(dMat)
}

my_CalcDist = function (x1, y1, x2, y2) 
{
    xd = x2 - x1
    yd = y2 - y1
    Distance = sqrt(xd * xd + yd * yd)
    return(Distance)
}