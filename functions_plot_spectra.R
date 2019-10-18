library(scales)
plot_rgb = 
    function(channel_dt,
             rv = "Runx2",
             gv = "mancr",
             bv = "nucleus", 
             norm1 = FALSE){
        stopifnot(c(rv, gv, bv) %in% colnames(channel_dt))
        stopifnot(c("i", "j") %in% colnames(channel_dt))
        #calculate rgb colors
        if(norm1){
            channel_dt[, rgbColor := rgb(rescale(get(rv), c(0, 1)), 
                                         rescale(get(gv), c(0, 1)), 
                                         rescale(get(bv), c(0, 1)))]    
        }else{
            channel_dt[, rgbColor := rgb(get(rv), 
                                         get(gv), 
                                         get(bv))]
        }
        
        
        #create legend placeholder
        leg_col = c("red", "green", "blue")
        names(leg_col) = c(rv, gv, bv)
        legend_dt = data.table::data.table(label = c(rv, gv, bv))
        
        p_rgb = ggplot() + 
            geom_point(data = legend_dt, 
                       aes(x = mean(channel_dt$i), 
                           y = mean(channel_dt$j), 
                           color = label)) +
            scale_color_manual(values = leg_col) +
            geom_raster(data = channel_dt, 
                        aes(x = i, 
                            y = j, 
                            fill = rgbColor)) + 
            scale_fill_identity() +
            scale_y_reverse() +
            theme(panel.background = element_blank()) +
            labs(title = "unmixed RGB", x= "pixel", y = "pixel") +
            coord_fixed()
        p_rgb
        
    }

# plot_rgb(rgb_dt)


plot_roi = 
    function(full_dt, 
             channel_dt, 
             xmin, xmax, ymin, ymax, 
             view_size = 50,
             rv = "Runx2",
             gv = "mancr",
             bv = "nucleus"){
        stopifnot(c("i", "j", "val", "channel") %in% colnames(full_dt))
        stopifnot(c(rv, gv, bv) %in% colnames(channel_dt))
        stopifnot(c("i", "j") %in% colnames(channel_dt))
        if(xmin > xmax){
            tmp = xmin
            xmin = xmax
            xmax = tmp
        }
        if(ymin > ymax){
            tmp = ymin
            ymin = ymax
            ymax = tmp
        }
        roi = c(xmin, xmax, ymin, ymax)
        if(view_size < diff(roi[1:2]) | view_size < diff(roi[3:4])){
            view_size = max(diff(roi[1:2]), diff(roi[3:4])) * 1.5
            
        }
        roi_exp = round(c(
            scales::expand_range(roi[1:2], add = (view_size - diff(roi[1:2]))/2),
            scales::expand_range(roi[3:4], add = (view_size - diff(roi[3:4]))/2)
        ))
        p_rgb.roi = 
            plot_rgb(channel_dt[i > roi_exp[1] & i <= roi_exp[2] & j > roi_exp[3] & j <= roi_exp[4]],
                     rv = rv, gv = gv, bv = bv) + 
            annotate("rect", 
                     xmin = roi[1]+.5, xmax = roi[2]+.5, 
                     ymin = roi[3]+.5, ymax = roi[4]+.5,
                     fill = NA, color = "white") #+
        # coord_fixed(xlim = roi_exp[1:2], ylim = roi_exp[3:4])
        
        full_dt.roi = full_dt[i > roi_exp[1] & i <= roi_exp[2] & j > roi_exp[3] & j <= roi_exp[4]]
        
        
        brks = scales::rescale(0:2/2, range(full_dt.roi$val))
        
        p_chan.roi = ggplot(full_dt.roi, 
                            aes(x = i, y = j, fill = val)) +
            geom_raster() +
            annotate("rect", 
                     xmin = roi[1]+.5, xmax = roi[2]+.5, 
                     ymin = roi[3]+.5, ymax = roi[4]+.5,
                     fill = NA, color = "blue") +
            scale_x_continuous(breaks = roi_exp[1:2]) +
            scale_y_reverse(breaks = roi_exp[3:4]) +
            scale_fill_gradientn(colours = c("white", "black"), breaks = brks) +
            facet_wrap(~channel, ncol = 8) +
            theme(panel.background = element_blank(), 
                  strip.background = element_blank(), 
                  strip.text = element_text(size = 6),
                  axis.text = element_text(size = 6), 
                  legend.position = "bottom", panel.spacing = unit(0, "npc")) +
            coord_fixed() +
            labs(x = "", y = "", fill = "raw value", title = "all channels")
        
        #line plots
        full_dt.roi[, wl := as.numeric(sub("nm", "", channel))]
        full_dt.roi[, in_roi := ifelse(i > roi[1] & i <= roi[2] & j > roi[3] & j <= roi[4], "In ROI", "outside")]
        p_chan_lines = ggplot(full_dt.roi, aes(x = wl, y = val, group = paste(i, j))) + geom_path() +
            facet_wrap(~in_roi, ncol = 1)
        
        cowplot::plot_grid(p_rgb.roi, p_chan.roi, p_chan_lines, rel_widths = c(2,1.5,1), nrow = 1)
        
    }
