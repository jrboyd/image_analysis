#to be run after generating csvs
library(data.table)
library(magrittr)
library(ggplot2)

micronpp = 0.21

res_dir = "results_pipeline_HLB_v6"
# res_dir = "results_pipeline_HLB_v1"
todo = file.path(dir(res_dir, full.names = TRUE), "csv/flatten_nucleus.tiff_data.csv")
# stopifnot(all(file.exists(todo)))
todo = todo[file.exists(todo)]

res_nam = todo %>% sub(".+results_", "", .) %>% sub("/csv.+", "", .)

null_dt = data.table(id = 1L, p1 = "a", p2_group = "a", p1_group = "a", p2 = "a", value = 1, rankOrder = 1)
null_dt = null_dt[-1]

load_FISH_res = function(f){
    dt = fread(f)
    colnames(dt)[2] = "id"
    dt = dt[`num of probeA probes` == 2 & `num of probeB probes` == 2]
    
    if(nrow(dt) == 0) return(null_dt)
    
    k = grepl("probe[AB][1-9] probe[AB][1-9]", colnames(dt))
    if(sum(k) == 0) return(null_dt)
    cn = c("id", colnames(dt)[k])
    dist_dt = dt[, mget(cn)]
    dist_dt = suppressWarnings({data.table::melt(dist_dt, id.vars = "id")[]})
    dist_dt = dist_dt[!is.na(value)]
    # dist_dt[, variable := gsub("probe", "", variable)]
    # dist_dt[, variable := gsub("[0-9]+", "", variable)]
    dist_dt[, c("p1", "p2") := tstrsplit(variable, " ")]
    dist_dt = unique(dist_dt)
    # dist_dt[, p1 := gsub("[0-9]+", "", p1)]
    # dist_dt[, p2 := gsub("probeB[0-9]+", "probeB", p2)]
    dist_dt = dist_dt[grep("probeA", p1)]
    dist_dt = dist_dt[grep("probeB", p2)]
    dist_dt[, p1_group := sub("[0-9]", "", p1)]
    dist_dt[, p2_group := sub("[0-9]", "", p2)]
    
    dist_dt = dist_dt[, .(p1_group, p2, value, rankOrder = rank(value)), by = .(id, p1, p2_group)]
    dist_dt
}

total_cells = sapply(todo, function(f){
    # print(f)
    nrow(fread(f))
})
names(total_cells) = res_nam

all_res = lapply(todo, function(f){
    print(f)
    load_FISH_res(f)
})
names(all_res) = res_nam

sapply(all_res, is.null)


valid_cells = sapply(all_res, nrow)/4

all_dt = rbindlist(all_res, use.names = TRUE, idcol = "source")
#all_dt = all_dt[rankOrder == 1]

all_dt[, c("cell", "p1_group", "p2_group", "field") := tstrsplit(basename(source), "_")]
# all_dt[, c("source", "well") := tstrsplit(source, "_")]
all_dt[,p1 := sub("[a-zA-Z]", "", "p1")]
all_dt[,p2 := sub("[a-zA-Z]", "", "p2")]

all_dt[grepl("10A_2P4-488_129.+1$", source) ]


all_dt[, mean(value) , by = .(source)]
all_dt[, comb := paste(p1_group, p2_group)]
# bad = all_dt[value > 20, .(source, id)]

all_dt$cell = factor(all_dt$cell)
all_dt$field = factor(all_dt$field, levels = sort(as.numeric(unique(all_dt$field))))

all_dt$is_outlier = FALSE
all_dt[rankOrder == 1, is_outlier := abs((value - mean(value)) / sd(value)) > 2, by = .(cell, comb)]
# all_dt[rankOrder == 1, is_outlier := value > quantile(value, .75), by = .(cell, comb)]

all_dt[, x := paste(cell, field)]
all_dt$x = factor(all_dt$x, levels = unique(all_dt[order(field)][order(cell)]$x))

all_dt = all_dt[cell != "U2OS"]
all_dt$cell = droplevels(all_dt$cell)

all_dt[, value := value * micronpp]

pdf(paste0("plots_distances_", res_dir, ".pdf"), width = 10, height = 8)
ggplot(all_dt[rankOrder == 1 & !is_outlier][, .N, by = .(cell, comb)], 
       aes(x = cell, y = N, fill = cell)) + 
    geom_bar(stat = "identity") + 
    facet_wrap("comb", ncol = 2) + 
    labs(title = "Counts of spot pairs in valid nuclei", 
         subtitle = "expected 2 fish spots per probe\noutlier distances removed", x = "")

p_cdf = ggplot(all_dt[rankOrder == 1 & !is_outlier], 
               aes(x = value, fill = cell, color = cell)) + 
    stat_ecdf(geom = "step") + 
    facet_wrap("comb", ncol = 2) + 
    labs(title = "CDF of distances per cell type", 
         subtitle = "long distances", 
         x = "um",
         y = "Cumulative Density")

p_cdf

p_cdf + 
    coord_cartesian(xlim = c(0, 10)*micronpp) +
    labs(subtitle = "short distances")

ggplot(all_dt[rankOrder == 1 & !is_outlier], 
       aes(x = value, fill = cell, group = cell)) +
    geom_density(alpha = .5) +
    scale_color_discrete(drop=FALSE) +
    coord_cartesian(xlim = c(0, 10)*micronpp) +
    labs(title = "Distribution of distances per cell type", subtitle = "short distances", x = "um") +
    # theme(strip.text = element_text(size = 8)) +
    facet_wrap("comb", ncol = 2)

ggplot(all_dt[rankOrder == 1 & !is_outlier], aes(x = value, fill = cell, group = cell)) +
    geom_density(alpha = .5) +
    scale_color_discrete(drop=FALSE) +
    # coord_cartesian(xlim = c(1.5, 10)*micronpp) +
    labs(title = "Distribution of distances per cell type", subtitle = "long distances, watch the y-axis", x = "um") +
    # theme(strip.text = element_text(size = 8)) +
    facet_wrap("comb", ncol = 2, scales = "free_y")

ggplot(all_dt[rankOrder == 1 & !is_outlier], aes(x = cell, y = value, fill = cell)) +
    geom_boxplot() +
    scale_fill_discrete(drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    theme(strip.text = element_text(size = 8)) +
    facet_wrap("comb", scales = "free", ncol = 2) +
    labs(y = "um", title = "Distribution of distances per cell type", x= "") 

ggplot(all_dt[rankOrder == 1 & !is_outlier], aes(x = x, y = value, fill = cell)) + 
    geom_boxplot(outlier.shape = NA) + #geom_jitter(width = .01) +
    facet_wrap("comb", scales = "free_y", ncol = 2) +
    theme(strip.text = element_text(size = 8), axis.text.x = element_text(angle = 30, size = 4)) +
    scale_fill_discrete(drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    labs(y = "um", title = "Distribution of distances per field", x= "field")
dev.off()




all_dt[, .N, by = .(cell, comb)]

all_dt[, mean(value), by = .(source)]
all_dt[, median(value), by = .(source)]

summary_dt = rbind(data.table(value = total_cells, sample = names(total_cells), stat = "total"),
                   data.table(value = valid_cells, sample = names(valid_cells), stat = "valid"))
# summary_dt[, group := sub("_[0-9]+$", "", basename(sample))]
summary_dt[, group := sub("_[0-9]+$", "", basename(sample))]
summary_dt = summary_dt[, .(value = sum(value)), by = .(group, stat)]
summary_dt[, c("cell", "p1", 'p2') := tstrsplit(group, "_")]
summary_dt[, pair := paste(p1, p2)]
# summary_dt[stat == "total", stat := "total nuclei"]
# summary_dt[stat == "valid", stat := "2and2 spots"]

summary_dt$cell = factor(summary_dt$cell)


lab_dt = dcast(summary_dt[, .(stat, value, cell, pair)], "pair+cell~stat")
lab_dt[, fraction := paste(format(round(100*valid / total, 1), nsmall = 1), "%")]

# summary_dt = summary_dt[, total_nuclei := total - valid]
pdf(paste0("summary_", res_dir, ".pdf"), width = 12, height = 4)
p = ggplot() + 
    geom_bar(data = summary_dt[stat == "total"], aes(x = cell, fill = stat, y = value), stat = "identity") +
    geom_bar(data = summary_dt[stat == "valid"], aes(x = cell, fill = stat, y = value), stat = "identity") +
    geom_label(data = lab_dt, aes(x = cell, y = valid, label = fraction)) +
    labs(y = "image analysis efficiency", title = "Total nuclei identified and fraction with valid probes") +
    # geom_bar(data = summary_dt, aes(x = cell, fill = stat, y = value), stat = "identity") + 
    facet_wrap("pair", strip.position = "bottom", nrow = 1)
p
dev.off()
p

summary_dt = rbind(data.table(value = total_cells, sample = names(total_cells), stat = "total"),
                   data.table(value = valid_cells, sample = names(valid_cells), stat = "valid"))
# summary_dt[, group := sub("_[0-9]+$", "", basename(sample))]
summary_dt[, group := basename(sample)]
summary_dt = summary_dt[, .(value = sum(value)), by = .(group, stat)]
summary_dt[, c("cell", "p1", "p2", "field") := tstrsplit(basename(group), "_")]
summary_dt[, pair := paste(p1, p2)]

lab_dt = dcast(summary_dt[, .(group, stat, value)], "group~stat")
lab_dt[, fraction := paste0(format(round(100*valid / total, 1), nsmall = 1), "%")]
lab_dt[, c("cell", "p1", "p2", "field") := tstrsplit(basename(group), "_")]
lab_dt[, pair := paste(p1, p2)]

agg_dt = copy(summary_dt)
agg_dt[, group := sub("_[0-9]+$", "", group)]
agg_dt = agg_dt[, .(value = sum(value)), by = .(group, stat, pair)]
agg_dt[, c("cell", "p1", "p2") := tstrsplit(basename(group), "_")]
agg_dt[, pair := paste(p1, p2)]

agg_lab_dt = dcast(agg_dt[, .(group, stat, value)], "group~stat")
agg_lab_dt[, fraction := paste0(format(round(100*valid / total, 1), nsmall = 1), "%")]
agg_lab_dt[, c("cell", "p1", "p2") := tstrsplit(basename(group), "_")]
agg_lab_dt[, pair := paste(p1, p2)]

# summary_dt = summary_dt[, total_nuclei := total - valid]
pdf(paste0("summary_", res_dir, ".pdf"), width = 18, height = 8)
p = ggplot() +
    geom_bar(data = agg_dt[stat == "total"], aes(x = pair, fill = stat, y = value), stat = "identity") +
    geom_bar(data = agg_dt[stat == "valid"], aes(x = pair, fill = stat, y = value), stat = "identity") +
    geom_label(data = agg_lab_dt, aes(x = pair, y = valid, label = fraction)) +
    geom_text(data = agg_lab_dt, aes(x = pair, y = total, label = total), nudge_y = 1) + 
    labs(y = "image analysis efficiency", title = "Total nuclei identified and fraction with valid probes") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    facet_wrap("cell")
print(p)
for(cl in unique(summary_dt$cell)){
    p = ggplot() + 
        geom_bar(data = summary_dt[stat == "total" & cell == cl], aes(x = group, fill = stat, y = value), stat = "identity") +
        geom_bar(data = summary_dt[stat == "valid" & cell == cl], aes(x = group, fill = stat, y = value), stat = "identity") +
        geom_label(data = lab_dt[cell == cl], aes(x = group, y = valid, label = fraction)) +
        geom_text(data = lab_dt[cell == cl], aes(x = group, y = total, label = total), nudge_y = 1) +
        labs(y = "image analysis efficiency", title = paste(cl, "-Total nuclei identified and fraction with valid probes")) +
        # geom_bar(data = summary_dt, aes(x = cell, fill = stat, y = value), stat = "identity") + 
        facet_grid("pair", scales = "free_x") + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    print(p)
}
dev.off()
