#to be run after generating csvs
library(data.table)
library(magrittr)
library(ggplot2)
res_dir = "results_pipeline_NCI_v1"
# res_dir = "results_pipeline_HLB_v1"
todo = file.path(dir(res_dir, full.names = TRUE), "csv/flatten_nucleus.tiff_data.csv")
stopifnot(all(file.exists(todo)))

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
all_dt = all_dt[rankOrder == 1]

# all_dt[, c("project", "cell", "p1_group", "p2_group", "field") := tstrsplit(source, "_")]
all_dt[, c("source", "well") := tstrsplit(source, "_")]
# all_dt[,p1 := sub("[a-zA-Z]", "", "p1")]
# all_dt[,p2 := sub("[a-zA-Z]", "", "p2")]

all_dt


all_dt[, mean(value) , by = .(source)]
all_dt[, comb := paste(p1_group, p2_group)]
# bad = all_dt[value > 20, .(source, id)]

ggplot(all_dt, aes(x = source, y = value, fill = cell)) + geom_boxplot() + facet_wrap("comb", scales = "free", ncol = 2, strip.position = "right")
ggplot(all_dt, aes(x = cell, y = value, fill = cell)) + 
    geom_boxplot() + theme(strip.text = element_text(size = 8)) +
    facet_wrap("comb", scales = "free", ncol = 2, strip.position = "right")

all_dt[, .N, by = .(cell, comb)]

ggplot(all_dt[value < 20], aes(group = source, x = value, color = source)) + geom_density()
all_dt[, mean(value), by = .(source)]
all_dt[, median(value), by = .(source)]

summary_dt = rbind(data.table(value = total_cells, sample = names(total_cells), stat = "total"),
                   data.table(value = valid_cells, sample = names(valid_cells), stat = "valid"))
summary_dt[, group := sub("_[0-9]+$", "", sample)]
summary_dt = summary_dt[, .(value = sum(value)), by = .(group, stat)]
summary_dt[, c("project", "cell", "p1", 'p2') := tstrsplit(group, "_")]
summary_dt[, pair := paste(p1, p2)]
# summary_dt[stat == "total", stat := "total nuclei"]
# summary_dt[stat == "valid", stat := "2and2 spots"]

ggplot(summary_dt, aes(x = cell, fill = stat, y = value)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + facet_wrap("pair")

# dist_dt = dist_dt[rankOrder == 1]
all_dt[id == 4]
dist_dt[id == 6]
dist_dt[id == 7]
dist_dt[id == 85]
