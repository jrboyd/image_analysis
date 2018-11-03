#to be run after generating csvs
library(data.table)
library(magrittr)
library(ggplot2)
res_dir = "results_pipeline_NCI_v1"
todo = file.path(dir(res_dir, full.names = TRUE), "csv/flatten_nucleus.tiff_data.csv")
stopifnot(all(file.exists(todo)))

cn = todo %>% sub(".+results_", "", .) %>% sub("/csv.+", "", .)

load_FISH_res = function(f){
    dt = fread(f)
    colnames(dt)[2] = "id"
    dt = dt[`num of probeA probes` == 2 & `num of probeB probes` == 2]
    
    k = grepl("probe[AB][1-9] probe[AB][1-9]", colnames(dt))
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

total_cells = sapply(todo, function(f)nrow(fread(f)))
names(total_cells) = cn

all_res = lapply(todo, load_FISH_res)
names(all_res) = c

valid_cells = sapply(all_res, nrow)/4

all_dt = rbindlist(all_res, use.names = TRUE, idcol = "source")
all_dt = all_dt[rankOrder == 1]

all_dt[, c("well", "field") := tstrsplit(source, "_")]


all_dt[, mean(value) , by = .(source)]
bad = all_dt[value > 20, .(source, id)]

ggplot(all_dt[value < 20], aes(x = source, y = value)) + geom_boxplot()

ggplot(all_dt[value < 20], aes(group = well, x = value, fill = well)) + geom_histogram()

# dist_dt = dist_dt[rankOrder == 1]
all_dt[id == 4]
dist_dt[id == 6]
dist_dt[id == 7]
dist_dt[id == 85]
