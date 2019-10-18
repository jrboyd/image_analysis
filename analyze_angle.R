library(magrittr)
library(data.table)
library(ggplot2)

get_ang_dt = function(f){
    attrib = f %>% dirname %>% dirname %>% basename %>% strsplit(., "_") %>% unlist
    pA = attrib[2]
    pB = attrib[3]
    
    load(f)
    if(length(dMatChannels) == 0){
        return(data.table(pair = character(), angle = numeric(), probeA = character(), probeB = character()))
    }
    ang_dt = rbindlist(lapply(dMatChannels, function(x)(data.table(angle = x[,4]))), use.names = TRUE, idcol = "pair")
    ang_dt[, pair := gsub("probeA", pA, pair)]
    ang_dt[, pair := gsub("probeB", pB, pair)]
    ang_dt[, c("probeA", "probeB") := tstrsplit(pair, "_")] 
    ang_dt[]
}

all_f = "results_pipeline_HLB_v6/" %>% dir(., full.names = TRUE) %>% paste0(., "/RData/dMatChannels.RData")
stopifnot(all(file.exists(all_f)))

ang_dt = pbapply::pblapply(all_f, get_ang_dt)
ang_dt = rbindlist(ang_dt)

theme_set(theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1)))

ggplot(ang_dt[probeA == probeB, .N, by = .(pair)], aes(x = pair, y = N)) + geom_bar(stat = "identity")

ggplot(ang_dt[probeA == probeB], aes(x = angle)) + geom_density() + facet_wrap("pair", ncol = 3)
ggplot(ang_dt[probeA == probeB], aes(x = angle)) + geom_histogram(bins = 60) + facet_wrap("pair", ncol = 3, scales = "free_y") 

ggplot(ang_dt[probeA != probeB], aes(x = angle)) + geom_density() + facet_wrap("pair", ncol = 3) 
ggplot(ang_dt[probeA != probeB], aes(x = angle)) + geom_histogram(bins = 60) + facet_wrap("pair", ncol = 3, scales = "free_y") 


dist_AeqB = hist(ang_dt[probeA == probeB & angle >=0 & angle <= 90]$angle, plot = FALSE, breaks = seq(0, 90, by = 10))$counts
dist_AnqB = hist(ang_dt[probeA != probeB ]$angle, plot = FALSE, breaks = seq(-90, 90, by = 10))$counts

chit = function(x){
    chisq.test(cbind(x, rep(sum(x) / length(x), length(x))))
}

chit(dist_AeqB)
chit(dist_AnqB)

bad = dist_AnqB
bad[2] = bad[2] * 3
chit(bad)
