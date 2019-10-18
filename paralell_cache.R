library(BiocFileCache)
library(parallel)

bfc = BiocFileCache(".cache_paralell", ask = FALSE)

options("mc.cores" = 1)
res = mclapply(1:40, function(i){
    f = bfcnew(bfc, paste("entry", i))
    saveRDS(runif(10), f)
    i
})
table(sapply(res, class))


options("mc.cores" = 4)

res = mclapply(1:40, function(i){
    f = bfcnew(bfc, paste("entry", i))
    saveRDS(runif(10), f)
    i
})

files = sapply(1:40, function(i){
    bfcnew(bfc, paste("entry", i))
})

options("mc.cores" = 40)
mclapply(1:40, function(i){
    f = files[i]
    saveRDS(runif(10), f)
    i
})
