library(data.table)
library(ggplot2)
res_files = c(
    "results2/motif_res_MCF10A-blocked_Runx1-4336BF_pooled_data.save",
    "results2/motif_res_MCF10A-released_Runx1_pooled_data.save",
    "results2/motif_res_MCF10A-dmso_Runx1_pooled_data.save",
    "results2/motif_res_AF-MCF10A_RUNX1_pooled_data.save",
    "results2/motif_res_AF-MCF10AT1_RUNX1_pooled_data.save"
)
names(res_files) = c("10A-blocked", "10A-released", "10A-dmso", "AF-MCF10A", "AF-MCF10AT1")
f = res_files[1]
k = 1:4
k = T
all_qdt = lapply(res_files[k], function(f){
    load(f)
    nrow(qdt)
    qdt
})

lapply(all_qdt, nrow)

# plot(all_qdt[[1]]$flex_frag_len, all_qdt[[2]]$flex_frag_len)

all_flex_frag_len = lapply(all_qdt, function(x)x[, .(name, flex_frag_len)])
dt = rbindlist(all_flex_frag_len, use.names = TRUE, idcol = "sample")

ddt = dcast(dt, "name~sample")
ggplot(ddt, aes(x = `10A-blocked`, y = `10A-dmso`)) + geom_point()
GGally::ggpairs(ddt[,-1], aes(alpha = .2, shape = "21"))


wh = dt[flex_frag_len < 100 | flex_frag_len >= 250, unique(name)]
dt = dt[!name %in% wh]

ddt = dcast(dt, "name~sample")
ggplot(ddt, aes(x = `10A-blocked`, y = `10A-dmso`)) + geom_point()
GGally::ggpairs(ddt[,-1], aes(alpha = .2, shape = "21"))
