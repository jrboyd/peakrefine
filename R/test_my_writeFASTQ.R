source("dev_chipsim_functions.R")
message("writing fastq")
load("dev_my_writeFASTQ.save")
lim = 2e4
read = read[seq(lim)]
quality = quality[seq(lim)]
name = name[seq(lim)]

st_my = system.time({
    my_writeFASTQ.dt
    (file = sub(".fastq", "_my.fastq", file),
                  read = read, quality = quality, name = name)
})
st_def = system.time({
    ChIPsim::writeFASTQ(file = sub(".fastq", "_def.fastq", file),
               read = read, quality = quality, name = name)
})

def_md5 = tools::md5sum("tmp_def.fastq")
my_md5 = tools::md5sum("tmp_my.fastq")

stopifnot(def_md5 == my_md5)

message("my method")
print(st_my)
message("default method")
print(st_def)

message("ratio")
print(st_def[3] / st_my[3])
