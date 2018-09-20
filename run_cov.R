#test by test file
library(seqsetvis)
for(f in dir("tests/testthat/", pattern = "test", full.names = TRUE)){
    print(f)
    testthat::test_file(f, reporter = testthat::RstudioReporter)
}

# parse_version = function(){
#     ver = read.table("DESCRIPTION", sep = "\n", stringsAsFactors = FALSE)[,1]
#     ver = ver[grepl("Version: ", ver)]
#     ver = sub("Version: ", "", ver)
#     ver
# }
#
# roxygen2::roxygenise()
#
# devtools::build(vignettes = TRUE)
# tbal = paste0("../seqsetvis_", parse_version(), ".tar.gz")
# if(!file.exists(tbal)) stop("can't find ", tbal)
# system(paste("R CMD check", tbal))
# system(paste("R CMD BiocCheck", tbal))
# # browseURL("/slipstream/home/joeboyd/R/x86_64-pc-linux-gnu-library/3.4/seqsetvis/doc/index.html")
#
# #build properly using Rbuildignore
# # install.packages()
# devtools::install(build_vignettes = TRUE)
# devtools::build_vignettes()
# # devtools::run_examples()
# browseURL("/slipstream/home/joeboyd/R/x86_64-pc-linux-gnu-library/3.4/seqsetvis/doc/seqsetvis_overview.html")
#
# a = covr::file_coverage(dir("R", full.names = T, pattern = "R$"), dir("tests/testthat/", full.names = T, pattern = "R$"))
# covr::report(a)

# covr::codecov()
# covr::codecov(token = "ec9fa5ec-1e23-4e3e-82bb-260ed1ee514a")
