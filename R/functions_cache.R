#' Conditionally run and cache FUN if results aren't already cached.  Returns result of FUN.
#'
#' @param bfc BiocFileCache
#' @param rname
#' @param FUN
#' @param force_overwrite
#'
#' @return
#' @export
#' @import BiocFileCache
#' @examples
bfcif = function(bfc, rname, FUN, force_overwrite = FALSE){
    # is rname in cache?
    if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname")) == 0){
        cache_path = BiocFileCache::bfcnew(bfc, rname = rname)

    }else{
        cache_path = BiocFileCache::bfcrpath(bfc, rname)
    }
    # does cached file exist?
    if(file.exists(cache_path) && !force_overwrite){
        load(BiocFileCache::bfcrpath(bfc, rname))
    }else{
        res = FUN()
        save(res, file = cache_path)
    }
    # return either new results or cached results
    res
}
