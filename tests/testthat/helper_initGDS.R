library(RAIDS)
library(withr)
library(gdsfmt)


## Create local GDS Sample File with some content
## The user should close the GDS file
## The file will be removed automatically
local_GDS_Sample_file <- function(path, env = parent.frame()) {
    GDS_file_tmp  <- createfn.gds(filename=path)
    defer(unlink(x=path, force=TRUE), envir = env)

    add.gdsn(GDS_file_tmp, "Ref.count", rep(10L, 12))
    add.gdsn(GDS_file_tmp, "Alt.count", rep(12L, 12))
    add.gdsn(GDS_file_tmp, "Total.count", rep(22L, 12))
    add.gdsn(GDS_file_tmp, "lap", rep(0.5, 12))
    sync.gds(GDS_file_tmp)

    return(GDS_file_tmp)
}

## Create local 1KG GDS File with some content
## The user should close the GDS file
## The file will be removed automatically
local_GDS_1KG_file <- function(path, env = parent.frame()) {
    GDS_file_tmp  <- createfn.gds(filename=path)
    defer(unlink(x=path, force=TRUE), envir = env)

    ## Create sample information initial
    add.gdsn(GDS_file_tmp, "sample.id", c("HTT101", "HTT102", "HTT103"))

    samples <- data.frame(sex=c(1, 1, 2), pop.group=c("GBR", "GIH", "GBR"),
                         superPop=c("EUR", "SAS", "EUR"), batch=rep(0, 3),
                         stringsAsFactors = FALSE)

    add.gdsn(GDS_file_tmp, "sample.annot", samples)

    add.gdsn(GDS_file_tmp, "sample.ref", c(1,0, 1))

    sync.gds(GDS_file_tmp)

    return(GDS_file_tmp)
}
