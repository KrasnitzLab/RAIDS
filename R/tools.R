#' @title Merge the pruning files by chromosome in one file
#'
#' @description TODO
#'
#' @param PATHPRUNED
#'
#' @param filePref
#'
#' @param fileOUT
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
#'
#' @keywords internal
groupChrPruning <- function(PATHPRUNED, filePref, fileOUT) {

    prunedList <- list()

    # Read the content of each file (one file per chromosome)
    for(i in seq_len(22)) {
        fileChr <- file.path(PATHPRUNED, paste0(filePref, i, ".rds"))

        if(file.exists(fileChr)) {
            prunedList[[i]] <- readRDS(fileChr)
        } else {
            stop("Problem with the file: ", fileChr)
        }
    }

    ## Merge the content of all files
    pruned <- do.call(c, prunedList)

    ## Save all the information into one file
    saveRDS(pruned, fileChr <- file.path(PATHPRUNED, fileOUT))
}
