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

groupChrPruning <- function(PATHPRUNED,
                            filePref,
                            fileOUT){

    prunedList <- list()

    for(i in seq_len(22)){
        fileChr <- file.path(PATHPRUNED, paste0(filePref, i, ".rds"))
        if(file.exists(fileChr)){
            prunedList[[i]] <- readRDS(fileChr)
        }else{
            stop(paste0("RPoblem with the file ", fileChr))
        }


    }
    pruned <- do.call(c, prunedList)
    saveRDS(pruned, fileChr <- file.path(PATHPRUNED, fileOUT))

}
