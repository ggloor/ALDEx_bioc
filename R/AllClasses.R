### =========================================================================
### All classes
### =========================================================================


### -------------------------------------------------------------------------
### aldex.clr
###

setClass("aldex.clr",
    slots=c(
        reads="data.frame",
        conds="character",
        mc.samples="numeric",
        denom="vector",
        verbose="logical",
        useMC="logical",
        analysisData="list"
        )
)

validReads <- function(object) {
    if (length(duplicated(rownames(object)))==0) {
        TRUE
    }
    else {
        paste("Unable to create aldex.clr object.
            Duplicated row names:",rownames(object)[duplicated(row.names(object))])
    }
    if (length(duplicated(colnames(object)))==0) {
        TRUE
    }
    else {
        paste("Unable to create aldex.clr object.
            Duplicated column names:",colnames(object)[duplicated(row.names(object))])
    }
}

setValidity("aldex.clr",validReads)
