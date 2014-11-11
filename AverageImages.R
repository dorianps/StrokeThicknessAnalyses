AverageImages <- function(foldermask) {
##### useful ANTsR function that matches the one in shell
    imageList <- Sys.glob(foldermask)

    if (length(imageList) == 0) {  # no files found
        stop(paste("No files found for", foldermask))
    }

    suppressWarnings(rm(temp))  # delete variables if already exist
    suppressWarnings(rm(temp))

    message(paste("Averaging", length(imageList), "files"))
    for (i in 1:length(imageList)) {
        temp <- as.array(antsImageRead(imageList[i], 3))
        if (i == 1) { # this is first image
            total <- temp
            standarddim <- dim(total)
            cat(".")			# add a dot to show still working
            next
        }
        else {					# from second image on
            if (!identical(standarddim, dim(temp))) {  # stop if image dimension wrong
                stop(paste("Dimension mismatch!", imageList[1], "with dimensions", paste(dim(temp), collapse="x"), "different from first image", paste(standarddim, collapse="x")))
            }
            
            total <- total+temp
            cat(".")			# add a dot to show still working
        }
    }

    message("")
    temp <- antsImageRead(imageList[1], 3)
    temp <- as.antsImage(total / length(imageList))
}
