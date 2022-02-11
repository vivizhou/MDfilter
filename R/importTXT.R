#' Import txt file from ANNOVAR annotation
#'
#' This function loads the txt file from ANNOVAR annotation as a data frame.The structure of the data will be printed at the console automatically. Please if numeric variables were imported correctly.
#'
#'
#' @param filePath file path of the txt file from ANNOVAR annotation
#'
#' @return A data frame of ANNOVAR annotation data.
#' @export
importTXT <- function(filePath = NULL) {

  if(!is.null(filePath) ) {
    mydata <- read.table(filePath, header = TRUE, sep = "\t", quote = "",
                         stringsAsFactors = FALSE, na = ".")

  } else {
    stop("Please specify file paths")
  }

  #The code below removed quoting sign "" in some cells
    for (i in 1: ncol(mydata)) {
      if (is.character(mydata[, i])) {
        mydata[, i] <- gsub('[\"]','',mydata[, i])
      }
    }

  #print data structure
  cat(str(mydata))

  return(mydata)
}


