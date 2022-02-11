#' Clean data
#'
#' This function prepare the data imported from ?import.txt.file for variant filtering.
#'
#' @param dataName Data imported by ?import.txt.file()
#' @param includeIDs IDs of the subjects included for analyses. This must match IDs used in the imported data file.If left empty, IDs will be extracted from the first to the last
#' @param excludeIDs IDs to be excluded from the data imported by ?import.txt.file()
#' @param removeMultiallelicSites True=remove multi-allelic sites; FALSE=keep multiallelic sites
#' @param remove.nocall.pct Variants with percentage of no calls larger than this value will be removed.
#' @param QD.threshold Variants with QD smaller than this value will be removed.
#' @param saveRData TRUE=Save the returned data frame into an RData file.
#' @param savePath Specify save path. If not specified the data will be saved at the working directory.
#' @param label Specify the label of the saved file.

#' @return A data frame of cleaned annotated data.
#'
#' @export
cleanData <- function(dataName,
                      includeIDs = NULL,
                      excludeIDs = NULL,
                      removeMultiallelicSites = FALSE,
                      remove.nocall.pct = 0.2, #percentage of samples with no call
                      QD.threshold = 5,
                      saveRData = FALSE,
                      savePath = NULL,
                      label = NULL) {

  #get subset of the data

  if (!is.null(excludeIDs)) {
    data_sub <- dataName[, ! colnames(dataName) %in% excludeIDs]
  } else {
    data_sub <- dataName
  }

  if (is.null(includeIDs)) {
    stop("Please specify IDs to be included (this must match the IDs in the imported data file.")
  } else {
    checkID <- which(!includeIDs %in% colnames(data_sub))
    if (length(checkID) > 0) {
      stop(paste0(paste0(includeIDs[checkID],collapse = ","), " does not match ID in the imported data file"))
    }
    subject_id <- includeIDs[includeIDs %in% colnames(rawData)]
    if (!is.null(excludeIDs)){
      subject_id <- subject_id[!subject_id %in% excludeIDs]
    }
  }

  cat("Included subject IDs: ")
  cat(subject_id, "\n")


  #extract individual genotype "GT" from genotype columns

  cat("Extract individual genotype [GT] from genotype columns \n ")
  data_sub_1 <- data_sub
  for (i in 1:length(subject_id)) {
    col_to_extract <- subject_id[i]
    col_GT <- paste0(subject_id[i], "_GT")
    data_sub_1[,col_GT] <- sapply(data_sub_1[,col_to_extract],FUN = function(x) {
      stringr::str_sub(x, start = 1, end = 3)
    })
  }

  cat("Individual genotypes saved in new columns: ")
  cat(colnames(data_sub_1)[(ncol(data_sub) + 1):ncol(data_sub_1)], "\n")

  # Remove variants with only "0/0" and "./." genotypes across all individuals
  data_sub_1 <- data_sub_1[
    !apply(data_sub_1[,paste0(subject_id, "_GT"), drop=FALSE],1, FUN = function(x) {
      all(x %in% c( "0/0","./."))
    })
    ,]

  #Extract [QD] from INFO column

  data_sub_3 <- data_sub_1

  data_sub_3$QD_fromINFO <- sapply(data_sub_1$INFO, FUN=function(x) {
    string_todetect<- x
    length_string <- stringr::str_length(string_todetect)
    QD_INFO_index <- stringr::str_locate( string_todetect, "QD")
    endpoint_index <- QD_INFO_index[1] - 1 +  stringr::str_locate(stringr::str_sub(string_todetect, start=QD_INFO_index[1], end = length_string), ";")[1] -1
    startpoint_index <- QD_INFO_index[1] + 3
    QD_number <- stringr::str_sub(string_todetect, start = startpoint_index, end = endpoint_index)
    return(QD_number)
  })
  data_sub_3$QD_fromINFO <- as.numeric(data_sub_3$QD_fromINFO)
  cat("QD extracted from INFO column and saved in a new column: QD_fromINFO", "\n")

  cat(paste0("Total number of variants in raw data: ", nrow(data_sub_3)), "\n")

  if (removeMultiallelicSites) {
    data_sub_3 <- data_sub_3[!unlist(sapply(data_sub_3[,"ALT"],FUN=function(alt) grepl(",", alt))),]
    n_multiallelic <- nrow(data_sub_3)

    cat("Number of variants left after removing multiallelic sites: ", n_multiallelic, "\n")

  }

  if (QD.threshold > 0) {
    data_sub_3 <- data_sub_3[as.numeric(data_sub_3$QD_fromINFO) >= QD.threshold,]
    n_QD <- nrow(data_sub_3)
    cat("Number of variants left after removing variants with QD lower than ", QD.threshold, ": ", n_QD, "\n")
  }


  if (remove.nocall.pct >0 ) {
    data_sub_3 <- data_sub_3[apply(data_sub_3[, c(paste0(subject_id,"_GT")), drop = FALSE], 1, FUN= function(row) mean(row == "./.")) < remove.nocall.pct, ]
    n_remove.nocall <- nrow(data_sub_3)
    cat("Number of variants left after removing variants with percentage of no calls more than ", remove.nocall.pct, ": ", n_remove.nocall, "\n")
  }

  #save data
  if (saveRData) {
    assign(label, data_sub_3)
    if (is.null(savePath)) {

    savePath <- paste0(getwd(),"/")
    warning(paste0("Save path not specified. Data will be saved at the working directory ",getwd()))


    }
    suppressWarnings(dir.create(savePath))
    save(list = paste0(label), file = paste0(savePath, label, ".RData"))
    cat("Saved object path: ")
    cat(paste0(savePath, label, ".RData"), "\n")
  }
  return(data_sub_3)
}
