#' Variant filtering
#'
#' This function filters the variants based on specified criteria and will print the filtering log on the console.
#'
#' @param data.name Data cleaned by ?clean.raw.data()
#' @param affected.id IDs of the affected subjects,which should match the IDs in the data file (column name)
#' @param unaffected.id IDs of the unaffected subjects,which should match the IDs in the data file (column name)
#' @param inheritance.pattern Model of inheritance pattern. #"AR" = autosomal recessive/"AD"=autosomal dominant/"XLD"=X-linked dominant/"XLR"=X-linked recessive. You can specify "AD" for de novo variant filtering. If XLR is specified, female and male IDs should also be provided.
#' @param femaleIDs Female and male IDs should also be provided.
#' @param maleIDs Female and male IDs should also be provided.
#' @param gene.name Specify the gene list to run candidate gene analysis.Example:gene.name = c("ALPL", "COL1A1", "PLS3")
#' @param frequency.col Specify the column names of the frequencies to be used for filtering. For example: frequency.col = c("ExAC_All", "X1000G_ALL", "ExAC_NFE", "X1000G_EUR")
#' @param freq.threshold The threshold of the frequency used for filtering.
#' @param protein.altering Only include protein altering variants.
#' @param include.UTRs On top of protein altering variants, TRUE=also include UTRs
#' @param include.synonymous On top of protein altering variants, TRUE=also to include synonymous variants
#' @param include.nonframeshift TRUE=also to include nonframeshift variants
#' @param CADD.threshold The threshold of the CADD score used for filtering
#' @param save.txt TRUE=Save the returned data frame into a txt file.
#' @param save.path Specify save path. If not specified the data will be saved at the working directory.
#' @param label Specify the label of the saved file.
#' @param save.genelist TRUE=Save the filtered gene list into a txt file.

#' @return A data frame of cleaned annotated data.
#'
#' @export
filterVariant <- function(data.name,
                           affected.id = NULL,
                           unaffected.id = NULL,
                           inheritance.pattern = NULL,
                           femaleIDs = NULL,
                           maleIDs = NULL, # specify female and male IDs if choose X-linked recessive(XLR)
                           gene.name = NULL, # specify gene list for candidate gene analysis
                           frequency.col = NULL,
                           freq.threshold = 1,
                           protein.altering = FALSE, # if true only include exonic, splicing, stopgain, stoploss; excluding synonymous variants and UTRs unless specified
                           include.UTRs = FALSE,
                           include.synonymous = FALSE,
                           include.nonframeshift = TRUE,
                           CADD.threshold = 0,
                           save.path = NULL,
                           save.txt = FALSE,
                           label = "test",
                           save.genelist = FALSE
                           ) {
  subdata <- data.name
  if (any(is.na(affected.id) | is.na(unaffected.id))) {
    stop("Please do not include missing values in affected/unaffected IDs.")
  }

  for (freqSource in frequency.col) {
    freqSourceIndex <- which(colnames(subdata) %in% freqSource)
    if (length(freqSourceIndex) < 1) {
      stop("Specified alelle frequency column ", freqSource, " does not exist in the data file.")
    }
  }

  if (inheritance.pattern == "XLR") {
    if (is.null(femaleIDs) & is.null(maleIDs)) {
      stop("Please specify male and female IDs when using X-linked recessive model")
    }
  }

  filterlog <- matrix("filtering log", ncol=1)
  subdata <- data.name
  n_start <- nrow(subdata)
  cat("Starting with number of variants: ", n_start, "\n")
  filterlog <- rbind(filterlog, paste0("Starting with number of variants: ", n_start))

  # Extract variant in the candidate genes
  if (!is.null(gene.name)) {
    #require data
    data(gene_coordinates)

    cat("Candidate genes selected:\n", paste0(gene.name, collapse = ","), " \n")
    filterlog <- rbind(filterlog, paste0("Candidate genes selected:\n", paste0(gene.name, collapse = ",")))

    #Obtain the gene coordinates of the candidate genes from the "gene_coordinates" file

    gene_coords <- gene_coordinates[grep(paste0(paste0("gene_name ", gene.name, ";"), collapse = "|"), ignore.case = TRUE, gene_coordinates$V9),]
    gene_coords <- unique(gene_coords[,c("V1", "V4", "V5")])

    #Find the variants with the coordinates

    subdata <- subdata[subdata$Chr %in% unique(gene_coords$V1),]
    subdata <- subdata[apply(subdata,1, FUN=function(x) {
      gene_coords_cut <- gene_coords[gene_coords$V1 %in% x["Chr"],, drop= FALSE]
      mark <- FALSE
      if (!((is.integer(as.numeric(x["End"]))|is.integer(as.numeric(x["Start"]))|
             is.na(as.numeric(x["End"]))|is.na(as.numeric(x["Start"]))))) {
      if (!(as.numeric(x["End"]) < gene_coords_cut$V4[1] | as.numeric(x["Start"]) > gene_coords_cut$V5[nrow(gene_coords_cut)])) {
        for ( i in 1:nrow(gene_coords_cut)) {
          if (!mark) {
            mark <- (as.numeric(x["Start"]) >= gene_coords_cut$V4[i] & as.numeric(x["Start"]) <= gene_coords_cut$V5[i]) | (as.numeric(x["End"]) >= gene_coords_cut$V4[i] & as.numeric(x["End"]) <= gene_coords_cut$V5[i])
          } else {
            break
          }
        }
      }}
      return(mark)
    }),]

    n_gene_filter <- nrow(subdata)

    cat("Number of variants in candidate gene list: ", n_gene_filter, "\n")
    filterlog <- rbind(filterlog, paste0("Number of variants in candidate gene list: ", n_gene_filter))

  }

  if (!is.null(inheritance.pattern)){

    cat("Variants will be filtered based on ", inheritance.pattern, " model\n")
    filterlog <- rbind(filterlog, paste0("Variants was filtered based on", inheritance.pattern, "model"))

    if (inheritance.pattern == "AD") {
      cat("variants are 0/1 in affected members including ", paste0(affected.id, collapse = ","),"\n")
      cat("Variants are 0/0 in unaffected members including ", paste0(unaffected.id, collapse = ","),"\n")
      filterlog <- rbind(filterlog, paste0("variants are 0/1 in affected members including ",paste0(affected.id, collapse = ",")))
      filterlog <- rbind(filterlog, paste0("variants are 0/0 in affected members including ",paste0(unaffected.id, collapse = ",")))

      if (!is.null(affected.id)) {
        subdata <- subdata[apply(subdata[,paste0(affected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% c("0/1"))),] #heterozygous so 1/1 removed
      }
      if(!is.null(unaffected.id)) {
        subdata <- subdata[apply(subdata[,paste0(unaffected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% "0/0")),]
      }
    }

    if (inheritance.pattern == "AR") {
      cat("variants are 1/1 in affected members including ", paste0(affected.id, collapse = ","),"\n")
      cat("Variants are 0/0 or 0/1 in unaffected members including ", paste0(unaffected.id, collapse = ","),"\n")
      filterlog <- rbind(filterlog, paste0("variants are 1/1 in affected members including ",paste0(affected.id, collapse = ",")))
      filterlog <- rbind(filterlog, paste0("variants are 0/0 or 0/1 in affected members including ",paste0(unaffected.id, collapse = ",")))

      if (!is.null(affected.id)) {
        subdata <- subdata[apply(subdata[,paste0(affected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% c("1/1"))),]
      }
      if(!is.null(unaffected.id)) {
        subdata <- subdata[apply(subdata[,paste0(unaffected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% c("0/0","0/1"))),]
      }
    }

    if (inheritance.pattern == "XLR"){
      subdata <- subdata[subdata$Chr == "chrX",]
      n_chrX <- nrow(subdata)
      cat("Number of variants in X chromosome: ", n_chrX, "\n")
      filterlog <- rbind(filterlog, paste0("Number of variants in X chromosome: ", n_chrX))

      if (!is.null(affected.id)) {
        female.affected.id <- NULL
        male.affected.id <- NULL
        if (!is.null(femaleIDs) & length(femaleIDs) >0) {
          female.affected.id <- femaleIDs[femaleIDs %in% affected.id]

        }
        if (!is.null(maleIDs) & length(maleIDs) >0) {
          male.affected.id <- maleIDs[maleIDs %in% affected.id]
        }

        if (!is.null(female.affected.id) & length(female.affected.id) >0) {
          cat("variants are 1/1 in affected female members including ", paste0(female.affected.id, collapse = ","),"\n")
          filterlog <- rbind(filterlog, paste0("variants are 1/1 in affected female members including ", paste0(female.affected.id, collapse = ",")))

          subdata <- subdata[apply(subdata[,paste0(female.affected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% c("1/1"))),] #heterozygous so 1/1 removed

        }

        if (!is.null(male.affected.id) & length(male.affected.id) >0) {
          cat("variants are 0/1 in affected male members including ", paste0(male.affected.id, collapse = ","),"\n")
          filterlog <- rbind(filterlog, paste0("variants are 0/1 in affected male members including ", paste0(male.affected.id, collapse = ",")))
          subdata <- subdata[apply(subdata[,paste0(male.affected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% c("0/1"))),] #heterozygous so 1/1 removed

        }
      }

      if (!is.null(unaffected.id)) {
        female.unaffected.id <- NULL
        male.unaffected.id <- NULL
        if (!is.null(femaleIDs) & length(femaleIDs) >0) {
          female.unaffected.id <- femaleIDs[femaleIDs %in% unaffected.id]

        }
        if (!is.null(maleIDs) & length(maleIDs) >0) {
          male.unffected.id <- maleIDs[maleIDs %in% unaffected.id]
        }


        if (!is.null(female.unaffected.id) & length(female.unaffected.id) >0) {
          cat("variants are 0/0 or 0/1 in affected female members including ", paste0(female.unaffected.id, collapse = ","),"\n")
          filterlog <- rbind(filterlog, paste0("variants are 0/0 or 0/1 in affected female members including ", paste0(female.unaffected.id, collapse = ",")))
          subdata <- subdata[apply(subdata[,paste0(female.unaffected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% c("0/0", "0/1"))),]

        }

        if (!is.null(male.unaffected.id) & length(male.unaffected.id) >0) {

        cat("variants are 0/0 in affected male members including ", paste0(male.unffected.id, collapse = ","),"\n")

        filterlog <- rbind(filterlog, paste0("variants are 0/0 in affected male members including ", paste0(male.unffected.id, collapse = ",")))

        subdata <- subdata[apply(subdata[,paste0(male.unaffected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% c("0/0"))),]
        }

      }
    }

    if (inheritance.pattern == "XLD") {
      subdata <- subdata[subdata$Chr == "chrX",]
      n_chrX <- nrow(subdata)
      cat("Number of variants in X chromosome: ", n_chrX, "\n")
      filterlog <- rbind(filterlog, paste0("Number of variants in X chromosome: ", n_chrX))


      if (!is.null(affected.id)) {
        cat("variants are 0/1 or 1/1 in affected members including ", paste0(affected.id, collapse = ","),"\n")
        filterlog <- rbind(filterlog, paste0("variants are 0/1 or 1/1 in affected members including ",paste0(affected.id, collapse = ",")))

        subdata <- subdata[apply(subdata[,paste0(affected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% c("0/1", "1/1"))),] #heterozygous so 1/1 removed

      }
      if(!is.null(unaffected.id)) {
        cat("variants are 0/0 in affected female members including ", paste0(unaffected.id, collapse = ","),"\n")
        filterlog <- rbind(filterlog, paste0("variants are 0/0 in affected members including ",paste0(unaffected.id, collapse = ",")))

        subdata <- subdata[apply(subdata[,paste0(unaffected.id, "_GT"), drop = FALSE],1, FUN= function(row) all(row %in% "0/0")),]
      }
    }

    n_segregation <- nrow(subdata)
    cat("Number of variants left after filtering on segregation info: ", n_segregation, "\n")

  }


  if (protein.altering) {
    func_include_UTR <- c("exonic", "exonic;splicing", "splicing", "", "UTR5", "UTR3", "UTR5;UTR3")
    subdata <- subdata[subdata$Func.refgene %in% func_include_UTR,]
    n_func_include_UTR <- nrow(subdata)
    cat("Number of exonic, splicing and UTR variants: " ,n_func_include_UTR, "\n")
    filterlog <- rbind(filterlog, paste0("Number of exonic, splicing and UTR variants: " ,n_func_include_UTR))


    # Include/exclude synonymous
    if (include.synonymous) {
      cat("Synonymous variants are included.\n")
      filterlog <- rbind(filterlog, paste0("Synonymous variants are included."))

    } else {
      non_splicing_synonymous <- !(subdata$ExonicFunc.refgene %in% c("synonymous SNV")) | subdata$Func.refgene %in% c("splicing", "exonic;splicing")
      subdata <- subdata[non_splicing_synonymous,]
      n_synonymous_rmd <- nrow(subdata)
      cat("Number of variants after removing non-splicing synonymous variants: ", n_synonymous_rmd, "\n")
      filterlog <- rbind(filterlog, paste0("Number of variants after removing non-splicing synonymous variants: ", n_synonymous_rmd))

    }

    # Include/exclude UTRs
    if(include.UTRs){
      cat("UTRs are included.\n")
      filterlog <- rbind(filterlog, paste0("UTRs are included."))

    } else {
      Func_include <- c("exonic", "exonic;splicing", "splicing", "")
      subdata <- subdata[subdata$Func.refgene %in% Func_include,]
      n_functional <- nrow(subdata)
      cat("Number of variants after removing UTRs: ", n_functional, "\n")
      filterlog <- rbind(filterlog, paste0("Number of variants after removing UTRs: ", n_functional))

    }

    if(include.nonframeshift) {
      cat("Nonframeshift indels if present are kept.\n")
      filterlog <- rbind(filterlog, paste0("Nonframeshift indels if present are kept."))

    } else {
      subdata <- subdata[!subdata$ExonicFunc.refGene.3bp.splice %in% c("nonframeshift"),]
      n_nonframeshift_rmd <- nrow(subdata)
      cat("Number of variants left after removing nonframeshift indels: ", n_nonframeshift_rmd, "\n")
      filterlog <- rbind(filterlog, paste0("Number of variants left after removing nonframeshift indels: ", n_nonframeshift_rmd))

    }

  }
  if (!is.null(frequency.col)) {
    for (freqSource in frequency.col) {
      freqSourceIndex <- which(colnames(subdata) %in% freqSource)
      subdata <- subdata[subdata[,freqSourceIndex] <= freq.threshold
                         | is.na(subdata[,freqSourceIndex]),]
      n_freq_filtered <- nrow(subdata)
      cat("Number of variants left after filtering allele frequencies on ", freqSource, " with frequency cutoff of ", freq.threshold,": ", n_freq_filtered, "\n")
      filterlog <- rbind(filterlog, paste0("Number of variants left after filtering allele frequencies on ", freqSource, " with frequency cutoff of ", freq.threshold,": ", n_freq_filtered))

    }
  }


  if (CADD.threshold > 0) {
    cat("Removing variants with CADD score smaller than ", CADD.threshold, "\n")
    subdata <- subdata[subdata$CADD_phred >= CADD.threshold | is.na(subdata$CADD_phred),]
    n_cadd <- nrow(subdata)
    cat("Number of variants left: ", n_cadd, "\n")
    filterlog <- rbind(filterlog, paste0("Number of variants left: ", n_cadd))

  }





  if (save.txt) {
    if (!is.null(save.path)) {
      save_path <- save.path
    }else{
      save_path <- paste0(getwd(), "/")
      warning(paste0("Save path not specified. Data will be saved at the working directory ",getwd()))
    }
    if(!dir.exists(save_path)){
      dir.create(save_path, recursive = TRUE)
    }
    cat("A subset of variants saved at: ", paste0(save_path, "Variants_list_", label, ".txt"), ".\n")
    write.table(subdata, file = paste0(save_path, "Variants_list_", label, ".txt"), sep = "\t", row.names = FALSE)

    cat("The filtering log saved at: ", paste0(save_path, "Filtering_log_",label, ".txt"), ".\n")
    write.table(filterlog, file = paste0(save_path, "Filtering_log_",label, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

  }

  if (save.genelist) {
    if (!is.null(save.path)) {
      save_path <- save.path
    }else{
      save_path <- paste0(getwd(), "/")
      warning(paste0("Save path not specified. Data will be saved at the working directory ",getwd()))
    }

    if(!dir.exists(save_path)){
      dir.create(save_path, recursive = TRUE)
    }
    cat("The gene list saved at: ", paste0(save_path, "Genelist_",label, ".txt"))
    genelist <- unique(subdata[,"Gene.refgene", drop = FALSE])
    write.table(genelist, file = paste0(save_path, "Genelist_",label, ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE,
                quote = FALSE)
  }
  return(subdata)
}
