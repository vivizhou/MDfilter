#' Generate pedigree plots
#'
#' This function generates pedigree plots to visualize variant information and individual genotypes.
#'
#' @param variant.list Filtered variant list from ?filter.variant()
#' @param includeIDs IDs of the subjects included for analyses. This must match IDs used in the imported data file.If left empty, IDs will be extracted from the first to the last
#' @param gene.name specify the list of genes in which the variants are to be visualized in the pedigree plots
#' @param pedigree.data Pedigree data load from excel sheet. check ?kinship2::pedigree() for correct format.
#' @param frequency.col Column names of the frequencies to be visualized in pedigree plot. For example: frequency.col = c("ExAC_All", "X1000G_ALL", "ExAC_NFE", "X1000G_EUR")
#' @param save.path Specify save path. If not specified the data will be saved at the working directory.
#' @param label Specify the label of the folder to save the files.

#' @return Pedigree plot saved as PNG files.
#'
#' @export
#'
plot.pedigreeGT <- function(variant.list,
                            pedigree.data,
                            includeIDs,
                            gene.name = NULL,
                            save.path = NULL,
                            frequency.col = c("ExAC_All", "X1000G_ALL"),
                            label = "test") {
  if (!is.null(save.path)) {
    save_path <- save.path
  }else{
    save_path <- paste0(getwd(), "/")
    warning(paste0("Save path not specified. Data will be saved at the working directory ",getwd()))
  }

  path_name <- paste0(save.path, "/",label)
  if(!dir.exists(path_name)){
    dir.create(path_name)
  }

  if (!is.null(gene.name)) {
    variant.list <- variant.list[variant.list$Gene.refgene %in% gene.name,]
  }

  for (i in (1:nrow(variant.list))) {
  variant <- variant.list[i,]
  variant_variant <- variant$Variant
  variant_pos <- paste0(variant$Chr,"_", variant$Start)
  variant_mut <- paste0("Ref: ",variant$Ref, "; Alt: ", variant$Alt)
  variant_gene <- paste0(variant$Gene.refgene)
  variant_AAchange <- paste0(variant$AAChange.refgene)
  variant_frequency <- c()
  for (freqSource in frequency.col) {
    variant_frequency <- c(variant_frequency,
                           paste0(freqSource,": ",variant[paste0(freqSource)]))
  }
  variant_frequency <- paste0(variant_frequency, collapse = ";")
  variant_CADD <- variant$CADD_phred
  variant_Func1 <- paste0(variant$Func.refgene)
  variant_Func2 <- paste0(variant$ExonicFunc.refgene.VI)
  GT_info <- as.data.frame(t(as.vector((variant[paste0(includeIDs, "_GT")]))))
  GT <- data.frame(id = (includeIDs),
                   GT = ((GT_info)),
                   stringsAsFactors = FALSE)
  colnames(GT)[2] <- "GT"
  pedpre_GT <- dplyr::left_join(pedigree.data, GT, by = "id")
  pedpre_GT$father[pedpre_GT$father == 0] <- NA
  pedpre_GT$mother[pedpre_GT$mother == 0] <- NA
  ped <- kinship2::pedigree(id = pedpre_GT$id,
                  dadid = pedpre_GT$father,
                  momid = pedpre_GT$mother,
                  sex = pedpre_GT$gender,
                  status = pedpre_GT$status,
                  as.matrix(data.frame(Affected = pedpre_GT$affected, DNA = pedpre_GT$avail)))
  strid <- paste(pedpre_GT$id, pedpre_GT$GT, sep = "\n")
  png(paste0(path_name, "/", variant_pos, ".png"), width = 786, height = 551)#, width = 1350, height = 580
  plot.pedigree(ped, id = strid, mar = c(12, 2.1, 2.1, 2.1))

  mtext(text = paste0(variant_v, "\n",
                        variant_pos, ";", variant_mut, "\n",
                        variant_gene, "\n",
                        "Frequencies::",variant_frequency, "\n",
                        "CADD=",variant_CADD, "\n",
                        "Func.refgene: ", variant_Func1,", ", "ExonicFunc: ", variant_Func2, "\n",
                        "AAchange:", variant_AAchange),
                        side = 1, at = 0, adj=0)
  kinship2::pedigree.legend(ped, location = "bottomright", radius = .1 )
  dev.off()

  }


}
