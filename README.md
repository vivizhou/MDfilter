# MDfilter
## Title: Mendelian disease variant filtering
Version: 0.0.0.1
Author: Wei Zhou

The package provides a few functions for family-based genetic variant filtering to identify pathogenic variant for Mendelian diseases.


Install the package:
```r
devtools::install_github("vivizhou/MDfilter")
library(MDfilter)
?importTXT
?cleanData
?filterVariant
?plotPedigree
```


### Example usage:
#### Import txt file from ANNOVAR annotation on vcf files

```{r, results='hide'}
file.path <- "data.vcf.full_annotation.txt"
rawData <- importTXT(file.path)
```
#### Import pedigree file

```{r, results='hide'}
pedigreePlot <- readxl::read_excel("pedigree.xlsx", na = "NA", sheet = 1)
head(pedigreePlot)

```

Get included IDs (match the variant data file), affected IDs, and unaffected IDs. Also specify female and male IDs if X-linked recessive is assumed. 

```{r, results='hide'}
familyIDs <- pedigreePlot$id
includedSubjectIDs <- familyIDs[familyIDs %in% colnames(rawData)]
affectedIDs <- pedigreePlot[(pedigreePlot$affected == 1) & !is.na(pedigreePlot$affected), "id", drop = TRUE]
paste0("Affected IDs: ", affectedIDs)
unaffectedIDs <- pedigreePlot[pedigreePlot$affected == 0 & !is.na(pedigreePlot$affected == 0), "id", drop = TRUE]
paste0("Unaffected IDs: ", unaffectedIDs)

```

#### Clean imported data file

```{r, results='hide'}
cleanedData <- cleanData(rawData, includeIDs = includedSubjectIDs,
                              removeMultiallelicSites = FALSE,
                              remove.nocall.pct = 0.2,
                              QD.threshold = 5,
                              saveRData = FALSE, savePath = "DataExample/",
                              label = "testCleanedData")

```


#### Filter the variants

```{r, results='hide'}
filteredVariants <- filterVariant(data.name = cleanedData,
                                   affected.id = affectedIDs,
                                   unaffected.id = unaffectedIDs,
                                   inheritance.pattern = "dominant", 
                                   femaleIDs = NULL,
                                   maleIDs = NULL, 
                                   gene.name = NULL,
                                   frequency.col = c("ExAC_All", "X1000G_ALL", "ExAC_NFE", "X1000G_EUR"),
                                   freq.threshold = 0.005,
                                   protein.altering = TRUE, 
                                   include.UTRs = FALSE,
                                   include.synonymous = FALSE,
                                   include.nonframeshift = FALSE,
                                   CADD.threshold = 15,
                                   save.path = "DataExample/",
                                   save.txt = TRUE,
                                   label = "test",
                                   save.genelist = TRUE)
```

#### Plot pedigree file for APC gene

```{r, results='markup'}
plotPedigree(variant.list = filteredVariants,
              pedigree.data = pedigreePlot,
              includeIDs = includedSubjectIDs,
              gene.name = "APC",
              save.path = "DataExample",
              savePNG = TRUE,
              label = "testPlot")
```
