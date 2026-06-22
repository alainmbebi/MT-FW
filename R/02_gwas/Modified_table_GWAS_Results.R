#This is a modified code from GAPIT to include the exploratory threshold

#' Modified_table_GWAS_Results
#'
#' Create a table of significant GWAS results.
#' @param folder Folder containing GWAS results.
#' @param fnames The files to read.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param expl.threshold Exploratory threshold (e.g. 4.0).
#' @param nrowstoread Number of rows to read.
#' @param useHBPvalues Logical, if TRUE, H.B.P.Values will be used.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas".
#' @return A table of significant GWAS results.
#' @export
Modified_table_GWAS_Results <- function(
    folder         = "GWAS_Results/",
    fnames         = list_Result_Files(folder),
    nrowstoread    = 1000,
    threshold      = 6,
    sug.threshold  = NULL,
    expl.threshold = NULL,
    skyline        = NULL
) {

  output <- NULL

  for (i in fnames) {

    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1] + 13,
                    gregexpr(".csv", i)[[1]][1] - 1)
    trait <- gsub("\\(Kansas\\)|\\(NYC\\)", "", trait)
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1] - 1)
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1] + 1, nchar(trait))
    sky   <- substr(i, gregexpr("\\(", i)[[1]][1] + 1,
                    gregexpr("\\)", i)[[1]][1] - 1)

    oi <- read.csv(paste0(folder, i), nrows = nrowstoread) %>%
      mutate(
        Model        = model,
        Type         = sky,
        Trait        = trait,
        negLog10_P   = -log10(P.value),
        negLog10_HBP = -log10(H.B.P.Value)
      )

    # ----------------------------------------------------------
    # Classify SNPs into Significant / Suggestive / Exploratory
    # Priority: Significant > Suggestive > Exploratory
    # ----------------------------------------------------------
    oi <- oi %>%
      mutate(
        Threshold = case_when(
          negLog10_P >= threshold                                          ~ "Significant",
          !is.null(sug.threshold)  & negLog10_P >= sug.threshold  &
            negLog10_P < threshold                                         ~ "Suggestive",
          !is.null(expl.threshold) & negLog10_P >= expl.threshold &
            (is.null(sug.threshold) | negLog10_P < sug.threshold)         ~ "Exploratory",
          TRUE                                                             ~ "Suggestive"
        )
      )

    # ----------------------------------------------------------
    # Filter rows to retain based on lowest active threshold
    # ----------------------------------------------------------
    lowest_threshold <- dplyr::case_when(
      !is.null(expl.threshold) ~ expl.threshold,
      !is.null(sug.threshold)  ~ sug.threshold,
      TRUE                     ~ threshold
    )

    oi <- oi %>% filter(negLog10_P >= lowest_threshold)

    output <- bind_rows(output, oi)
  }

  # Remove nobs column if present
  if (sum(colnames(output) == "nobs") > 0) {
    output <- dplyr::select(output, -nobs)
  }

  # Skyline filter
  if (!is.null(skyline)) {
    if (skyline == "NYC") {
      output <- output %>%
        filter(!paste(Model, Type) %in% c("FarmCPU Kansas", "BLINK Kansas"))
    }
    if (skyline == "Kansas") {
      output <- output %>%
        filter(!paste(Model, Type) %in% c("FarmCPU NYC", "BLINK NYC"))
    }
  }

  # Remove duplicates and sort
  output <- output %>%
    arrange(desc(P.value)) %>%
    filter(!duplicated(paste(SNP, Model, P.value))) %>%
    # Order Threshold factor for clean downstream use
    mutate(
      Threshold = factor(Threshold,
                         levels = c("Significant", "Suggestive", "Exploratory"))
    ) %>%
    arrange(desc(negLog10_P))

  output
}
