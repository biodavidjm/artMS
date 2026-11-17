# Everything about enrichments (including plots)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# @title Enrich for Protein Complexes using CORUM
#
# @description Enrich for Protein Complexes using CORUM
# @param df (data.frame) Data.frame with columns: `Protein` and `Conditions`
# @param backgroundNumber (int) Background number of genes
# @return (data.frame) A text delimited data.frame with protein complex 
# enrichment results
# @keywords internal, enrichment, protein, complexes
.artms_enrichForComplexes <- function(df, 
                                      backgroundNumber) {
  
  if(any(missing(df) | missing(backgroundNumber)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  listOfConditions <- unique(df$Comparisons)
  
  complexEnrichmentConditions <- data.frame()
  
  for (i in seq_len(length(listOfConditions))) {
    condi <- listOfConditions[i]
    tmp <- unique(df$Protein[which(df$Comparisons == condi)])
    tmpEnrich <-
      .artms_foldComplexEnrichment(tmp, 
                                   artms_data_corum_mito_database, 
                                   backgroundNumber)
    # Check point
    checkpc <- dim(tmpEnrich)[1]
    if (checkpc > 0) {
      tmpEnrich$Comparisons <- condi
      complexEnrichmentConditions <-
        rbind(complexEnrichmentConditions, tmpEnrich)
    }
  }
  return(complexEnrichmentConditions)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Enrichment of changes in protein abundance or PTMs
#'
#' @description Enrichment analysis of the selected proteins
#' @param dataset (data.frame) with a `Gene` and `Comparison or Label` (with
#' the name of the comparisons specified in the contrast file) columns
#' @param species (char) Specie, only supported "human" or "mouse"
#' @param background (vector) Background genes for the enrichment analysis.
#' @param heatmaps (logical) if `TRUE` generates heatmaps (pdf),
#' `FALSE` (default) otherwise.
#' @param output_name (char) Name of the annotation files, which will be used
#' as well for the heatmaps (if `heatmaps` is selected)
#' Default `output_name = "enrichment.txt"`
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame) Results from the enrichment analysis using Gprofiler
#' and heatmaps (if selected)
#' @keywords enrichment
#' @examples \dontrun{
#' # The data must be annotated (Protein and Gene columns)
#' data_annotated <- artmsAnnotationUniprot(
#'                       x = artms_data_ph_msstats_results,
#'                       columnid = "Protein",
#'                       species = "human")
#' # And then the enrichment
#' enrich_set <- artmsEnrichLog2fc(
#'                    dataset = data_annotated,
#'                    species = "human",
#'                    background = unique(data_annotated$Gene))
#'}
#' @export
artmsEnrichLog2fc <- function(dataset,
                              species,
                              background,
                              heatmaps = FALSE,
                              output_name = "enrichment.txt",
                              verbose = TRUE) {
                               
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package <gprofiler2> required to run this function. ",
         "Please, install it (install.packages('gprofiler2')) and run it again")
  }

  if(any(missing(dataset) |
         missing(species) |
         missing(background)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")

  if (heatmaps) {
    if (!grepl("\\.txt", output_name)) {
      stop("The <output_name> argument does not have extension '.txt'")
    }
  }
  
  # First check the "comparisons" column. if it is a msstats results file
  # should have instead the "Label" column.
  if (!any(grepl("Comparisons", colnames(dataset)))) {
    if (any(grepl("Label", colnames(dataset)))) {
      dataset <- artmsChangeColumnName(dataset, "Label", "Comparisons")
    } else{
      stop("This dataset does not have a <Label> or <Comparisons> column")
    }
  }
  
  # Selecting unique genes in each comparison
  pretmp <- dataset[c('Gene', 'Comparisons')]
  pretmp <- unique(pretmp)
  
  tmp = split(pretmp$Gene, pretmp$Comparisons, drop = TRUE)
  
  if (species == "human") {
    enrichgenes <- artmsEnrichProfiler(x = tmp,
                                       categorySource = c(
                                         'GO:BP',
                                         'GO:MF',
                                         'GO:CC',
                                         'KEGG',
                                         'REAC',
                                         'CORUM',
                                         'HPA'
                                       ),
                                       species = 'hsapiens',
                                       background = background,
                                       verbose = verbose) # 'HP'
  } else if (species == "mouse") {
    enrichgenes <-
      artmsEnrichProfiler(
        tmp,
        categorySource = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'CORUM'),
        species = 'mmusculus',
        background, verbose = verbose
      )
  } else{
    stop(" The species (",
         species,
         ") not supported for the enrichment analysis!! ")
  }
  
  if (dim(enrichgenes)[1] == 0) {
    message("--- No significant results from the enrichment analysis ")
  } else{
    if (is.null(heatmaps)) {
      enrichgenes2plot <- enrichgenes[which(enrichgenes$term.size < 500), ]
      for (i in unique(enrichgenes2plot$domain)) {
        if(verbose) message("\t--- Plotting '", i, "' annotations... ")
        tmp <-
          enrichgenes2plot[which(enrichgenes2plot$domain == i), ]
        outfile <- gsub(".txt", paste0("_", i, ".txt"), output_name)
        .artms_EnrichmentPlotHeatmaps(tmp, outfile)
        if(verbose) message("out! ")
      }
    }
  }
  return(enrichgenes)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# @title Plot Corum Enrichment results
#
# @description Heatmap of significant enrichment
# @param x (data.frame) output from `artms_enrichForComplexes`
# @param outfile (char) output file name (must have the extenstion `.pdf`)
# @param theTitle (char) Plot's title
# @return (pdf) heatmap of the significantly enriched protein complexes
# @keywords internal, plot, heatmap, enrichment
.artms_plotCorumEnrichment <- function(x, 
                                       outfile, 
                                       theTitle) {
  Comparisons = p_value = var = NULL

  if(any(missing(x) | 
         missing(outfile) |
         missing(theTitle)))
    stop("Missed (one or many) required argument(s)
        Please, check the help of this function to find out more")
  
  checkPoint <- length(unique(x$Comparisons))
  
  if (checkPoint >= 1) {
    # Dealing with the smallest pvalues (pvalue = 0)
    dftemp <- x[-which(x$pvalue == 0), ]
    if (dim(dftemp)[1] > 0) {
      theMinimal <- min(dftemp$pvalue) / 10
      # Make the value replacement
      x$p_value <- x$pvalue
      x$p_value[x$p_value == 0] <- theMinimal
    } else{
      x$p_value <- x$pvalue
    }
    
    x$p_value <- -log10(x$p_value)
    
    ##LEGACY
    # toplot <- data.table::dcast(data = x,
    #                             ComplexName ~ Comparisons,
    #                             value.var = "p_value",
    #                             fun.aggregate = sum,
    #                             fill = 0)
    
    toplot <- x %>% 
      tidyr::pivot_wider(id_cols = ComplexName, 
                         names_from = Comparisons, 
                         values_from = p_value, 
                         values_fn = list(p_value = sum), 
                         values_fill = list(p_value = 0))

    toplot <- as.data.frame(toplot)
    
    rownames(toplot) <- toplot$ComplexName
    toplotmatrix <- subset(toplot, select = -c(ComplexName))

    xmatrix <- data.matrix(toplotmatrix)
    
    # HEATMAP
    # palette.breaks <- seq(1, 4, 0.1)
    # color.palette  <- grDevices::colorRampPalette(c("white", "steelblue"))(length(palette.breaks))
    
    if(var(xmatrix[,1]) == 0){
      message("----Heatmap of Corum enrichment is not possible")
    }else{
      pheatmap::pheatmap(
        xmatrix,
        filename = outfile,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        cellheight = 10,
        cellwidth = 25,
        main = theTitle,
        fontsize = 6,
        fontsize_row = 8,
        fontsize_col = 12,
        border_color = 'black',
        fontfamily = "Helvetica",
        treeheight_row = FALSE,
        treeheight_col = FALSE,
        # color = color.palette
      )
    }

  } else{
    message("---(-) Not enough enriched comparisons to plot the heatmap ")
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Enrichment analysis using g:Profiler (gprofiler2)
#'
#' @description This function simplifies the enrichment analysis performed by
#' the excellent tool g:Profiler, using the `gost()` function of the
#' `gprofiler2` package (which replaces the deprecated `gProfileR` package).
#' An internet connection is required.
#' @param x (list, data.frame) List of protein ids. It can be anything:
#' either a list of ids, or you could also send a data.frame and it will find
#' the columns with the IDs. Is not cool? Multiple list can be also sent
#' simultaneously, as for example running:
#' `tmp <- split(enrichment$Gene, enrichment$cl_number, drop= TRUE)`
#' @param categorySource (vector) Resources providing the terms on which
#' the enrichment will be performed. The resources supported by g:Profiler are:
#' - GO (GO:BP, GO:MF, GO:CC): Gene Ontology. `GO` selects the three of them
#' - KEGG: Biological pathways
#' - REAC: Biological pathways (Reactome)
#' - WP: Biological pathways (WikiPathways)
#' - TF: Regulatory motifs in DNA (TRANSFAC TFBS)
#' - MIRNA: Regulatory motifs in DNA (miRTarBase microRNAs)
#' - CORUM: protein complexes database
#' - HP: Human Phenotype Ontology
#' - HPA: Protein databases (Human Protein Atlas)
#' `OMIM` is no longer available in g:Profiler and it is ignored.
#' Electronic GO annotations (IEA) are excluded.
#' @param species (char) Specie code: Organism names are constructed by
#' concatenating the first letter of the name and the family name.
#' Example: human - ’hsapiens’, mouse - ’mmusculus’. Check g:Profiler to find
#' out more about supported species.
#' @param background (vector) gene list to use as background for the enrichment
#' analysis. Default: `NA` (no custom background: all known genes are used)
#' @details This function calls `gprofiler2::gost()` with these arguments:
#' - ordered_query = FALSE
#' - significant = TRUE
#' - exclude_iea = TRUE
#' - measure_underrepresentation = FALSE
#' - evcodes = TRUE (needed to report the `intersection` column)
#' - user_threshold = 0.05
#' - correction_method = "g_SCS" (named "analytical" in gProfileR)
#' - domain_scope = "known", or "custom" when a `background` is provided
#'
#' The column names of the results are kept as in previous versions of artMS
#' (based on `gProfileR`): `query.number`, `significant`, `p.value`,
#' `term.size`, `query.size`, `overlap.size`, `recall`, `precision`,
#' `term.id`, `domain`, `term.name`, `intersection`.
#' Note that `domain` now uses the g:Profiler source names (e.g. `GO:BP`,
#' `KEGG`, `REAC`).
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame) The significant enrichment results (zero rows if
#' there are none)
#' @keywords enrichment
#' @examples \dontrun{
#' # annotate the MSstats results to get the Gene name
#' data_annotated <- artmsAnnotationUniprot(
#'                                      x = artms_data_ph_msstats_results,
#'                                      columnid = "Protein",
#'                                      species = "human")
#'
#' # Filter the list of genes with a log2fc > 2
#' filtered_data <- 
#' unique(data_annotated$Gene[which(data_annotated$log2FC > 2)])
#'
#' # And perform enrichment analysis
#' data_annotated_enrich <- artmsEnrichProfiler(
#'                                    x = filtered_data,
#'                                    categorySource = c('KEGG'),
#'                                    species = "hsapiens",
#'                                    background = unique(data_annotated$Gene))
#'}
#' @export
artmsEnrichProfiler <- function(x,
                                categorySource = c('GO'),
                                species,
                                background = NA,
                                verbose = TRUE) {
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package <gprofiler2> required to run this function. ",
         "Please, install it (install.packages('gprofiler2')) and run it again")
  }

  if(any(missing(x) | missing(species)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")

  # 'GO' selects the three Gene Ontology branches
  if ("GO" %in% categorySource) {
    categorySource <- unique(c(setdiff(categorySource, "GO"),
                               "GO:BP", "GO:MF", "GO:CC"))
  }
  # OMIM is no longer a g:Profiler data source
  if ("OMIM" %in% categorySource) {
    if(verbose) message("---(!) 'OMIM' is no longer available in ",
                        "g:Profiler: ignored")
    categorySource <- setdiff(categorySource, "OMIM")
  }

  # No custom background unless the user provides one
  custom_bg <- NULL
  if (length(background) > 0 && !all(is.na(background))) {
    custom_bg <- unique(as.character(background[!is.na(background)]))
  }

  if(verbose) message("---+ Enrichment analysis using gProfiler...", appendLF = FALSE)
  scope <- if (is.null(custom_bg)) "known" else "custom"
  gp <- suppressWarnings(
    gprofiler2::gost(query = x,
                     organism = species,
                     ordered_query = FALSE,
                     multi_query = FALSE,
                     significant = TRUE,
                     exclude_iea = TRUE,
                     measure_underrepresentation = FALSE,
                     evcodes = TRUE,
                     user_threshold = 0.05,
                     correction_method = "g_SCS",
                     domain_scope = scope,
                     custom_bg = custom_bg,
                     sources = categorySource)
  )
  enrichData <- .artms_gostToGprofiler(gp$result)

  if(verbose) message("done! ")
  return(enrichData)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# @title Convert gprofiler2::gost() results to the gProfileR column names
#
# @description Earlier versions of artMS used `gProfileR::gprofiler()`. This
# function keeps those column names so that downstream code and the output
# files keep the same format. List columns (e.g. `parents`) are dropped
# because they cannot be written to text files.
# @param res (data.frame) the `result` element returned by `gost()`, or `NULL`
# when there are no significant results
# @return (data.frame) with the columns 'query.number', 'significant',
# 'p.value', 'term.size', 'query.size', 'overlap.size', 'recall', 'precision',
# 'term.id', 'domain', 'term.name', 'intersection' (zero rows if `res` is NULL)
# @keywords internal, enrichment
.artms_gostToGprofiler <- function(res) {
  legacy_names <- c("query.number", "significant", "p.value", "term.size",
                    "query.size", "overlap.size", "recall", "precision",
                    "term.id", "domain", "term.name", "intersection")
  gost_names <- c("query", "significant", "p_value", "term_size",
                  "query_size", "intersection_size", "recall", "precision",
                  "term_id", "source", "term_name", "intersection")
  if (is.null(res) || nrow(res) == 0) {
    empty <- as.data.frame(matrix(nrow = 0, ncol = length(legacy_names)))
    colnames(empty) <- legacy_names
    return(empty)
  }
  missing_cols <- setdiff(gost_names, colnames(res))
  if (length(missing_cols) > 0) {
    stop("Unexpected gprofiler2 output. Missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  out <- as.data.frame(res)[gost_names]
  colnames(out) <- legacy_names
  return(out)
}

# Little function to
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# @title Simplify the gProfiler output
#
# @description Simplify the output from `artmsEnrichProfiler` resulted from
# running `gProfileR`
# @param gp (data.frame) with the results
# @return (data.frame) with the following columns:
#       'query.number', 'domain', 'p.value', 'query.size', 'overlap.size',
#       'term.size', 'recall', 'precision', 'term.id', 'term.name',
#       'intersection'
# @keywords internal, cleaning
.artms_cleanGPROFILER <- function(gp) {
  sendBack <- gp[c(
    'query.number',
    'domain',
    'p.value',
    'query.size',
    'overlap.size',
    'term.size',
    'recall',
    'precision',
    'term.id',
    'term.name',
    'intersection'
  )]
  return(sendBack)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# @title Plot and save heatmaps of the significant enrichment results
#
# @description plot and save heatmaps of the significant enrichment results
# @param x (data.frame) output from gprofiler
# @param out_file (char) output file name (must have `.txt` extension)
# @return (pdf) A heatmap of the most significant enrichments
# @keywords internal, plot, heatmap, enrichments
.artms_EnrichmentPlotHeatmaps <- function(x, 
                                          out_file) {
  
  term.name = query.number = p.value = NULL
  
  ##LEGACY
  # x <- dcast(x,
  #            term.name ~ query.number,
  #            value.var = 'p.value',
  #            max,
  #            fill = 1)
  
  x <- x %>% tidyr::pivot_wider(id_cols = term.name, 
                                names_from = query.number, 
                                values_from = p.value,
                                values_fn = list(p.value = max),
                                values_fill = list(p.value = 1))
  
  # Let's stop this thing if there is not enough terms (we need at least 2)
  if (dim(x)[1] < 2) {
    return(message(" Not enough terms in this domain "))
  } else{
    row.names(x) = x$term.name
    x$term.name = c()
    # At least one comparison ready
    if (dim(x)[2] > 0) {
      idx <- apply(x, 1, function(y) {
        return(min(y) <= .05)
      })
      x <- x[which(idx), ]
      x <- -log10(x)
      
      # setting up breaks for colors of -log10 p-values
      HEAT_STYLE = "Q" ##HEAT_STYLE="SCORE"
      P_T = 0.05
      BREAKS = c(0, -log10(c(P_T, P_T / 10, P_T / 100, P_T / 1000)))
      LABELS = c(0, "  5E-2", "  5E-3", "< 5E-4", "")
      
      ## format stuff for heatmap
      colors = c(colorRampPalette(
        brewer.pal(n = 7, name = "Purples"))(length(BREAKS) - 1))
      LOWER_T = BREAKS[length(BREAKS)]
      term_groups_selected_w_display = x
      term_groups_selected_w_display[term_groups_selected_w_display > LOWER_T] =
        LOWER_T
      if (length(x) > 1) {
        #Do you need a main title: 
        #main=paste( gsub(".txt","",basename(out_file)) , 
        #"(color: -log10 p-value )"),
        pheatmap(
          term_groups_selected_w_display,
          cluster_cols = FALSE,
          cellheight = 10,
          cellwidth = 10,
          scale = "none",
          filename = gsub('.txt', '_heatmap.pdf', out_file),
          fontsize = 6,
          fontsize_row = 8,
          fontsize_col = 8,
          border_color = NA,
          color = colors,
          breaks = BREAKS,
          legend_breaks = BREAKS,
          legend_labels = LABELS,
          fontfamily = "Helvetica"
        )
      } else{
        message(" [artMS currently doesn't support heatmaps of a single set] ")
        #pheatmap(term_groups_selected_w_display, 
        #cluster_cols= FALSE,cluster_rows= FALSE, 
        #cellheight=10, cellwidth=10, scale="none", 
        #filename=gsub('.txt','_heatmap.pdf',out_file), 
        #fontsize=6, fontsize_row=8, fontsize_col=8,
        # border_color=NA, color = colors, 
        # fontfamily="Helvetica")
      }
    } else{
      message(" [Not enough significant terms for this domain] ")
    }
  }
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# @title Enrichment analysis of Protein Complexes (based on CORUM database)
#
# @description Enrichment analysis of Protein Complexes
# (based on CORUM database)
# @param mylist (vector) of protein ids
# @param corum (data.frame) The corum database (with the corum format
# and labels)
# @param background (int) Total number of proteins (number) to use as
# background
# @return (data.frame) with the list of protein complexes with a fold change
# larger than 1
# @keywords internal, enrichment, protein, complexes
# .artms_foldComplexEnrichment()
.artms_foldComplexEnrichment <- function(mylist, 
                                         corum, 
                                         background) {
  corum$Num.Uniprot.IDs <-
    vapply(corum$subunits.UniProt.IDs, function(x)
      length(unlist(strsplit(
        as.character(x), ";"
      ))), FUN.VALUE = 0)
  
  # Matches between my list and the protein complexes
  matchThis <- function(cor, mylist) {
    corv <- unlist(strsplit(as.character(cor), ";"))
    length(table(mylist[mylist %in% corv]))
  }
  corum$MyListOverlap <-
    vapply(corum$subunits.UniProt.IDs, function(x)
      matchThis(x, mylist), FUN.VALUE = 0)
  
  # Fold Complex Enrichment of the proteins observed in the list of proteins
  # over the expected. If it is greater than 1, it indicates that the category
  # is overrepresented in the experiment.
  # Conversely, the category is underrepresented if it is less than 1.
  
  # For the human complex enrichment we used this number for the background
  # numberHumanProteosome <- 8522
  # background <- numberHumanProteosome
  # This number is based on this reference:
  # Mol Cell Proteomics. 2012 Mar;11(3):M111.014050.
  # doi: 10.1074/mcp.M111.014050. Epub 2012 Jan 25.
  # Comparative proteomic analysis of eleven common cell lines reveals
  # ubiquitous but varying expression of most proteins.
  # Geiger T1, Wehner A, Schaab C, Cox J, Mann M.
  
  querysize <- length(unique(mylist)) #uploaded proteins
  corum$querysize <- querysize
  corum$expected <- (corum$Num.Uniprot.IDs / background) * querysize
  corum$FCEnrichment <- corum$MyListOverlap / corum$expected
  corum$pvalue <-
    phyper(
      corum$MyListOverlap,
      corum$Num.Uniprot.IDs,
      (background - corum$Num.Uniprot.IDs),
      querysize,
      lower.tail = FALSE
    )
  # Return only complexes with a FC > 1
  toreturn <-
    corum[which(corum$FCEnrichment >= 1 & corum$pvalue < 0.05), ]
  return(toreturn)
}
