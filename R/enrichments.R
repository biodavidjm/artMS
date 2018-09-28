# Everything about enrichments (including plots)

# ------------------------------------------------------------------------------
# @title Enrich for Protein Complexes using CORUM
#
# @description Enrich for Protein Complexes using CORUM
# @param df (data.frame) Data.frame with columns: `Protein` and `Conditions`
# @param backgroundNumber (int) Background number of genes
# @return (data.frame) A text delimited data.frame with protein complex 
# enrichment results
# @keywords internal, enrichment, protein, complexes
.artms_enrichForComplexes <- function(df, backgroundNumber) {
  listOfConditions <- unique(df$Comparisons)
  
  complexEnrichmentConditions <- NULL
  
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


# ------------------------------------------------------------------------------
#' @title Enrichment of changes in protein abundance or PTMs
#'
#' @description Enrichment analysis of the selected proteins
#' @param dataset (data.frame) with a `Gene` and `Comparison or Label` (with
#' the name of the comparisons specified in the contrast file) columns
#' @param specie (char) Specie, only supported "human" or "mouse"
#' @param background (vector) Background genes for the enrichment analysis.
#' @param heatmaps (logical) if `TRUE` generates heatmaps (pdf),
#' `FALSE` (default) otherwise.
#' @param output_name (char) Name of the annotation files, which will be used
#' as well for the heatmaps (if `heatmaps` is selected)
#' Default `output_name = "enrichment.txt"`
#' @return (data.frame) Results from the enrichment analysis using Gprofiler
#' and heatmaps (if selected)
#' @keywords enrichment
#' @examples
#' # The data must be annotated (Protein and Gene columns)
#' data_annotated <- artms_annotationUniprot(
#'                       data = artms_data_ph_msstats_results,
#'                       columnid = "Protein",
#'                       sps = "human")
#' # And then the enrichment
#' enrich_set <- artms_enrichLog2fc(
#'                    dataset = data_annotated,
#'                    specie = "human",
#'                    background = unique(data_annotated$Gene))
#' @export
artms_enrichLog2fc <- function(dataset,
                               specie,
                               background,
                               heatmaps = FALSE,
                               output_name = "enrichment.txt") {
  # Check point
  if (heatmaps) {
    if (!grepl("\\.txt", output_name)) {
      stop("WRONG output_name VALUE: DOES NOT HAVE EXTENSION .txt)")
    }
  }
  
  # First check the "comparisons" column. if it is a msstats results file
  # should have instead the "Label" column.
  if (!any(grepl("Comparisons", colnames(dataset)))) {
    if (any(grepl("Label", colnames(dataset)))) {
      dataset <- artms_changeColumnName(dataset, "Label", "Comparisons")
    } else{
      stop("THE DATASET DOES NOT HAVE A <Label> or <Comparisons> column")
    }
  }
  
  # Selecting unique genes in each comparison
  pretmp <- dataset[c('Gene', 'Comparisons')]
  pretmp <- unique(pretmp)
  
  tmp = split(pretmp$Gene, pretmp$Comparisons, drop = TRUE)
  
  if (specie == "human") {
    enrichgenes <-
      artms_enrichProfiler(
        tmp,
        categorySource = c(
          'GO:BP',
          'GO:MF',
          'GO:CC',
          'KEGG',
          'REAC',
          'CORUM',
          'HPA',
          'OMIM'
        ),
        specie = 'hsapiens',
        background
      ) # 'HP'
  } else if (specie == "mouse") {
    enrichgenes <-
      artms_enrichProfiler(
        tmp,
        categorySource = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'CORUM'),
        specie = 'mmusculus',
        background
      )
  } else{
    stop("\nSORRY, this specie (",
         specie,
         ") is not supported in the enrichment!!\n")
  }
  
  if (dim(enrichgenes)[1] == 0) {
    cat("--- NO SIGNIFICANT RESULTS from the ENRICHMENT ANALYSIS\n")
  } else{
    if (is.null(heatmaps)) {
      enrichgenes2plot <- enrichgenes[which(enrichgenes$term.size < 500), ]
      for (i in unique(enrichgenes2plot$domain)) {
        cat("\t--- Plotting '", i, "' annotations... ")
        tmp <-
          enrichgenes2plot[which(enrichgenes2plot$domain == i), ]
        outfile <- gsub(".txt", paste0("_", i, ".txt"), output_name)
        .artms_EnrichmentPlotHeatmaps(tmp, outfile)
        cat("out!\n")
      }
    }
  }
  return(enrichgenes)
}


# ------------------------------------------------------------------------------
# @title Plot Corum Enrichment results
#
# @description Heatmap of significant enrichment
# @param df (data.frame) output from `artms_enrichForComplexes`
# @param outfile (char) output file name (must have the extenstion `.pdf`)
# @param theTitle (char) Plot's title
# @return (pdf) heatmap of the significantly enriched protein complexes
# @keywords internal, plot, heatmap, enrichment
.artms_plotCorumEnrichment <- function(df, outfile, theTitle) {
  checkPoint <- length(unique(df$Comparisons))
  if (checkPoint >= 1) {
    # Dealing with the smallest pvalues (pvalue = 0)
    dftemp <- df[-which(df$pvalue == 0), ]
    if (dim(dftemp)[1] > 0) {
      theMinimal <- min(dftemp$pvalue) / 10
      # Make the value replacement
      df$p_value <- df$pvalue
      df$p_value[df$p_value == 0] <- theMinimal
    } else{
      df$p_value <- df$pvalue
    }
    
    df$p_value <- -log10(df$p_value)
    
    toplot <-
      data.table::dcast(
        data = df,
        ComplexName ~ Comparisons,
        value.var = "p_value",
        fun.aggregate = sum,
        fill = 0
      )
    rownames(toplot) <- toplot$ComplexName
    toplotmatrix <- subset(toplot, select = -c(ComplexName))
    x <- data.matrix(toplotmatrix)
    # HEATMAP
    palette.breaks <- seq(1, 3, 0.1)
    color.palette  <-
      colorRampPalette(c("white", "steelblue"))(length(palette.breaks))
    pheatmap::pheatmap(
      x,
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
      color = color.palette
    )
  } else{
    cat("---(-) Not enough enriched comparisons to plot the heatmap\n")
  }
}

# ------------------------------------------------------------------------------
#' @title Enrichment analysis using GprofileR
#'
#' @description This function simplifies the enrichment analysis performed by
#' the excellent tool GprofileR.
#' @param x (list, data.frame) List of protein ids. It can be anything:
#' either a list of ids, or you could also send a data.frame and it will find
#' the columns with the IDs. Is not cool? Multiple list can be also sent
#' simultaneously, as for example running:
#' `tmp <- split(enrichment$Gene, enrichment$cl_number, drop= TRUE)`
#' @param categorySource (vector) Resources providing the terms on which
#' the enrichment will be performed. The supported resources by gprofiler are:
#' - GO (GO:BP, GO:MF, GO:CC): Gene Ontology (see more below)
#' - KEGG: Biological pathways
#' - REAC: Biological pathways (Reactome)
#' - TF: Regulatory motifs in DNA (TRANSFAC TFBS)
#' - MI: Regulatory motifs in DNA (miRBase microRNAs)
#' - CORUM: protein complexes database
#' - HP: Human Phenotype Ontology
#' - HPA: Protein databases (Human Protein Atlas)
#' - OMIM: Online Mendelian Inheritance in Man annotations:
#' - BIOGRID: BioGRID protein-protein interactions
#' The type of annotations for Gene Ontology:
#' - Inferred from experiment (IDA, IPI, IMP, IGI, IEP)
#' - Direct assay (IDA) / Mutant phenotype (IMP]
#' - Genetic interaction (IGI) / Physical interaction (IPI)
#' - Traceable author (TAS) / Non-traceable author (NAS) /
#' Inferred by curator (IC)
#' - Expression pattern (IEP) / Sequence or structural similarity (ISS)
#' / Genomic context (IGC)
#' - Biological aspect of ancestor (IBA) / Rapid divergence (IRD)
#' - Reviewed computational analysis (RCA) / Electronic annotation (IEA)
#' - No biological data (ND) / Not annotated or not in background (NA)
#' @param specie (char) Specie code: Organism names are constructed by
#' concatenating the first letter of the name and the family name.
#' Example: human - ’hsapiens’, mouse - ’mmusculus’. Check gProfileR to find out
#' more about supported species.
#' @param background (vector) gene list to use as background for the enrichment
#' analysis. Default: `NA`
#' @details This function uses the following `gprofiler` arguments as default:
#' - ordered_query = FALSE
#' - significant = TRUE
#' - exclude_iea = TRUE
#' - underrep = FALSE
#' - evcodes = FALSE
#' - region_query = FALSE
#' - max_p_value = 0.05
#' - min_set_size = 0
#' - max_set_size = 0
#' - min_isect_size = 0
#' - correction_method = "analytical" #Options: "gSCS", "fdr", "bonferroni"
#' - hier_filtering = "none"
#' - domain_size = "known" # annotated or known
#' - numeric_ns = ""
#' - png_fn = NULL
#' - include_graph = TRUE
#' @return The enrichment results as provided by gprofiler
#' @keywords enrichment
#' @examples
#' # annotate the MSstats results to get the Gene name
#' data_annotated <- artms_annotationUniprot(
#'                                      data = artms_data_ph_msstats_results,
#'                                      columnid = "Protein",
#'                                      sps = "human")
#'
#' # Filter the list of genes with a log2fc > 2
#' filtered_data <- 
#' unique(data_annotated$Gene[which(data_annotated$log2FC > 2)])
#'
#' # And perform enrichment analysis
#' data_annotated_enrich <- artms_enrichProfiler(
#'                                    x = filtered_data,
#'                                    categorySource = c('KEGG'),
#'                                    specie = "hsapiens",
#'                                    background = unique(data_annotated$Gene))
#' @export
artms_enrichProfiler <-
  function(x,
           categorySource = c('GO'),
           specie,
           background = NA) {
    gProfileR::set_base_url("http://biit.cs.ut.ee/gprofiler")
    cat("---+ Enrichment analysis using gProfiler...")
    enrichData <- gprofiler(
      x,
      organism = specie,
      # "scerevisiae","hsapiens", "mmusculus"
      ordered_query = FALSE,
      significant = TRUE,
      exclude_iea = TRUE,
      # do you want to exclude electronic annotations (IEA)?
      underrep = FALSE,
      evcodes = FALSE,
      region_query = FALSE,
      max_p_value = 0.05,
      min_set_size = 0,
      max_set_size = 0,
      min_isect_size = 0,
      correction_method = "analytical",
      #Options: "gSCS", "fdr", "bonferroni"
      hier_filtering = "none",
      domain_size = "known",
      # annotated or known
      custom_bg = background,
      numeric_ns = "",
      png_fn = NULL,
      include_graph = TRUE,
      src_filter = categorySource
    )
    cat("done!\n")
    return(enrichData)
  }

# Little function to
# ------------------------------------------------------------------------------
# @title Simplify the gProfiler output
#
# @description Simplify the output from `artms_enrichProfiler` resulted from
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

# ------------------------------------------------------------------------------
# @title Plot and save heatmaps of the significant enrichment results
#
# @description plot and save heatmaps of the significant enrichment results
# @param dat (data.frame) output from gprofiler
# @param out_file (char) output file name (must have `.txt` extension)
# @return (pdf) A heatmap of the most significant enrichments
# @keywords internal, plot, heatmap, enrichments
.artms_EnrichmentPlotHeatmaps <- function(dat, out_file) {
  # formatting data to heatmap compatible format
  x <-
    dcast(dat,
          term.name ~ query.number,
          value.var = 'p.value',
          max,
          fill = 1)
  # Let's stop this thing if there is not enough terms (we need at least 2)
  if (dim(x)[1] < 2) {
    return(cat(" NOT ENOUGH TERMS FOR THIS DOMAIN "))
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
        cat(" [SORRY!! We currently don't support heatmaps of a single set] ")
        #pheatmap(term_groups_selected_w_display, 
        #cluster_cols= FALSE,cluster_rows= FALSE, 
        #cellheight=10, cellwidth=10, scale="none", 
        #filename=gsub('.txt','_heatmap.pdf',out_file), 
        #fontsize=6, fontsize_row=8, fontsize_col=8,
        # border_color=NA, color = colors, 
        # fontfamily="Helvetica")
      }
    } else{
      cat(" [SORRY, NOT ENOUGH SIGNIFICANT TERMS FOR THIS DOMAIN] ")
    }
  }
  
}

# ------------------------------------------------------------------------------
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
.artms_foldComplexEnrichment <- function(mylist, corum, background) {
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
