# ------------------------------------------------------------------------------
#' @title Annotate table with Gene Symbol and Name based on Uniprot ID(s)
#'
#' @description Annotate gene name and symbol based on uniprot ids. It will 
#' take the column from your data.frame specified by the `columnid` argument,
#' search for the gene symbol, name, and entrez based on the species (`species`
#' argument) and merge the information back to the input data.frame
#' @param x (data.frame) to be annotated (or file path and name)
#' @param columnid (char) The column with the uniprotkb ids
#' @param species (char) The species name. Check `?artmsMapUniprot2Entrez`
#' to find out more about supported species.
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame) with two new columns: `Gene` and `Protein.name`
#' @keywords annotation, uniprot
#' @examples
#' # This example adds annotations to the evidence file available in
#' # artMS, based on the column 'Proteins'.
#'
#' evidence_anno <- artmsAnnotationUniprot(x = artms_data_ph_evidence,
#'                                          columnid = 'Proteins',
#'                                          species = 'human')
#' @export
artmsAnnotationUniprot <- function(x, 
                                   columnid, 
                                   species,
                                   verbose = TRUE) {
  
  if(any(missing(x) | 
         missing(columnid) |
         missing(species)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  if(!is.character(columnid)) stop("Argument 'columnid' must be a character")
  if(!is.character(species)) stop("Argument <species> must be a character")
  
  x <- .artms_checkIfFile(x)
  
  # Deal with isoforms UNIPROT-1 etc
  if(!columnid %in% colnames(x)){
    stop("The column id <",columnid,"> was not found in this dataset")
  }
  x$EntryRoot <- x[[columnid]]
  
  # Check first if it is uniprot:
  isUniprot <- TRUE
  if( length(x$EntryRoot[grep("\\w{2}_\\d{1,}\\.\\d{1,}", x$EntryRoot)]) > 100 ){
    message("\nWarning!
    It looks like this dataset is using RefSeq as protein ID.
    RefSeq is not (yet) supported for gene annotation.\n")
    isUniprot <- FALSE
  }
  
  x$EntryRoot <- gsub("(\\w+)(-)(\\d+)", "\\1", x$EntryRoot)
  
  # Deal with protein groups and PTMs: select the first ID 
  # 1. strsplit based on ;
  x$EntryRoot <- unlist(lapply(x$EntryRoot, function(x) unlist(strsplit(x, "\\;"))[1]))
  # 2. now based on PTM: either STYK
  x$EntryRoot <- unlist(lapply(x$EntryRoot, function(x) unlist(strsplit(x, "_[S|T|Y|K]"))[1]))
  
  theUniprots <- as.character(unique(x$EntryRoot))
  
  if( isFALSE(artmsIsSpeciesSupported(species = species)) | isFALSE(isUniprot) ){
    if(isFALSE(artmsIsSpeciesSupported(species = species))) 
      message("---(-) Species ", species," not supported: no info about proteins will be provided")
    keepSearchName <- paste0(columnid, "Copy")
    y <- artmsChangeColumnName(x, columnid, keepSearchName)
    y <- artmsChangeColumnName(y, "EntryRoot", "Protein")
    y$Gene <- y$Protein
    y$ProteinName <- NA
    y$EntrezID <- NA
    send_back <- y[c("Protein", "Gene", "ProteinName", "EntrezID")]
    send_back <- unique(send_back)
    y <- subset(y, select = -c(Gene, ProteinName, EntrezID))
    send_back <- merge(send_back, y, by = "Protein")
    return(send_back)
  }else{
    preload <- artmsMapUniprot2Entrez(uniprotkb = theUniprots, 
                                      species = species)
    
    dc_merged <-merge(x, 
                      preload,
                      by.x = "EntryRoot",
                      by.y = "UNIPROT",
                      all.x = TRUE)
    
    # Move Gene name to the left:
    gene_first <- preload[, c("SYMBOL", "UNIPROT", "GENENAME", "ENTREZID")]
    dc_merged <- subset(dc_merged, select = -c(SYMBOL, GENENAME, ENTREZID))
    send_back <-
      merge(
        gene_first,
        dc_merged,
        by.x = "UNIPROT",
        by.y = "EntryRoot",
        all.y = TRUE
      )
    keepSearchName <- paste0(columnid,"Ref")
    send_back <- artmsChangeColumnName(send_back, columnid, keepSearchName)
    names(send_back)[grep("^UNIPROT$", names(send_back))] <- "Protein"
    names(send_back)[grep("^SYMBOL$", names(send_back))] <- "Gene"
    names(send_back)[grep("^GENENAME$", names(send_back))] <- "ProteinName"
    names(send_back)[grep("^ENTREZID$", names(send_back))] <- "EntrezID"
    # Some uniprot entries might not have yet a gene name, 
    # which will be an empty value. Replace with Entry name:
    send_back$Gene[which(send_back$Gene == "")] <- NA
    send_back$Gene[is.na(send_back$Gene)] <- send_back$Protein[is.na(send_back$Gene)]
    return(send_back)
  }
}


# ------------------------------------------------------------------------------
# @title Select the entry 
#
# @description Selet the Uniprot Entry ID from a full Uniprot id. For example:
#
# From `sp|P55011|S12A2_HUMAN` will select `P55011`
#
# @param x (vector) of protein ids. If the id is a uniprot full entry, it will
# extract the ENTRY id
# @return (vector) with only the Entry ID (if found)
# @keywords internal, Protein ID, Uniprot, Entry
.artms_selectEntryFromFullUniprot <- function(x){
  isolateIt <-  unlist(strsplit(x, ";"))
  replaceIt <- as.character(sapply(isolateIt, function(y) gsub("(^sp|tr)(\\|)(.*)(\\|.*)", "\\3", y)))
  if(length(replaceIt) > 1){
    replaceIt <- paste(replaceIt, collapse = ";")
  }
  return(replaceIt)
}


# ------------------------------------------------------------------------------
#' @title Leave only the Entry ID from a typical full Uniprot IDs in a 
#' given column 
#'
#' @description Downloading a Reference Uniprot fasta database includes several 
#' Uniprot IDs for every protein. If the regular expression available in 
#' Maxquant is not activated, the full id will be used in the Proteins,
#' Lead Protein, and Leading Razor Protein columns. This script leaves only the
#' Entry ID.
#' 
#' For example, values in a Protein column like this:
#' 
#' `sp|P12345|Entry_name;sp|P54321|Entry_name2` 
#' 
#'  will be replace by
#'  
#' `P12345;P54321``
#' 
#' @param x (data.frame) that contains the `columnid`
#' @param columnid (char) Column name with the full uniprot ids
#' @return (data.frame) with only Entry IDs.
#' 
#' @keywords annotation, ids
#' @examples
#' # Example of data frame with full uniprot ids and sequences
#' p <- c("sp|A6NIE6|RN3P2_HUMAN;sp|Q9NYV6|RRN3_HUMAN", 
#'        "sp|A7E2V4|ZSWM8_HUMAN", 
#'        "sp|A5A6H4|ROA1_PANTR;sp|P09651|ROA1_HUMAN;sp|Q32P51|RA1L2_HUMAN", 
#'        "sp|A0FGR8|ESYT2_HUMAN")
#' s <- c("ALENDFFNSPPRK", "GWGSPGRPK", "SSGPYGGGGQYFAK", "VLVALASEELAK")
#' evidence <- data.frame(Proteins = p, Sequences = s, stringsAsFactors = FALSE)
#' 
#' # Replace the Proteins column with only Entry ids
#' evidence <- artmsLeaveOnlyUniprotEntryID(x = evidence, columnid = "Proteins")
#' @export
artmsLeaveOnlyUniprotEntryID <- function(x, columnid){
  
  if(!(is.data.frame(x) | is.data.table(x))){
    stop("<x> must be a data.frame or data.table")
  }
  
  if(!is.vector(columnid)){
    stop("<columnid> must be a valid column name")
  }
  
  if(!(columnid %in% colnames(x))){
    stop("The <columnid> is not a column of <x>. Check the name again")
  }
  
  x[[columnid]] <- unlist(lapply(x[[columnid]], 
                                 function(z) sapply(z, .artms_selectEntryFromFullUniprot)))
  
  return(x)
}


# ------------------------------------------------------------------------------
#' @title Map GENE SYMBOL, NAME, AND ENTREZID to a vector of Uniprot IDS
#'
#' @description Map GENE SYMBOL, NAME, AND ENTREZID to a vector of Uniprot IDS
#' @param uniprotkb (vector) Vector of UniprotKB IDs
#' @param species (char) The species name. Species currently supported 
#' as part of artMS: check `?artmsIsSpeciesSupported()` to find out the
#' list of supported species`
#' @return (data.frame) with ENTREZID and GENENAMES mapped on UniprotKB ids
#' @keywords annotation, ids
#' @examples
#' # Load an example with human proteins
#' exampleID <- c("Q6P996", "B1N8M6")
#' artmsMapUniprot2Entrez(uniprotkb = exampleID, 
#'                        species = "HUMAN")
#'                        
#' @export
artmsMapUniprot2Entrez <- function(uniprotkb, 
                                   species) {
  
  if(any(missing(uniprotkb) | 
         missing(species)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  if(!is.vector(uniprotkb)) stop("Argument <uniprotkb> is not a vector")
  if(!is.character(species)) stop("Argument <species> is not a vector")

  if(isFALSE(artmsIsSpeciesSupported(species = species))){
    stop("Species ", species," not supported")
  }else{
    thePack <- artmsIsSpeciesSupported(species = species)
  }

  suppressMessages(
    mappings <-
      AnnotationDbi::select(
        eval(as.symbol(thePack)),
        uniprotkb,
        c("UNIPROT", "SYMBOL",
          "GENENAME", "ENTREZID"),
        keytype = "UNIPROT"
      )
  )

  mappings <- unique(mappings)
  # It migth come with 1 uniprot to many gene names: 
  # take the first one, which should be the main gene name
  mappings <- mappings[!duplicated(mappings$UNIPROT),]
  return(mappings)
}


# ------------------------------------------------------------------------------
#' @title Check if a species is supported and available
#'
#' @description Given a species name, it checkes whether is supported, and 
#' if supported, check whether the annotation package is installed.
#' @param species (char) The species name. Species currently supported 
#' as part of artMS:
#' - HUMAN
#' - MOUSE
#' 
#' And the following species can be used as well, but the user needs to 
#' install the corresponding org.db package:
#' - ANOPHELES (`install.packages(org.Ag.eg.db)`)
#' - BOVINE (`install.packages(org.Bt.eg.db)`)
#' - WORM (`install.packages(org.Ce.eg.db)`)
#' - CANINE (`install.packages(org.Cf.eg.db)`)
#' - FLY (`install.packages(org.Dm.eg.db)`)
#' - ZEBRAFISH (`install.packages(org.Dr.eg.db)`)
#' - CHICKEN (`install.packages(org.Gg.eg.db)`)
#' - RHESUS (`install.packages(org.Mmu.eg.db)`)
#' - CHIMP (`install.packages(org.Pt.eg.db)`)
#' - RAT (`install.packages(org.Rn.eg.db)`)
#' - YEAST (`install.packages(org.Sc.sgd.db)`)
#' - PIG (`install.packages(org.Ss.eg.db)`)
#' - XENOPUS (`install.packages(org.Xl.eg.db)`)
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (string) Name of the package for the given species
#' @keywords annotation, species
#' @examples
#' # Should return TRUE
#' artmsIsSpeciesSupported(species = "HUMAN")
#' artmsIsSpeciesSupported(species = "CHIMP")
#' @export
artmsIsSpeciesSupported <- function(species, 
                                    verbose = TRUE) {
  
  if(any(missing(species)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  if(!is.character(species)) stop("Argument <species> is not a vector")
  
  species <- toupper(species)
  
  LongName <- c(
    "ANOPHELES", 
    "BOVINE",
    "WORM",
    "CANINE",
    "FLY",
    "ZEBRAFISH",
    "CHICKEN",
    "HUMAN",
    "MOUSE",
    "RHESUS",
    "CHIMP",
    "RAT",
    "YEAST",
    "PIG",
    "XENOPUS")
  PackageName <- c(
    "org.Ag.eg.db", 
    "org.Bt.eg.db",
    "org.Ce.eg.db",
    "org.Cf.eg.db",
    "org.Dm.eg.db",
    "org.Dr.eg.db",
    "org.Gg.eg.db",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "org.Mmu.eg.db",
    "org.Pt.eg.db",
    "org.Rn.eg.db",
    "org.Sc.sgd.db",
    "org.Ss.eg.db",
    "org.Xl.eg.db")
  
  OrgDB <- data.frame("LongName" = LongName, 
                      "PackageName" = PackageName, 
                      stringsAsFactors = FALSE)
  
  if(species %in% OrgDB$LongName){
    thePack <- OrgDB$PackageName[which(OrgDB$LongName == species)]
    
    if( !(thePack %in% rownames(installed.packages())) ){
      if(verbose) message("---(-) The package <",thePack,"> is not installed in your system.
           Just run: install.packages('", thePack,"') and try again")
      return(FALSE)
    } else{
      return(thePack)
    }
  }else{
    return(FALSE)
  }
}




