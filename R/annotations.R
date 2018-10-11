# ------------------------------------------------------------------------------
#' @title Annotate table with Gene Symbol and Name based on Uniprot ID(s)
#'
#' @description Annotate gene name and symbol based on uniprot ids. It will 
#' take the column from your data.frame specified by the `columnid` argument,
#' search for the gene symbol, name, and entrez based on the species (`sps`
#' argument) and merge the information back to the input data.frame
#' @param x (data.frame) to be annotated (or file path and name)
#' @param columnid (char) The column with the uniprotkb ids
#' @param sps (char) The species name. Check `?artmsMapUniprot2Entrez`
#' to find out more about supported species.
#' @return (data.frame) with two new columns: `Gene` and `Protein.name`
#' @keywords annotation, uniprot
#' @examples
#' # This example adds annotations to the evidence file available in
#' # artMS, based on the column 'Proteins'.
#'
#' evidence_anno <- artms_annotationUniprot(x = artms_data_ph_evidence,
#'                                          columnid = 'Proteins',
#'                                          sps = 'human')
#' @export
artms_annotationUniprot <- function(x, columnid, sps) {
  
  if(any(missing(x) | 
         missing(columnid) |
         missing(sps)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  if(!is.character(columnid)) stop("Argument 'columnid' must be a character")
  if(!is.character(sps)) stop("Argument <sps> must be a character")
  
  x <- .artms_checkIfFile(x)
  
  theUniprots <- as.character(unique(x[[columnid]]))
  
  preload <- artmsMapUniprot2Entrez(uniprotkb = theUniprots, 
                                             species = sps)
  
  dc_merged <-merge(x, 
                    preload,
                    by.x = columnid,
                    by.y = "UNIPROT",
                    all.x = TRUE)
  
  # Move Gene name to the left:
  gene_first <- preload[, c("SYMBOL", "UNIPROT", "GENENAME")]
  dc_merged <- subset(dc_merged, select = -c(SYMBOL, GENENAME))
  send_back <-
    merge(
      gene_first,
      dc_merged,
      by.x = "UNIPROT",
      by.y = columnid,
      all.y = TRUE
    )
  names(send_back)[grep("^UNIPROT$", names(send_back))] <- "Protein"
  names(send_back)[grep("^SYMBOL$", names(send_back))] <- "Gene"
  names(send_back)[grep("^GENENAME$", names(send_back))] <- "Protein.names"
  # Some uniprot entries might not have yet a gene name, 
  # which will be an empty value. Replace with Entry name:
  send_back$Gene[which(send_back$Gene == "")] <- NA
  send_back$Gene[is.na(send_back$Gene)] <-
    send_back$Protein[is.na(send_back$Gene)]
  return(send_back)
}

# ------------------------------------------------------------------------------
#' @title Map GENE SYMBOL, NAME, AND ENTREZID to a vector of Uniprot IDS
#'
#' @description Map GENE SYMBOL, NAME, AND ENTREZID to a vector of Uniprot IDS
#' @param uniprotkb (vector) Vector of UniprotKB IDs
#' @param species (char) The species name. Species currently supported 
#' as part of artMS:
#' - HUMAN
#' - MOUSE
#' 
#' And the following species can be used as well, but the user needs to 
#' install the corresponding org.db package:
#' - ANOPHELES (`install.packages(org.Ag.eg.db)`)
#' - ARABIDOPSIS (`install.packages(org.At.tair.db)`)
#' - BOVINE (`install.packages(org.Bt.eg.db)`)
#' - WORM (`install.packages(org.Ce.eg.db)`)
#' - CANINE (`install.packages(org.Cf.eg.db)`)
#' - FLY (`install.packages(org.Dm.eg.db)`)
#' - ZEBRAFISH (`install.packages(org.Dr.eg.db)`)
#' - ECOLI_STRAIN_K12 (`install.packages(org.EcK12.eg.db)`)
#' - ECOLI_STRAIN_SAKAI (`install.packages(org.EcSakai.eg.db)`)
#' - CHICKEN (`install.packages(org.Gg.eg.db)`)
#' - RHESUS (`install.packages(org.Mmu.eg.db)`)
#' - MALARIA (`install.packages(org.Pf.plasmo.db)`)
#' - CHIMP (`install.packages(org.Pt.eg.db)`)
#' - RAT (`install.packages(org.Rn.eg.db)`)
#' - YEAST (`install.packages(org.Sc.sgd.db)`)
#' - PIG (`install.packages(org.Ss.eg.db)`)
#' - XENOPUS (`install.packages(org.Xl.eg.db)`)
#' @return (data.frame) with EntrezID and GENENAMES mapped on UniprotKB ids
#' @keywords annotation, ids
#' @examples
#' # Load an example
#' exampleID <- c("Q6P996", "B1N8M6")
#' artmsMapUniprot2Entrez(uniprotkb = exampleID, 
#'                        species = "HUMAN")
#' @export
artmsMapUniprot2Entrez <- function(uniprotkb, 
                                   species) {
  
  if(any(missing(uniprotkb) | 
         missing(species)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  if(!is.vector(uniprotkb)) stop("Argument <uniprotkb> is not a vector")
  if(!is.character(species)) stop("Argument <species> is not a vector")
  
  species <- toupper(species)
  
  LongName <- c(
    "ANOPHELES", 
    "ARABIDOPSIS",
    "BOVINE",
    "WORM",
    "CANINE",
    "FLY",
    "ZEBRAFISH",
    "ECOLI_STRAIN_K12",
    "ECOLI_STRAIN_SAKAI",
    "CHICKEN",
    "HUMAN",
    "MOUSE",
    "RHESUS",
    "MALARIA",
    "CHIMP",
    "RAT",
    "YEAST",
    "PIG",
    "XENOPUS")
  PackageName <- c(
    "org.Ag.eg.db", 
    "org.At.tair.db",
    "org.Bt.eg.db",
    "org.Ce.eg.db",
    "org.Cf.eg.db",
    "org.Dm.eg.db",
    "org.Dr.eg.db",
    "org.EcK12.eg.db",
    "org.EcSakai.eg.db",
    "org.Gg.eg.db",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "org.Mmu.eg.db",
    "org.Pf.plasmo.db",
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
    
    if( !(thePack %in% rownames(installed.packages())) )
      stop("---(-) The package <",thePack,"> is not installed in your system.
           Just run: install.packages('",thePack,"') and try again")
    
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
  }else{
    stop("Specie ", species, " not supported. 
         Please, check help to find out more about supported species")
  }
  
  mappings <- unique(mappings)
  # It migth come with 1 uniprot to many gene names: 
  # take the first one, which should be the main gene name
  mappings <- mappings[!duplicated(mappings$UNIPROT),]
  return(mappings)
}
