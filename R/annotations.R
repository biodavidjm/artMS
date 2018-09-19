# ------------------------------------------------------------------------------
#' @title Annotate table with Gene Symbol and Name based on Uniprot ID(s)
#' 
#' @description Annotate gene name and symbol based on uniprot ids. It will take
#' the column from your data.frame specified by the `columnid` argument, 
#' search for the gene symbol, name, and entrez based on the specie (`sps`
#' argument) and merge the information back to the input data.frame
#' @param data (data.frame) to be annotated (or file path and name)
#' @param columnid (char) The column with the uniprotkb ids
#' @param sps (char) The specie name. Check `?artms_mapUniprot2entrezGeneName` 
#' to find out more about supported species.
#' @return (data.frame) with two new columns: `Gene` and `Protein.name`
#' @keywords annotation, uniprot
#' @examples{
#' # This example adds annotations to the evidence file available in 
#' # artMS, based on the column "Proteins".
#' 
#' evidence_anno <- artms_annotationUniprot(data = artms_data_ph_evidence,
#'                                          columnid = "Proteins", 
#'                                          sps = "human")
#' }
#' @export
artms_annotationUniprot <- function(data, columnid, sps) {
  
  data <- .artms_checkIfFile(data)
  
  theUniprots <- unique(data[[columnid]])
  preload <- artms_mapUniprot2entrezGeneName(uniprotkb = theUniprots, specie = sps)
  
  dc_merged <- merge(data, preload, by.x = columnid, by.y = 'UNIPROT', all.x = T)
  # Move Gene name to the left:
  gene_first <- preload[,c('SYMBOL', 'UNIPROT', 'GENENAME')]
  dc_merged <- subset(dc_merged, select=-c(SYMBOL, GENENAME))
  send_back <- merge(gene_first, dc_merged, by.x = 'UNIPROT', by.y = columnid, all.y = T)
  names(send_back)[grep('^UNIPROT$', names(send_back) )] <- 'Protein'
  names(send_back)[grep('^SYMBOL$', names(send_back) )] <- 'Gene'
  names(send_back)[grep('^GENENAME$', names(send_back) )] <- 'Protein.names'
  # Some uniprot entries might not have yet a gene name, which will be an 
  # empty value. Replace with Entry name:
  send_back$Gene[which(send_back$Gene == "")] <- NA
  send_back$Gene[is.na(send_back$Gene)] <- send_back$Protein[is.na(send_back$Gene)]
  return(send_back)
}

# ------------------------------------------------------------------------------
#' @title Map UniprotKB to Entrez and Gene Name
#' 
#' @description Map UniprotKB to Entrez and Gene Name
#' @param uniprotkb (vector) Vector of UniprotKB IDs
#' @param specie (char) The specie name. Currently supporting:
#' - ANOPHELES
#' - ARABIDOPSIS
#' - BOVINE
#' - WORM
#' - CANINE
#' - FLY
#' - ZEBRAFISH
#' - ECOLI_STRAIN_K12
#' - ECOLI_STRAIN_SAKAI
#' - CHICKEN
#' - HUMAN
#' - MOUSE
#' - RHESUS
#' - MALARIA
#' - CHIMP
#' - RAT
#' - PIG
#' - XENOPUS
#' @return (data.frame) with EntrezID and GENENAMES mapped on UniprotKB ids
#' @keywords annotation, ids
#' artms_mapUniprot2entrezGeneName()
#' @export
artms_mapUniprot2entrezGeneName <- function(uniprotkb, specie){
  
  specie <- toupper(specie)
  
  if(specie == "ANOPHELES"){
    mappings <- AnnotationDbi::select(org.Ag.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "ARABIDOPSIS"){
    mappings <- AnnotationDbi::select(org.At.tair.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "BOVINE"){
    mappings <- AnnotationDbi::select(org.Bt.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "WORM"){
    mappings <- AnnotationDbi::select(org.Ce.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "CANINE"){
    mappings <- AnnotationDbi::select(org.Cf.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "FLY"){
    mappings <- AnnotationDbi::select(org.Dm.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "ZEBRAFISH"){
    mappings <- AnnotationDbi::select(org.Dr.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "ECOLI_STRAIN_K12"){
    mappings <- AnnotationDbi::select(org.EcK12.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "ECOLI_STRAIN_SAKAI"){
    mappings <- AnnotationDbi::select(org.EcSakai.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "CHICKEN"){
    mappings <- AnnotationDbi::select(org.Gg.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "HUMAN"){
    mappings <- AnnotationDbi::select(org.Hs.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "MOUSE"){
    mappings <- AnnotationDbi::select(org.Mm.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "RHESUS"){
    mappings <- AnnotationDbi::select(org.Mmu.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "MALARIA"){
    mappings <- AnnotationDbi::select(org.Pf.plasmo.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "CHIMP"){
    mappings <- AnnotationDbi::select(org.Pt.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "RAT"){
    mappings <- AnnotationDbi::select(org.Rn.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "YEAST"){
    mappings <- AnnotationDbi::select(org.Sc.sgd.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "PIG"){
    mappings <- AnnotationDbi::select(org.Ss.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else if(specie == "XENOPUS"){
    mappings <- AnnotationDbi::select(org.Xl.eg.db, uniprotkb, c("UNIPROT", "SYMBOL","GENENAME", "ENTREZID"), keytype = "UNIPROT")
  }else{
    cat("ERROR: Specie not supported.")
    stop("PLEASE, CHECK HELP TO FIND OUT MORE ABOUT SUPPORTED SPECIES")
  }
  mappings <- unique(mappings)
  # It migth come with 1 uniprot to many gene names: take the first one, 
  # which should be the main gene name
  mappings <- mappings[!duplicated(mappings$UNIPROT),]
  return(mappings)
}




