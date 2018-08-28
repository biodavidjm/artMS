# ------------------------------------------------------------------------------
#' @title Map UniprotKB to Entrez and Gene Name
#' 
#' @description Map UniprotKB to Entrez and Gene Name
#' @param theUniprots List of UniprotKB IDs
#' @param specie The specie name. Currently supporting:
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
#' @return A data.frame with EntrezID and GENENAMES mapped on UniprotKB ids
#' @keywords annotation, ids
#' artms_mapUniprot2entrezGeneName()
#' @export
artms_mapUniprot2entrezGeneName <- function(theUniprots, specie){
  
  specie <- toupper(specie)
  
  if(specie == "ANOPHELES"){
    mappings <- AnnotationDbi::select(org.Ag.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "ARABIDOPSIS"){
    mappings <- AnnotationDbi::select(org.At.tair.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "BOVINE"){
    mappings <- AnnotationDbi::select(org.Bt.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "WORM"){
    mappings <- AnnotationDbi::select(org.Ce.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "CANINE"){
    mappings <- AnnotationDbi::select(org.Cf.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "FLY"){
    mappings <- AnnotationDbi::select(org.Dm.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "ZEBRAFISH"){
    mappings <- AnnotationDbi::select(org.Dr.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "ECOLI_STRAIN_K12"){
    mappings <- AnnotationDbi::select(org.EcK12.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "ECOLI_STRAIN_SAKAI"){
    mappings <- AnnotationDbi::select(org.EcSakai.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "CHICKEN"){
    mappings <- AnnotationDbi::select(org.Gg.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "HUMAN"){
    mappings <- AnnotationDbi::select(org.Hs.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "MOUSE"){
    mappings <- AnnotationDbi::select(org.Mm.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "RHESUS"){
    mappings <- AnnotationDbi::select(org.Mmu.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "MALARIA"){
    mappings <- AnnotationDbi::select(org.Pf.plasmo.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "CHIMP"){
    mappings <- AnnotationDbi::select(org.Pt.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "RAT"){
    mappings <- AnnotationDbi::select(org.Rn.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "YEAST"){
    mappings <- AnnotationDbi::select(org.Sc.sgd.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "PIG"){
    mappings <- AnnotationDbi::select(org.Ss.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else if(specie == "XENOPUS"){
    mappings <- AnnotationDbi::select(org.Xl.eg.db, uniprots, c("UNIPROT", "SYMBOL","GENENAME"), keytype = "UNIPROT")
  }else{
    cat("ERROR: Specie not supported.")
    stop("PLEASE, CHECK HELP TO FIND OUT MORE ABOUT SUPPORTED SPECIES")
  }
  return(mappings)
}

