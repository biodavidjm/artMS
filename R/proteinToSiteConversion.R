# ------------------------------------------------------------------------------
#' @title Converts the `Proteins` column of the evidence file to site-specific
#' `Uniprot_PTM` notation
#'
#' @description It enables the site-specific quantification of PTMs by
#' converting the `Proteins` column of the evidence file to an `Uniprot_PTM`
#' or `RefSeq_PTM` notation. 
#' In this way, each of the modified peptides can be quantified
#' independently across conditions. 
#' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#' WARNING: we have detected a version of MaxQuant (1.6.3.4) outputs a`
#' "Modified sequence" column of the evidence file that has two important 
#' changes for the annotation of phosphorylation:
#' - Uses `p` instead of `(ph)`
#' - The modified residue (i.e. STY) is the residue on the left of the `p`, 
#' instead of the residue to the right of `(ph)` as usual.
#' These changes significantly affect the algorithm used by this artMS function
#' to extract the PTM site information. We are working on debugging the issue
#' and use an alternative function independent of the "Mofified sequence" 
#' column. If your MaxQuant version is affected, then the column `mod_sites` of
#' the output `-mapping.txt` will have many NA values for peptides with `p`
#' modifications (and no NA values should be present)
#' @param evidence_file (char) The evidence file name and location
#' @param ref_proteome_file (char) The reference proteome used as database
#' to search the `evidence.txt` file with MaxQuant. It will be used to map the
#' modified peptide to the protein sequence and find the site location.
#' Therefore, it does not use the MaxQuant's `Phospho (STY)Sites.txt`
#' @param output_file (char) Output file name 
#' (`-sites-evidence.txt` recommended)
#' @param mod_type (char) The posttranslational modification. Options:
#' - `UB`: Protein Ubiquitination
#' - `PH`: Protein Phosphorylation
#' - `AC`: Protein Acetylation
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (file) Return a new evidence file with the `Proteins` column
#' modified by adding the sequence site location(s) + postranslational
#' modification(s) to the uniprot entry / refseq id.
#' 
#' Examples: `A34890_ph3`; `Q64890_ph24_ph456`; `Q64890_ub34_ub129_ub234`;
#' `Q64890_ac35`.
#' @keywords evidence, convert, ptm, ph, ub, ac
#' @examples
#' # Testing warning if files are not submitted
#' artmsProtein2SiteConversion(evidence_file = NULL, ref_proteome_file = NULL, 
#' output_file = NULL)
#' @export
artmsProtein2SiteConversion <- function (evidence_file,
                                           ref_proteome_file,
                                           output_file,
                                           mod_type,
                                           verbose = TRUE) {
  if(is.null(evidence_file) & 
     is.null(ref_proteome_file) & 
     is.null(output_file)){
    return("Files must not be NULL")
  }
  
  if(any(missing(evidence_file) | 
         missing(ref_proteome_file) |
         missing(output_file) | 
         missing(mod_type)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  mod_type <- toupper(mod_type)
  if(!mod_type %in% c("PH", "UB", "AC"))
    stop("the mod_type ", mod_type, " is not supported")
  
  if(verbose) message(">> CONVERTING EVIDENCE TO PTM SITE-SPECIFIC ")
  
  # When the evidence file is not a data.frame, check that it does not have
  # the same as the output file
  if(is.vector(evidence_file)){
    if(evidence_file == output_file) 
      stop("<output_file> cannot be the same as <evidence_file>")
  }
  
  if(verbose)
    message(">> PROCESSING THE EVIDENCE FILE FOR A SITE SPECIFIC ANALYSIS ")
  
  if (mod_type == 'UB') {
    if(verbose) message('--- SELECTING << UB >> MODIFIED PEPTIDES ')
    maxq_mod_residue = 'K\\(gl\\)'
    mod_residue = 'K'
  } else if (mod_type == 'PH') {
    if(verbose) message('--- SELECTING << PH >> MODIFIED PEPTIDES ')
    maxq_mod_residue = '(S|T|Y)\\(ph\\)'
    mod_residue = 'S|T|Y'
  } else if (mod_type == 'AC') {
    if(verbose) message('--- SELECTING << AC >> MODIFIED PEPTIDES ')
    maxq_mod_residue = 'K\\(ac\\)'
    mod_residue = 'K'
  } else{
    stop(
      mod_type, " is not supported. 
      Check help to get the list of supported PTMs"
    )
  }
  
  ## map mod sites in data to index
  if(verbose) message("--- OPENING EVIDENCE FILE ")
  ## read in maxq. data
  maxq_data <- .artms_checkIfFile(evidence_file)
  maxq_data <- as.data.table(maxq_data)
  
  # remove contaminants, keep unique sequences, fix names
  maxq_data <-
    maxq_data[grep("CON__|REV__", maxq_data$Proteins, invert = TRUE), ]
  
  # In case the evidence file comes as data.frame
  if("Modified sequence" %in% colnames(maxq_data)){
    if(verbose) message("\t--(+) 'Modified sequence' column available") 
  }else if("Modified.sequence" %in% colnames(maxq_data)){
    setnames(maxq_data, 'Modified.sequence', 'Modified sequence')  
  }else{
    stop("Problem with evidence file: no <Modified sequence> column found")
  }
  
  unique_peptides_in_data <-
    unique(maxq_data[, c('Proteins', 'Modified sequence'), with = FALSE])
  setnames(unique_peptides_in_data, 'Modified sequence', 'sequence') 
  
  
  if(verbose) message("--- READING REFERENCE PROTEOME ")
  ## read in reference proteome
  ref_proteome <- read.fasta(
    file = ref_proteome_file,
    seqtype = "AA",
    as.string = TRUE,
    set.attributes = TRUE,
    legacy.mode = TRUE,
    seqonly = FALSE,
    strip.desc = FALSE
  )
  
  ## make mod-site index
  p_seqs <- c()
  p_names <- c()
  p_annots <- c()
  
  for (e in ref_proteome) {
    p_seqs <- c(p_seqs, e[1])
    p_names <- c(p_names, attr(e, 'name'))
    p_annots <- c(p_annots, attr(e, 'Annot'))
  }
  
  ref_table <-
    data.table(names = p_names,
               annots = p_annots,
               seqs = p_seqs)
  
  # REFSEQ database
  ref_table[, uniprot_ac := gsub('([a-z,0-9,A-Z,\\.,\\_]+)', '\\1', names)]
  
  # UNIPROT
  ref_table[, 
  uniprot_ac := gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)',
                       '\\2', names)]
  
  # get all indicies/locations of the mod_residue S|T|Y
  indices <-
    lapply(ref_table$seqs, function(x)
      as.vector(str_locate_all(x, pattern = mod_residue)[[1]][, 1]))
  # list which residue it actually finds at each location
  ptm_sites <-
    lapply(ref_table$seqs, function(x)
      as.vector(str_match_all(x, pattern = mod_residue)))
  # get the number of residue locations per string
  lengths <- unlist(lapply(indices, FUN = length))
  # repeate the uniprot_ac as many times as a residue (S|T|Y) is found
  keys <- rep(ref_table$uniprot_ac, lengths)
  # combine the list of the exploded proteins (above) and add in which sites
  # are found as well as the location the residue was found in
  protein_indices <-
    data.table(
      uniprot_ac = keys,
      ptm_site = unlist(ptm_sites),
      res_index = unlist(indices)
    )

  
  # EXTRACT PTM POSITIONS FROM MODIFIED PEPTIDES IN EVIDENCE FILE
  mod_sites <- c()
  mod_seqs <- c()
  
  if(verbose) message("--- EXTRACTING PTM POSITIONS FROM THE MODIFIED PEPTIDES 
      (it might take some time) ")
  for (i in seq_len(nrow(unique_peptides_in_data))) {
    entry <- unique_peptides_in_data[i, ]
    peptide_seq <- entry$sequence
    ## cleanup the sequence (removing all modifications) for matching the 
    ## protein sequence
    peptide_seq_clean <- gsub('[a-z,0-9,\\(,\\),_]', '', peptide_seq)
    mod_sites_in_peptide <-
      str_locate_all(string = peptide_seq, pattern = maxq_mod_residue)[[1]][, 1]
    
    if (length(mod_sites_in_peptide) > 0) {
      uniprot_acs <- entry$Proteins
      # separates the ambiguous cases (;) and appends the site info to all 
      # the proteins
      uniprot_acs <-
        str_split(string = uniprot_acs, pattern = ';')[[1]]
      
      for (uac in uniprot_acs) {
        # find the protein's full sequence based on the uniprot_ac
        protein_seq <- ref_table[uniprot_ac == uac, ]$seqs
        if (length(protein_seq) > 0) {
          ## get the position of the peptide in the protein sequence
          peptide_index_in_protein <-
            str_locate(protein_seq, peptide_seq_clean)[[1]][1]
          
          for (m in seq_len(length(mod_sites_in_peptide))) {
            mod_site <- mod_sites_in_peptide[m]
            peptide_seq_before_site <-
              str_sub(peptide_seq, 1, mod_site - 1)
            ## count all AA (not counting all modifications) before the
            ## modification to get the relative position of the modification
            ## in the peptide sequence
            residues_before_site <-
              str_count(string = peptide_seq_before_site, 
                        pattern = 'A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y')
            mod_site_index_in_protein <-
              peptide_index_in_protein + residues_before_site
            protein_mod_sites <- protein_indices[uniprot_ac == uac, ]
            if (!is.na(mod_site_index_in_protein)) {
              mod_res <-
                protein_mod_sites[res_index == mod_site_index_in_protein, 
                                  ptm_site]
              mod_site_id <-
                sprintf(
                  '%s_%s%s',
                  uac,
                  str_sub(
                    protein_seq,
                    mod_site_index_in_protein,
                    mod_site_index_in_protein
                  ),
                  mod_site_index_in_protein
                )
              mod_sites <- c(mod_sites, mod_site_id)
              mod_seqs <- c(mod_seqs, peptide_seq)
              # message("mod_site_id ", mod_site_id, " peptide_seq ",peptide_seq)
              # stopifnot(length(mod_sites) == length(mod_seqs))
            } else{
              message(
                sprintf(
                  'MISMATCH\t%s\n\tPEPTIDE_SEQ\t%s\n\tMOD_SITE\t%s\n\t
                  PEPTIDE_IDX_IN_PROTEIN\t%s\n\tRESIDUES_BEFORE_SITE
                  \t%s\n\tPROTEIN_SEQ\t%s\n',
                  mod_site_id,
                  peptide_seq,
                  mod_site,
                  peptide_index_in_protein,
                  residues_before_site,
                  protein_seq
                )
              )
            }
          }
        }else{
          message("- Protein <", uac, "> not in the sequence database\n")
        }
      }
    }
  }
  
  mod_site_mapping <- data.table(mod_sites, mod_seqs)
  mod_site_mapping_agg <-
    aggregate(
      mod_sites ~ mod_seqs,
      mod_site_mapping,
      FUN = function(x)
        paste(x, collapse = ',')
    )
  
  setnames(maxq_data, 'Modified sequence', 'mod_seqs')
  unmapped_mod_seqs <-
    maxq_data[!(mod_seqs %in% mod_site_mapping_agg$mod_seqs) &
                grepl('(gl)', mod_seqs) & !grepl('REV__|CON__', Proteins), ]
  unmapped_mod_seqs <-
    unique(unmapped_mod_seqs[, c('mod_seqs', 'Proteins'), with = FALSE])
  if (dim(unmapped_mod_seqs)[1] > 0) {
    if(verbose) message('>> UNABLE TO MAP \t')
    if(verbose) print(unmapped_mod_seqs)
  } else{
    if(verbose) message("--- ALL SEQUENCES MAPPED ")
  }
  
  final_data <-
    merge(maxq_data, mod_site_mapping_agg, by = 'mod_seqs')
  setnames(
    final_data,
    c('Proteins', 'mod_sites', 'mod_seqs'),
    c('Proteins_ref', 'Proteins', 'Modified sequence')
  )
  write.table(
    final_data,
    file = output_file,
    eol = '\n',
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  ## write a mapping table
  protein_seq_mapping <-
    unique(maxq_data[, c('Proteins', 'mod_seqs'), with = FALSE])
  setnames(protein_seq_mapping, 'Proteins', 'Protein')
  mapping_table <-
    merge(protein_seq_mapping,
          mod_site_mapping_agg,
          by = 'mod_seqs',
          all = TRUE)
  write.table(
    mapping_table,
    file = gsub('.txt', '-mapping.txt', output_file),
    eol = '\n',
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  if(verbose) 
    message(
      ">> FILES OUT:\n",
      "--- New evidence-site file: ",
      output_file,
      "\n",
      "--- Details of the Mappings: ",
      gsub('.txt', '-mapping.txt', output_file)
    )
  if(verbose) message(">> CONVERSION COMPLETED")
}
