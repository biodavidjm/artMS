
#' @import data.table
#' @import seqinr
#' @import stringr
#' @import bit64
#' @import getopt
#' @import reshape2
#' @import biomaRt
#' @import limma
#' @import ggplot2
#' @import bit64



# ~~~~~~~~~~~~~~~~~~~ BEGIN FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~


#' @title Convert Silac to Long Format
#' @description This function converts a silac labeled MaxQuant file and splits the H/L cases apart into separate rows.
#' @param filename The filepath to the MaxQuant silac data (txt tab delimited file).
#' @param output The filepath to the intended output file (txt tab delimited file).
#' @keywords silac
#' MQutil.SILACToLong()
MQutil.SILACToLong = function(filename, output) {
    library(data.table)
    library(reshape2)
    file = Sys.glob(filename)
    cat(sprintf("\tPROCESSING:\n\t%s\n", paste(file, collapse = "\n\t")))
    tmp = fread(file, integer64 = "double")
    
    # reshape the data and split the heavy and light data
    tmp_long = reshape2::melt(tmp, measure.vars = c("Intensity L", "Intensity H"))
    tmp_long[, `:=`(Intensity, NULL)]
    setnames(tmp_long, "value", "Intensity")
    setnames(tmp_long, "variable", "IsotopeLabelType")
    setnames(tmp_long, "Raw file", "Raw.file")
    levels(tmp_long$IsotopeLabelType) = c("L", "H")
    tmp_long[!is.na(tmp_long$Intensity) && tmp_long$Intensity < 1, ]$Intensity = NA
    write.table(tmp_long, file = output, sep = "\t", quote = F, row.names = F, col.names = T)
}

#' @title Concatenate MaxQuant search result files
#' @description Concatenates multiple MaxQuant files into one  big file.
#' @param filenames The filepaths to the MaxQuant searched results files (txt tab delimited file).
#' @param output The filepath to the intended output file (txt tab delimited file).
#' @keywords concatenate join merge
#' MQutil.concat()
MQutil.concat = function(filenames, output) {
    library(data.table)
    files = Sys.glob(filenames)
    cat(sprintf("\tPROCESSING:\n\t%s\n", paste(files, collapse = "\n\t")))
    
    res = NULL
    unique_files = c()
    
    for (f in 1:length(files)) {
        file = files[f]
        tmp = fread(file, stringsAsFactors = F, colClasses = c(Intensity = "character"), integer64 = "double")
        tmp$Intensity = as.numeric(tmp$Intensity)
        # tmp$Intensity L = as.numeric(tmp[,'Intensity L',with=F])
        
        unique_files_current = unique(tmp[["Raw file"]])
        if (!is.null(intersect(unique_files_current, unique_files)) && length(intersect(unique_files_current, unique_files)) > 
            0) 
            cat(sprintf("\tWARNING DUPLICATE RAW FILE ENTRIES IN FILE %s:\t%s\n", file, paste(intersect(unique_files_current, unique_files), 
                collapse = ",")))
        select_colnames = grep("Raw file|Intensity|Proteins|Modifications|Sequence|Modified sequence|Charge|Protein group IDs|id|Retention time|Reverse|Contaminant", 
            colnames(tmp), ignore.case = F)
        tmp = tmp[, select_colnames, with = F]
        res = rbind(res, tmp, fill = T)
        unique_files = c(unique_files, unique_files_current)
    }
    
    cat(sprintf("\tWRITING\t%s\n", output))
    write.table(res, file = output, eol = "\n", sep = "\t", quote = F, row.names = F, col.names = T)
    cat(sprintf("\tWRITING\t%s\n", gsub(".txt", "-keys.txt", output)))
    # print(unique_files)
    write.table(unique_files, file = gsub(".txt", "-keys.txt", output), eol = "\n", sep = "\t", quote = F, row.names = F, col.names = F)
}

MQutil.getKeys = function(filename, output) {
    library(data.table)
    file = Sys.glob(filename)
    cat(sprintf("\tPROCESSING:\n\t%s\n", paste(file, collapse = "\n\t")))
    # tmp = data.table(read.delim(file, stringsAsFactors=F))
    tmp = fread(file, stringsAsFactors = F, integer64 = "double")
    write.table(unique(tmp$Raw.file), file = output, eol = "\n", sep = "\t", quote = F, row.names = F, col.names = F)
    cat(sprintf("\tWRITTEN\t%s\n", output))
}


#' @title Convert Protein format to Site format
#' @description This function is a preprocessing function that is used with PTM data. If you want to run a site specific analysis, before running MSStats on a UB/PH/AC set, you need to convert the evidence file into a format that MSStats will be able to diffentiate the sites with. This outputs a new file that should be used as the input file for your MSStats group conparison analysis.
#' @param maxq_file The filepath to the MaxQuant searched results file (txt tab delimited file).
#' @param output The filepath to the intended output file (txt tab delimited file).
#' @param ref_proteome_file Filepath to the fasta file containing the reference proteom.
#' @param mod_type The kind of PTM site analysis to be perfomed. Pick between ub (default), ph, and ac.
#' @keywords site, ph, ub, ac, ubiquitination, phosphorylation, acetylation
#' MQutil.ProteinToSiteConversion()
MQutil.ProteinToSiteConversion <- function(maxq_file, ref_proteome_file, output_file, mod_type = "ub") {
    
    if (mod_type == "ub") {
        maxq_mod_residue = "K\\(gl\\)"
        mod_residue = "K"
    } else if (mod_type == "ph") {
        maxq_mod_residue = "(S|T|Y)\\(ph\\)"
        mod_residue = "S|T|Y"
    } else if (mod_type == "ac") {
        maxq_mod_residue = "K\\(ac\\)"
        mod_residue = "K"
    }
    
    ## read in ref. proteome
    ref_proteome = read.fasta(file = ref_proteome_file, seqtype = "AA", as.string = T, set.attributes = TRUE, legacy.mode = TRUE, 
        seqonly = FALSE, strip.desc = FALSE)
    
    ###################### make mod-site index
    
    p_seqs = c()
    p_names = c()
    p_annots = c()
    
    for (e in ref_proteome) {
        p_seqs = c(p_seqs, e[1])
        p_names = c(p_names, attr(e, "name"))
        p_annots = c(p_annots, attr(e, "Annot"))
    }
    
    ref_table = data.table(names = p_names, annots = p_annots, seqs = p_seqs)
    ref_table[, `:=`(uniprot_ac, gsub("([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)", "\\2", names))]
    
    # get all indicies/locations of the mod_resudue S|T|Y
    indices = lapply(ref_table$seqs, function(x) as.vector(str_locate_all(x, pattern = mod_residue)[[1]][, 1]))
    # list which residue it actually finds at each location
    ptm_sites = lapply(ref_table$seqs, function(x) as.vector(str_match_all(x, pattern = mod_residue)))
    # get the number of residue locations per string
    lengths = unlist(lapply(indices, FUN = length))
    # repeate the uniprot_ac as many times as a residue (S|T|Y) is found
    keys = rep(ref_table$uniprot_ac, lengths)
    # combine the list of the exploded proteins (above) and add in which sites are found as well as the location the residue was
    # found in
    protein_indices = data.table(uniprot_ac = keys, ptm_site = unlist(ptm_sites), res_index = unlist(indices))
    
    ################################# map mod sites in data to index
    
    ## read in maxq. data
    maxq_data = fread(maxq_file, integer64 = "double")
    # remove contaminants, keep unique sequences, fix names
    maxq_data = maxq_data[grep("CON__|REV__", maxq_data$Proteins, invert = T), ]
    unique_peptides_in_data = unique(maxq_data[, c("Proteins", "Modified sequence"), with = F])
    setnames(unique_peptides_in_data, "Modified sequence", "sequence")
    
    mod_sites = c()
    mod_seqs = c()
    
    for (i in 1:nrow(unique_peptides_in_data)) {
        entry = unique_peptides_in_data[i, ]
        peptide_seq = entry$sequence
        ## cleanup the sequence (removing all modifications) for matching the protein sequence
        peptide_seq_clean = gsub("[a-z,0-9,\\(,\\),_]", "", peptide_seq)
        mod_sites_in_peptide = str_locate_all(string = peptide_seq, pattern = maxq_mod_residue)[[1]][, 1]
        
        if (length(mod_sites_in_peptide) > 0) {
            uniprot_acs = entry$Proteins
            # separates the ambiguous cases (;) and appends the site info to all the proteins
            uniprot_acs = str_split(string = uniprot_acs, pattern = ";")[[1]]
            
            for (uac in uniprot_acs) {
                # find the protein's full sequence based on the uniprot_ac
                protein_seq = ref_table[uniprot_ac == uac, ]$seqs
                if (length(protein_seq) > 0) {
                  ## get the position of the peptide in the protein sequence
                  peptide_index_in_protein = str_locate(protein_seq, peptide_seq_clean)[[1]][1]
                  
                  for (m in 1:length(mod_sites_in_peptide)) {
                    mod_site = mod_sites_in_peptide[m]
                    peptide_seq_before_site = str_sub(peptide_seq, 1, mod_site - 1)
                    ## count all AA (not counting all modifications) before the modification to get the relative position of the modification in the
                    ## peptide sequence
                    residues_before_site = str_count(string = peptide_seq_before_site, pattern = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y")
                    # location of the modified site (end index)
                    mod_site_index_in_protein = peptide_index_in_protein + residues_before_site
                    protein_mod_sites = protein_indices[uniprot_ac == uac, ]
                    if (!is.na(mod_site_index_in_protein)) {
                      # cat(sprintf('%s\n',mod_site_id))
                      mod_res = protein_mod_sites[res_index == mod_site_index_in_protein, ptm_site]
                      mod_site_id = sprintf("%s_%s%s", uac, str_sub(protein_seq, mod_site_index_in_protein, mod_site_index_in_protein), 
                        mod_site_index_in_protein)
                      mod_sites = c(mod_sites, mod_site_id)
                      mod_seqs = c(mod_seqs, peptide_seq)
                      stopifnot(length(mod_sites) == length(mod_seqs))
                    } else {
                      cat(sprintf("MISMATCH\t%s\n\tPEPTIDE_SEQ\t%s\n\tMOD_SITE\t%s\n\tPEPTIDE_IDX_IN_PROTEIN\t%s\n\tRESIDUES_BEFORE_SITE\t%s\n\tPROTEIN_SEQ\t%s\n", 
                        mod_site_id, peptide_seq, mod_site, peptide_index_in_protein, residues_before_site, protein_seq))
                    }
                  }
                }
            }
        }
    }
    
    mod_site_mapping = data.table(mod_sites, mod_seqs)
    mod_site_mapping_agg = aggregate(mod_sites ~ mod_seqs, mod_site_mapping, FUN = function(x) paste(x, collapse = ","))
    
    setnames(maxq_data, "Modified sequence", "mod_seqs")
    unmapped_mod_seqs = maxq_data[!(mod_seqs %in% mod_site_mapping_agg$mod_seqs) & grepl("(gl)", mod_seqs) & !grepl("REV__|CON__", 
        Proteins), ]
    unmapped_mod_seqs = unique(unmapped_mod_seqs[, c("mod_seqs", "Proteins"), with = F])
    if (dim(unmapped_mod_seqs)[1] > 0) {
        cat("UNABLE TO MAP\n")
        print(unmapped_mod_seqs)
    } else {
        cat("ALL SEQUENCES MAPPED\n")
    }
    
    final_data = merge(maxq_data, mod_site_mapping_agg, by = "mod_seqs")
    setnames(final_data, c("Proteins", "mod_sites", "mod_seqs"), c("Proteins_ref", "Proteins", "Modified sequence"))
    write.table(final_data, file = output_file, eol = "\n", sep = "\t", quote = F, row.names = F, col.names = T)
    
    ## write a mapping table
    protein_seq_mapping = unique(maxq_data[, c("Proteins", "mod_seqs"), with = F])
    setnames(protein_seq_mapping, "Proteins", "Protein")
    mapping_table = merge(protein_seq_mapping, mod_site_mapping_agg, by = "mod_seqs", all = T)
    write.table(mapping_table, file = gsub(".txt", "-mapping.txt", output_file), eol = "\n", sep = "\t", quote = F, row.names = F, 
        col.names = T)
}


#' @title Annotate MSStats results file.
#' @description Annotates the proteins in the results or results-wide files after they've been through the MSStats pipeline. Multiple species can be searched at once, simply separate them by a '-'. (eg. human-mouse).
#' @param input_file The filepath to the MaxQuant searched results file (txt tab delimited file). This can either be in either long or wide format.
#' @param output_file The filepath to the intended output file (txt tab delimited file).
#' @param uniprot_ac_col Name of the column containing the uniprot ID's. Defaults to 'Protein'.
#' @param uniprot_dir Directory location of where the uniprot database flatfiles are located.
#' @param species The species contained in this search of uniprot files. Can include multiple species separated by a '-' .
#' @keywords annotate annotation uniprot
#' MQutil.annotate()
MQutil.annotate = function(input_file = opt$input, output_file = opt$output, uniprot_ac_col = "Protein", group_sep = ";", uniprot_dir = "~/Box Sync/db/mist/", 
    species = "HUMAN") {
    cat(">> ANNOTATING\n")
    results = read.delim(input_file, stringsAsFactors = F, sep = "\t")
    
    # remove unnamed proteins that are listed as ''
    if (length(which(results$Protein == "")) > 0) 
        results <- results[-which(results$Protein == ""), ]
    
    # read in all the annotation files from the uniprot_db directory
    species_split = unlist(strsplit(species, "-"))
    Uniprot = NULL
    for (org in species_split) {
        cat(sprintf("\tLOADING %s\n", org))
        tmp = read.delim2(sprintf("%s/uniprot_protein_descriptions_%s.txt", uniprot_dir, org), stringsAsFactors = F, quote = "")
        if (is.null(Uniprot)) {
            Uniprot = as.data.frame(tmp)
        } else {
            Uniprot = rbind(Uniprot, tmp)
        }
    }
    
    # get list of all unique prey entries in this file. Keep 'group_sep' in mind.
    preys <- unique(results$Protein)
    preys <- preys.original <- data.frame(prey = preys, idx = 1:length(preys), stringsAsFactors = F)
    # split apart all the preys and index them so we can piece them back together when merging
    preys <- do.call(rbind, apply(preys, 1, function(y) {
        data.frame(prey = unlist(strsplit(y[1], ";")), idx = as.numeric(y[2]), stringsAsFactors = F)
    }))
    
    
    if (!any(Uniprot$Entry %in% preys$prey)) {
        return(warning("NO ANNOTATIONS FOUND IN THE DATABASES. PLEASE CHECK THAT THE CORRECT DATABASES WERE ENTERED,\n  AND THAT THE PROTEIN COLUMN CONTAINS THE CORRECT NAMES TO SEARCH."))
    }
    # annotate all the preys wiht the Uniprot info
    preys <- merge(preys, Uniprot[, c("Entry", "Entry.name", "Protein.names", "Gene.names")], by.x = "prey", by.y = "Entry", all.x = T)
    # aggregate all the preys on the indexes so we can merge with the original data
    
    # merge protein name
    tmp <- aggregate(data = preys[, c("prey", "idx", "Entry.name")], . ~ idx, paste, collapse = ";")
    names(tmp) <- c("idx", "uniprot_ac", "Protein_name")
    preys.new <- merge(preys.original, tmp, by = "idx", all.x = T)
    
    # merge protein description
    tmp <- aggregate(data = preys[, c("idx", "Protein.names")], . ~ idx, paste, collapse = ";")
    names(tmp) <- c("idx", "Protein_desc")
    preys.new <- merge(preys.new, tmp, by = "idx", all.x = T)
    
    # merge protein description
    preys$Gene.names <- gsub(" .*", "", preys$Gene.names)
    tmp <- aggregate(data = preys[, c("idx", "Gene.names")], . ~ idx, paste, collapse = ";")
    names(tmp) <- c("idx", "Gene.names")
    preys.new <- merge(preys.new, tmp, by = "idx", all.x = T)
    
    # merge the annotations all back into the original data
    results_out <- merge(results, preys.new, by.x = "Protein", by.y = "prey", all.x = T)
    results_out$idx = c()
    
    # alert user of any unmapped proteins
    unmapped = unique(results_out[is.na(results_out$uniprot_ac), "Protein"])
    cat("UNMAPPED PROTEINS\t", length(unmapped), "\n")
    cat("\t", paste(unmapped, collapse = "\n\t"), "\n")
    write.table(results_out, file = output_file, sep = "\t", quote = F, row.names = F, col.names = T)
    cat(">> ANNOTATING COMPLETE!\n")
}


#' #' @title Convert default MSStat results file to Wide format.
#' #' @description Converts the normal MSStats output file into 'wide' format where each row represents a protein's results, and each column represents the comparison made by MSStats. The fold change and p-value of each comparison will be it's own column.
#' #' @param input_file The filepath to the MaxQuant searched results file (txt tab delimited file). This can either be in either long or wide format.
#' #' @param save_file The filepath to the intended output file (txt tab delimited file).
#' #' @param return_results Whether to return the results at the end of the function (default =F).
#' #' @keywords wide MSStats results
#' #' MQutil.resultsWide()
#' MQutil.resultsWide = function(input_file, save_file=F, return_results=F){
#'   input = checkIfFile(input_file)
#'   input_l = melt(data = input[,c('Protein', 'Label','log2FC','adj.pvalue'),with=F],id.vars=c('Protein', 'Label'))
#' 
#'   ## then cast to get combinations of LFCV/PVAl and Label as columns
#'   input_w = dcast.data.table( Protein ~ Label+variable, data=input_l, value.var=c('value'))
#'   if(save_file){
#'     write.table(input_w, file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
#'   }
#'   if(return_results){
#'     return(input_w)
#'   }
#' }


#' @title Map the results back into sites format after MSStats PTM site analysis.
#' @description Used with PTM datasets. Map back the sites to correct proteins after MSStats analysis. This file is created previously when running the conver-sites function.
#' @param input_file The filepath to the MaxQuant searched results file (txt tab delimited file). This can either be in either long or wide format.
#' @param mapping_file The filepath to the mapping file (ending in '-mapping.txt')
#' @param output_file The filepath to the intended output file (txt tab delimited file).
#' @keywords map mapping MSStats results ptm sites
#' MQutil.mapSitesBack()
MQutil.mapSitesBack = function(input_file, mapping_file, output_file) {
    input = fread(input_file, integer64 = "double")
    setnames(input, "Protein", "mod_sites")
    mapping = fread(mapping_file, integer64 = "double")
    mapping = unique(mapping[!is.na(mod_sites), c("Protein", "mod_sites"), with = F])
    mapping = aggregate(Protein ~ mod_sites, data = mapping, FUN = function(x) paste(x, collapse = ";"))
    out = merge(input, mapping, by = "mod_sites", all.x = T)
    write.table(out[, c(ncol(out), 1:(ncol(out) - 1)), with = F], file = output_file, eol = "\n", sep = "\t", quote = F, row.names = F, 
        col.names = T)
}


MQutil.plotHeat = function(input_file, output_file, labels = NULL, names = "Protein", cluster_cols = F, display = "log2FC", row_names = "Protein", 
    col_names = "Label") {
    heat_data = fread(input_file, integer64 = "double")
    # heat_data = mss_F[,c('uniprot_id','Label','log2FC')]
    
    ## good
    if (display == "log2FC") {
        heat_data_w = dcast(names ~ Label, data = heat_data, value.var = "log2FC")
    } else if (display == "pvalue") {
        heat_data$adj.pvalue = -log10(heat_data$adj.pvalue + 10^-16)
        heat_data_w = dcast(names ~ Label, data = heat_data, value.var = "adj.pvalue")
    }
    
    ## try
    
    # gene_names = uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
    rownames(heat_data_w) = heat_data_w$names
    heat_data_w = heat_data_w[, -1]
    heat_data_w[is.na(heat_data_w)] = 0
    max_val = ceiling(max(heat_data_w))
    min_val = floor(min(heat_data_w))
    extreme_val = max(max_val, abs(min_val))
    if (extreme_val%%2 != 0) 
        extreme_val = extreme_val + 1
    bin_size = 2
    signed_bins = (extreme_val/bin_size)
    colors_neg = rev(colorRampPalette(brewer.pal("Blues", n = extreme_val/bin_size))(signed_bins))
    colors_pos = colorRampPalette(brewer.pal("Reds", n = extreme_val/bin_size))(signed_bins)
    colors_tot = c(colors_neg, colors_pos)
    
    if (is.null(labelOrder)) {
        pheatmap(heat_data_w, scale = "none", cellheight = 10, cellwidth = 10, filename = out_file, color = colors_tot, breaks = seq(from = -extreme_val, 
            to = extreme_val, by = bin_size), cluster_cols = cluster_cols, fontfamily = "mono")
    } else {
        heat_data_w = heat_data_w[, labelOrder]
        pheatmap(heat_data_w, scale = "none", cellheight = 10, cellwidth = 10, filename = out_file, color = colors_tot, breaks = seq(from = -extreme_val, 
            to = extreme_val, by = bin_size), cluster_cols = cluster_cols, fontfamily = "mono")
    }
    
    heat_data_w
}

## still need to make fileType, conditions, cluster_cols, display configurable through command line
MQutil.plotHeatmap = function(input_file, output_file, fileType = "l", labels = "*", cluster_cols = F, display = "log2FC", lfc_lower = -2, 
    lfc_upper = 2, FDR = 0.05) {
    ## read input
    input = read.delim(input_file, stringsAsFactors = F)
    
    ## select data points by LFC & FDR criterium in single condition and adding corresponding data points from the other conditions
    sign_hits = significantHits(input, labels = labels, LFC = c(lfc_lower, lfc_upper), FDR = FDR)
    sign_labels = unique(sign_hits$Label)
    cat(sprintf("\tSELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n", lfc_lower, lfc_upper, FDR, nrow(sign_hits)/length(sign_labels)))
    
    ## REPRESENTING RESULTS AS HEATMAP plot heat map for all contrasts
    if (any(grepl("uniprot_genename", colnames(sign_hits)))) {
        heat_labels = paste(sign_hits$Protein, sign_hits$uniprot_genename, sep = " ")
    } else {
        heat_labels = sign_hits$Protein
    }
    
    heat_labels = gsub("\\sNA$", "", heat_labels)
    heat_data_w = plotHeat(mss_F = sign_hits, out_file = output_file, names = heat_labels, cluster_cols = cluster_cols, display = display)
}

MQutil.simplify = function(input_file, output_file) {
    input = fread(input_file, integer64 = "double")
    output = simplifyOutput(input)
    write.table(output, file = output_file, sep = "\t", quote = F, row.names = F, col.names = T, eol = "\n")
}


#' @title Convert MaxQuant file into a format to work with SAINT
#' @description Converts the MaxQuant evidence file to the 3 required files for SAINTexpress. One can choose to either use the spectral counts or the intensities for the analysis.
#' @param input_file The filepath to the MaxQuant searched results file (txt tab delimited file). This can either be in either long or wide format.
#' @param keys_file The filepath to the keys file used in the MSStats analysis.
#' @param ref_proteome_file The filepath to the reference proteom (txt tab delimited file).
#' @param quant_variable Pick which quantitative value to use: 'spectral_counts' or 'ms1'
#' @param output_file Path of where to output the results to.
#' @keywords SAINT spectral spectral_counts ms1 intensity
#' MQutil.MaxQToSaint()
MQutil.MaxQToSaint = function(data_file, keys_file, ref_proteome_file, quant_variable = "spectral_count", output_file) {
    cat(">> CONVERTING TO SAINT FORMAT\n")
    
    data = fread(data_file, integer64 = "double")
    keys = fread(keys_file, integer64 = "double")
    
    ## write baits in format hIP101-10 PV_2C_co_uni T
    
    saint_baits = keys[, c("BioReplicate", "Condition", "SAINT"), with = F]
    
    ## write interactions in format hIP101-10 PV_2C_co_uni Q9NTG7 1
    
    tryCatch(setnames(data, "Raw file", "RawFile"), error = function(e) cat("Raw.file not found\n"))
    tryCatch(setnames(keys, "Raw.file", "RawFile"), error = function(e) cat("Raw.file not found\n"))
    
    cat("\tVERIFYING DATA AND KEYS\n")
    if (any(!c("RawFile", "IsotopeLabelType", "Condition", "BioReplicate", "Run", "SAINT") %in% colnames(keys))) {
        stop("COLNAMES IN KEYS NOT CONFORM TO SCHEMA\n\tRawFile\tIsotopeLabelType\tCondition\tBioReplicate\tRun\tSAINT\n")
    }
    if (!"IsotopeLabelType" %in% colnames(data)) 
        data[, `:=`(IsotopeLabelType, "L")]
    data = mergeMaxQDataWithKeys(data, keys, by = c("RawFile", "IsotopeLabelType"))
    data_f = filterMaxqData(data)
    data_f = removeMaxQProteinGroups(data_f)  ## do we want this or not?
    
    cat("\tAGGREGATING ON", quant_variable, "VALUES...\n")
    ## aggregate over technical replicates if necessary
    if (quant_variable == "spectral_count") {
        setnames(data_f, "MS/MS Count", "spectral_counts")
        data_f_agg = aggregate(spectral_counts ~ BioReplicate + Condition + Proteins + Sequence + Charge, data = data_f, FUN = max)
        data_f_agg = aggregate(spectral_counts ~ BioReplicate + Condition + Proteins, data = data_f_agg, FUN = sum)
    } else if (quant_variable == "ms1") {
        data_f_agg = aggregate(Intensity ~ BioReplicate + Condition + Proteins + Sequence + Charge, data = data_f, FUN = max)
        data_f_agg = aggregate(Intensity ~ BioReplicate + Condition + Proteins, data = data_f_agg, FUN = sum)
    } else {
        stop("ERROR!! Wrong value for variable to quantify. Please use 'spectral_count' or 'ms1'.")
    }
    
    ## IP name, bait name, prey name, and spectral counts or intensity values
    saint_interactions = data_f_agg
    
    ## write preys in format Q9NTG7 43573.5 Q9NTG7
    
    ref_proteome = read.fasta(file = ref_proteome_file, seqtype = "AA", as.string = T, set.attributes = TRUE, legacy.mode = TRUE, 
        seqonly = FALSE, strip.desc = FALSE)
    p_lengths = c()
    p_names = c()
    for (e in ref_proteome) {
        p_lengths = c(p_lengths, nchar(e[1]))
        p_names = c(p_names, attr(e, "name"))
    }
    ref_table = data.table(names = p_names, lengths = p_lengths)
    ref_table[, `:=`(uniprot_ac, gsub("([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)", "\\2", names))]
    ref_table[, `:=`(uniprot_id, gsub("([a-z,0-9,A-Z]+\\|{1})([a-z,0-9,A-Z]+\\|{1})([A-Z,a-z,0-9,_]+)", "\\3", names))]
    
    unique_preys = data.table(uniprot_ac = unique(data_f_agg$Proteins))
    saint_preys = ref_table[, c("uniprot_ac", "lengths", "uniprot_id"), with = F]
    saint_preys = merge(unique_preys, saint_preys, by = "uniprot_ac", all.x = T)
    missing_lengths = nrow(saint_preys[is.na(saint_preys$uniprot_id), ])
    saint_preys[is.na(saint_preys$uniprot_id), ]$uniprot_id = saint_preys[is.na(saint_preys$uniprot_id), ]$uniprot_ac
    if (missing_lengths > 0) {
        cat(sprintf("\tWARNING!\tCOMPUTING %s MISSING LENGTHS WITH THE MEDIAN LENGTH FROM THE DATASET\n", missing_lengths))
        saint_preys[is.na(saint_preys$lengths), ]$lengths = median(saint_preys$lengths, na.rm = T)
    }
    
    
    # Check if output directory exists and create it if not.
    if (!dir.exists(dirname(output_file))) {
        dir.create(dirname(output_file), recursive = T)
    }
    
    ## WRITE
    write.table(saint_baits, file = gsub(".txt", "-saint-baits.txt", output_file), eol = "\n", sep = "\t", quote = F, row.names = F, 
        col.names = F)
    write.table(saint_preys, file = gsub(".txt", "-saint-preys.txt", output_file), eol = "\n", sep = "\t", quote = F, row.names = F, 
        col.names = F)
    write.table(saint_interactions, file = gsub(".txt", "-saint-interactions.txt", output_file), eol = "\n", sep = "\t", quote = F, 
        row.names = F, col.names = F)
}


MQutil.dataPlots = function(input_file, output_file) {
    data_mss = fread(input_file, integer64 = "double")
    unique_subjects = unique(data_mss$PROTEIN)
    condition_length = length(unique(data_mss$GROUP_ORIGINAL))
    min_abu = min(data_mss$ABUNDANCE, na.rm = T)
    max_abu = max(data_mss$ABUNDANCE, na.rm = T)
    
    pdf(output_file, width = condition_length * 1.5, height = 3)
    
    cat("PRINTING CONDITION PLOTS\n")
    for (subject in unique_subjects) {
        subject_data = data_mss[PROTEIN == subject, ]
        cat(sprintf("\t%s\n", subject))
        p = ggplot(data = subject_data, aes(x = SUBJECT_ORIGINAL, y = ABUNDANCE, colour = FEATURE))
        p = p + geom_point(size = 2) + facet_wrap(facets = ~GROUP_ORIGINAL, drop = T, scales = "free_x", ncol = condition_length) + 
            ylim(min_abu, max_abu) + theme(axis.text.x = element_text(angle = -90, hjust = 1)) + guides(colour = FALSE) + xlab(NULL) + 
            ggtitle(subject)
        print(p)
    }
    dev.off()
}


#' @title Get Spectral Counts
#' @description Outputs the spectral counts from the MaxQuant evidence file.
#' @param input_file The filepath to the MaxQuant searched results file (txt tab delimited file). This can either be in either long or wide format.
#' @param keys_file The filepath to the keys file used in the MSStats analysis.
#' @param output_file The filepath to the intended output file (txt tab delimited file).
#' @keywords spectral spectral_counts spectral
#' MQutil.spectralCounts()
MQutil.spectralCounts = function(input_file, keys_file, output_file) {
    data = fread(input_file, integer64 = "double")
    keys = fread(keys_file, integer64 = "double")
    
    tryCatch(setnames(data, "Raw file", "RawFile"), error = function(e) cat("Raw file not found. Try searching Raw.file instead\n"))
    tryCatch(setnames(keys, "Raw.file", "RawFile"), error = function(e) cat("Raw file not found. Try searching Raw file instead\n"))
    
    cat("\tVERIFYING DATA AND KEYS\n")
    if (!"IsotopeLabelType" %in% colnames(data)) 
        data[, `:=`(IsotopeLabelType, "L")]
    data = mergeMaxQDataWithKeys(data, keys, by = c("RawFile", "IsotopeLabelType"))
    data_sel = data[, c("Proteins", "Condition", "BioReplicate", "Run", "MS/MS Count"), with = F]
    setnames(data_sel, "MS/MS Count", "spectral_counts")
    data_sel = aggregate(spectral_counts ~ Proteins + Condition + BioReplicate + Run, data = data_sel, FUN = sum)
    data_sel = data.frame(data_sel, bait_name = paste(data_sel$Condition, data_sel$BioReplicate, data_sel$Run, sep = "_"))
    write.table(data_sel[, c("bait_name", "Proteins", "spectral_counts")], file = output_file, eol = "\n", sep = "\t", quote = F, 
        row.names = F, col.names = T)
}


#' @title Convert MaxQuant to MIST format
#' @description Converts MaxQuant evidence file into a file format compatible with the MiST pipeline using MS/MS.Count. Note that this is the MiST data file, and that an additional keys file will have to be constructed before running MiST. Multiple species can be searched at once, simply separate them by a '-'. (eg. HUMAN-MOUSE).
#' @param input_file The filepath to the MaxQuant searched results file (txt tab delimited file). This can either be in either long or wide format.
#' @param keys_file The filepath to the keys file used in the MSStats analysis.
#' @param output_file The filepath to the intended output file (txt tab delimited file).
#' @param species The species to in the data to be used to get weights from. Multiple species can be listed, separating with a '-'. (eg. HUMAN-MOUSE)
#' @param uniprot_dir The directory where the uniprot files are for the different species.
#' @keywords spectral MIST mist
#' MQutil.MISTformat()
# Convert MaxQuant file into a Protein Prospector like format to run through the Mist pipeline
MQutil.MISTformat = function(input_file, keys_file, output_file, species = "HUMAN", uniprot_dir = "~/Box Sync/db/mist/") {
    cat("\tREADING IN DATA AND KEYS\n")
    
    data <- data.table(read.delim(input_file, stringsAsFactors = F))
    keys = data.table(read.delim(keys_file, stringsAsFactors = F))
    
    tryCatch(setnames(data, "Raw file", "RawFile"), error = function(e) cat("Raw file in evidence not found: trying Raw.file instead\n"))
    tryCatch(setnames(data, "Raw.file", "RawFile"), error = function(e) cat("Raw.file in evidence not found: trying Raw file instead\n"))
    tryCatch(setnames(keys, "Raw file", "RawFile"), error = function(e) cat("Raw file in keys not found: trying Raw.file instead\n"))
    tryCatch(setnames(keys, "Raw.file", "RawFile"), error = function(e) cat("Raw.file in keys not found: trying Raw file instead\n"))
    tryCatch(setnames(data, "MS/MS Count", "ms_spectral_counts"), error = function(e) cat("MS/MS Count column not found in evidence file: trying MS.MS.Count instead\n"))
    tryCatch(setnames(data, "MS.MS.Count", "ms_spectral_counts"), error = function(e) cat("MS.MS.Count column not found in evidence file: trying MS/MS Count instead\n"))
    
    cat("\t\nVERIFYING DATA AND KEYS\n")
    if (!"IsotopeLabelType" %in% colnames(data)) 
        data[, `:=`(IsotopeLabelType, "L")]
    data = mergeMaxQDataWithKeys(data, keys, by = c("RawFile", "IsotopeLabelType"))
    data_sel = data[, c("Proteins", "Condition", "BioReplicate", "Run", "RawFile", "ms_spectral_counts"), with = F]
    
    data_sel = aggregate(ms_spectral_counts ~ Proteins + Condition + BioReplicate + Run + RawFile, data = data_sel, FUN = sum)
    data_sel = data.frame(data_sel, bait_name = paste(data_sel$Condition, data_sel$BioReplicate, data_sel$Run, sep = "_"))
    
    # clean up proteins & annotate ~~~~~~~~~~~~~~~~~~~~~~~ remove CON's
    if (length(grep("^CON__", data_sel$Proteins)) > 0) 
        data_sel = data_sel[-grep("^CON__", data_sel$Proteins), ]
    if (length(grep("^REV__", data_sel$Proteins)) > 0) 
        data_sel = data_sel[-grep("^REV__", data_sel$Proteins), ]
    # remove the party sets
    if (length(grep(";", data_sel$Proteins)) > 0) 
        data_sel = data_sel[-grep(";", data_sel$Proteins), ]  # NOTE!!! We lose a lot of entries this way... :\
    # keep only uniprot id data_sel$uniprot_id = gsub('^.*\\|','', data_sel$Proteins)
    data_sel$Proteins = gsub("(^.*\\|)([A-Z0-9]+)(\\|.*$)", "\\2", data_sel$Proteins)
    # remove blank protein names
    if (any(data_sel$Proteins == "")) {
        data_sel <- data_sel[-which(data_sel$Proteins == ""), ]
    }
    
    # Add this column that it is required for the preprocessing done by mist (this column is generated by Prospector and it is used
    # to eleminate peptides. In this case, does not apply but it has to be there)
    data_sel$ms_unique_pep = ""
    # re-order
    data_sel <- data_sel[, c("RawFile", "Proteins", "ms_unique_pep", "ms_spectral_counts")]
    # RENAMING!
    names(data_sel) = c("id", "ms_uniprot_ac", "ms_unique_pep", "ms_spectral_counts")
    
    # remove interactions with ms_spectral_counts=0
    if (any(data_sel$ms_spectral_counts == 0)) {
        data_sel <- data_sel[-which(data_sel$ms_spectral_counts == 0), ]
    }
    
    # annotate proteins and add Masses for Mist
    species_split = unlist(strsplit(species, "-"))
    Uniprot = NULL
    for (org in species_split) {
        cat(sprintf("\tLOADING %s\n", org))
        tmp = read.delim2(sprintf("%s/uniprot_protein_descriptions_%s.txt", uniprot_dir, org), stringsAsFactors = F, quote = "")
        if (is.null(Uniprot)) {
            Uniprot = as.data.frame(tmp)
        } else {
            Uniprot = rbind(Uniprot, tmp)
        }
    }
    results_annotated = merge(data_sel, Uniprot, all.x = T, by.x = "ms_uniprot_ac", by.y = "Entry")
    
    # Adding the Mass column: it is required for MIST for preprocessing!  For that, we will calculate the mass of the whole protein
    # just taking the average molecular weight of an amino acid: 110Da
    results_annotated$Mass <- results_annotated$Length * 110
    
    write.table(results_annotated, file = output_file, eol = "\n", sep = "\t", quote = F, row.names = F, col.names = T)
    cat("\nMIST FILES CREATED! Have a nice day :)\n")
    
    keysout <- subset(keys, select = c(RawFile, Condition))
    keysname <- "keys_mist.txt"
    write.table(keysout, file = keysname, eol = "\n", sep = "\t", quote = F, row.names = F, col.names = F)
    cat("MIST keys FILE ALSO CREATED. Truly enjoy your day ;-)\n\n")
}

# Convert MaxQuant file into a Protein Prospector like format to run through the Mist pipeline
MQutil.MISTINTformat = function(input_file, keys_file, output_file, species = "HUMAN", uniprot_dir = "~/Box Sync/db/mist/") {
    cat("\n\n>>Generating input files for MIST using INTENSITY VALUES\n\n")
    cat("\tREADING IN DATA AND KEYS\n")
    
    data <- data.table(read.delim(input_file, stringsAsFactors = F))
    keys = data.table(read.delim(keys_file, stringsAsFactors = F))
    
    tryCatch(setnames(data, "Raw file", "RawFile"), error = function(e) cat("Raw file in evidence not found: trying Raw.file instead\n"))
    tryCatch(setnames(data, "Raw.file", "RawFile"), error = function(e) cat("Raw.file in evidence not found: trying Raw file instead\n"))
    tryCatch(setnames(keys, "Raw file", "RawFile"), error = function(e) cat("Raw file in keys not found: trying Raw.file instead\n"))
    tryCatch(setnames(keys, "Raw.file", "RawFile"), error = function(e) cat("Raw.file in keys not found: trying Raw file instead\n"))
    tryCatch(setnames(data, "Intensity", "ms_intensity"), error = function(e) stop("\n\nINTENSITY NOT FOUND IN THE evidence FILE!!\n\n"))
    
    
    cat("\n\tVERIFYING DATA AND KEYS\n")
    if (!"IsotopeLabelType" %in% colnames(data)) 
        data[, `:=`(IsotopeLabelType, "L")]
    data = mergeMaxQDataWithKeys(data, keys, by = c("RawFile", "IsotopeLabelType"))
    data_sel = data[, c("Proteins", "Condition", "BioReplicate", "Run", "RawFile", "ms_intensity"), with = F]
    
    data_sel = aggregate(ms_intensity ~ Proteins + Condition + BioReplicate + Run + RawFile, data = data_sel, FUN = sum)
    data_sel = data.frame(data_sel, bait_name = paste(data_sel$Condition, data_sel$BioReplicate, data_sel$Run, sep = "_"))
    
    # clean up proteins & annotate ~~~~~~~~~~~~~~~~~~~~~~~ remove CON's
    if (length(grep("^CON__", data_sel$Proteins)) > 0) 
        data_sel = data_sel[-grep("^CON__", data_sel$Proteins), ]
    if (length(grep("^REV__", data_sel$Proteins)) > 0) 
        data_sel = data_sel[-grep("^REV__", data_sel$Proteins), ]
    # remove the party sets
    if (length(grep(";", data_sel$Proteins)) > 0) 
        data_sel = data_sel[-grep(";", data_sel$Proteins), ]  # NOTE!!! We lose a lot of entries this way... :\
    # keep only uniprot id data_sel$uniprot_id = gsub('^.*\\|','', data_sel$Proteins)
    data_sel$Proteins = gsub("(^.*\\|)([A-Z0-9]+)(\\|.*$)", "\\2", data_sel$Proteins)
    # remove blank protein names
    if (any(data_sel$Proteins == "")) {
        data_sel <- data_sel[-which(data_sel$Proteins == ""), ]
    }
    
    # Add this column that it is required for the preprocessing done by mist (this column is generated by Prospector and it is used
    # to eleminate peptides. In this case, does not apply but it has to be there)
    data_sel$ms_unique_pep = ""
    # re-order
    data_sel <- data_sel[, c("RawFile", "Proteins", "ms_unique_pep", "ms_intensity")]
    # RENAMING!
    names(data_sel) = c("id", "ms_uniprot_ac", "ms_unique_pep", "ms_intensity")
    
    # remove interactions with ms_intensity=0
    if (any(data_sel$ms_intensity == 0)) {
        data_sel <- data_sel[-which(data_sel$ms_intensity == 0), ]
    }
    
    # annotate proteins and add Masses for Mist
    species_split = unlist(strsplit(species, "-"))
    Uniprot = NULL
    for (org in species_split) {
        cat(sprintf("\tLOADING %s\n", org))
        tmp = read.delim2(sprintf("%s/uniprot_protein_descriptions_%s.txt", uniprot_dir, org), stringsAsFactors = F, quote = "")
        if (is.null(Uniprot)) {
            Uniprot = as.data.frame(tmp)
        } else {
            Uniprot = rbind(Uniprot, tmp)
        }
    }
    results_annotated = merge(data_sel, Uniprot, all.x = T, by.x = "ms_uniprot_ac", by.y = "Entry")
    
    # Adding the Mass column: it is required for MIST for preprocessing!  For that, we will calculate the mass of the whole protein
    # just taking the average molecular weight of an amino acid: 110Da
    results_annotated$Mass <- results_annotated$Length * 110
    
    write.table(results_annotated, file = output_file, eol = "\n", sep = "\t", quote = F, row.names = F, col.names = T)
    cat("\n>>MIST FILES CREATED!\n")
    
    keysout <- subset(keys, select = c(RawFile, Condition))
    keysname <- "keys_mistint.txt"
    write.table(keysout, file = keysname, eol = "\n", sep = "\t", quote = F, row.names = F, col.names = F)
    cat(">>MIST for INTENSITIES keys FILE ALSO CREATED. Truly enjoy your day ;-)\n\n")
}



# takes in a dcast data frame containing proteins (rows) and a value for each condition (columns)
MQutil.combine_sq_values <- function(dat, pos, neg) {
    # collapse POSTITIVE conditions together if there's more than one
    if (length(pos) > 1) {
        condition_names = names(dat)[grep(paste(paste("^", pos, "$", sep = ""), collapse = "|"), names(dat))]
        pos.counts <- dat[, c("Protein", condition_names)]
        # create a collapsed name for the new column
        collapsed_condition_name = paste(condition_names, collapse = "|")
        # Create new column of the collapsed counts
        pos.counts[, collapsed_condition_name] <- apply(pos.counts[, condition_names], 1, paste, collapse = "|")
        # keep only proteins and new column
        pos.counts <- pos.counts[, c("Protein", collapsed_condition_name)]
    } else {
        pos.counts <- dat[, c("Protein", pos)]
    }
    
    # collapse NEGATIVE conditions together if there's more than one
    if (length(neg) > 1) {
        condition_names = names(dat)[grep(paste(paste("^", neg, "$", sep = ""), collapse = "|"), names(dat))]
        neg.counts <- dat[, c("Protein", condition_names)]
        # create a collapsed name for the new column
        collapsed_condition_name = paste(condition_names, collapse = "|")
        # Create new column of the collapsed counts
        neg.counts[, collapsed_condition_name] <- apply(neg.counts[, condition_names], 1, paste, collapse = "|")
        # keep only proteins and new column
        neg.counts <- neg.counts[, c("Protein", collapsed_condition_name)]
    } else {
        neg.counts <- dat[, c("Protein", neg)]
    }
    
    # combine pos and neg values together into final representation
    all.counts <- merge(pos.counts, neg.counts, by = "Protein")
    final_name = paste(names(pos.counts)[2], names(neg.counts)[2], sep = "-")
    all.counts[, final_name] <- apply(all.counts[, -1], 1, paste, collapse = " - ")
    all.counts <- all.counts[, c("Protein", final_name)]
    
    return(all.counts)
}


#' @title Add Abundancd data to Results.
#' @description Aggregates the normalized abundance and replicate data from the samples. Uses the MSstat output file ...mss-sampleQuant.txt for the aggregations, and is applied directly to the MSstats results in wide format. The resulting file will have 'abundance' appended to the end of the file name.
#' @param sq_file The filepath to the Sample Quantification file (txt tab delimited file).
#' @param contrast_file The filepath to the Contrst file used to generate the MSstats results (txt tab delimited file).
#' @param results_file The filepath to the results in wide format file used (txt tab delimited file).
#' @keywords abundance msstats sample quatification samplequantification
#' MQutil.sampleQuant()
# Main wrapper that consolidates the abundance data for a results.wide file
MQutil.sampleQuant <- function(sq_file, contrast_file, results_file) {
    cat(">>SUMMARIZING ABUNDANCE DATA\n")
    
    # Read in sampleQuant data
    cat(">> LOADING DATA FILE\n")
    x <- read.delim(sq_file, sep = "\t", stringsAsFactors = F)
    
    # read in the results-wide or results-wide-annotated data
    results.wide = read.delim(results_file, stringsAsFactors = F)
    
    # identify protein column name. This can be different depending on if it's PTM data or not
    if (length(grep("mod_sites", names(results.wide))) > 0) {
        prot_col = "mod_sites"
    } else {
        prot_col = "Protein"
    }
    
    # CONTRASTS
    cat(">>  LOADING CONTRAST FILE\n")
    contrasts = read.delim(contrast_file, stringsAsFactors = F)
    # make sure the column names are in alphabetical order before continuing
    contrasts = as.matrix(contrasts[, order(dimnames(contrasts)[[2]], decreasing = F)])
    
    # Begin work ~~~~~~~~~~~~ put data into long format, remove missing (NA) values
    x <- melt(x, factorsAsStrings = F, na.rm = T)
    
    # Keep only condiditon names
    x$variable <- unlist(lapply(strsplit(as.character(x$variable), "_"), function(y) {
        l = length(y)
        return(paste(y[-l], collapse = "_"))
    }))
    
    # count how many times a protein shows up for a condition
    cat("\tSUMMARIZING REPLICATE COUNTS\n")
    x.counts <- dcast(x, Protein ~ variable, fill = NA_real_)
    
    # Combine all the normalized intensities for each replicate
    cat("\tSUMMARIZING REPLICATE INTENSITIES\n")
    x$value <- round(x$value, 2)
    x.intensities <- dcast(Protein ~ variable, data = x, paste, collapse = ";", fill = NA_character_)
    
    ########################################### ~~~~~~~~~~ GET CONTRAST DATA ~~~~~~~~~~~~~
    cat("\tAGGREGATING CONTRASTS\n")
    x.counts.list <- x.intensities.list <- list()
    for (i in 1:dim(contrasts)[1]) {
        # get which conditions are being contrasted this time
        conditions <- contrasts[i, grep("[^0]", contrasts[i, ])]
        # order the conditions so the positives are on the left, negatives on the right
        conditions <- conditions[order(conditions, decreasing = T)]
        
        # get the names associated with each side of the contrast
        pos <- names(conditions)[which(conditions > 0)]
        neg <- names(conditions)[which(conditions < 0)]
        
        # combine all the counts into a single summarized column
        x.counts.list[[i]] <- MQutil.combine_sq_values(dat = x.counts, pos = pos, neg = neg)
        # add an extra space to make excel not convert everything to dates
        x.counts.list[[i]][, 2] = paste0(" ", x.counts.list[[i]][, 2])
        # combine all the intensities into a single summarized column
        x.intensities.list[[i]] <- MQutil.combine_sq_values(dat = x.intensities, pos = pos, neg = neg)
    }
    
    # combine all the goods together
    x.counts = Reduce(function(...) merge(..., all = T), x.counts.list)
    x.intensities = Reduce(function(...) merge(..., all = T), x.intensities.list)
    
    names(x.counts)[-1] = paste0(names(x.counts)[-1], "_count")
    names(x.intensities)[-1] = paste0(names(x.intensities)[-1], "_intensity")
    # combing all the results together
    results.all <- merge(x.counts, x.intensities, by = "Protein")
    
    # add these to the original results-wide file
    results.all <- merge(results.wide, results.all, by.x = prot_col, by.y = "Protein")
    # write out the file
    out_file <- gsub(".txt", "-abundance.txt", results_file)
    write.table(results.all, out_file, quote = F, row.names = F, sep = "\t")
    
    cat(">> ABUNDANCE SUMMARIZATION COMPLETE. HAVE A NICE DAY :)\n")
    # return(results.all)
}



# create the replicate plots based on the pairings from the replicate plot file.
MQutil.replicatePlots <- function(input_file, keys_file, replicate_file, out_file) {
    cat(">> READIN IN FILES...\n")
    # read in data
    dat <- read.delim(input_file, stringsAsFactors = F)
    # keys
    keys <- read.delim(keys_file, stringsAsFactors = F)
    # profile plot list
    repplot <- read.delim(replicate_file, stringsAsFactors = F)
    
    # remove negatives from MaxQuant
    if (length(grep("__", dat$Proteins)) > 0) 
        dat <- dat[-grep("__", dat$Proteins), ]
    
    # for UB, jj suggests using unique peptide and charge for distinguishing numbers NOTE: dimensions betwen x and dat may differ
    # if there is data in dat that isn't in the keys file
    names(dat)[grep("Raw.file", names(dat))] <- "RawFile"
    # x <- merge(dat, keys[,c('RawFile','Condition','BioReplicate','IsotopeLabelType')], by=c('RawFile', 'IsotopeLabelType') ) ##
    # !!!!!!! DIfferent for SILAC
    x <- merge(dat, keys[, c("RawFile", "Condition", "BioReplicate")], by = c("RawFile"))
    
    # Put into a data matrix format
    x <- dcast(data = x, Proteins + Modified.sequence + Charge ~ Condition + BioReplicate, value.var = "Intensity", max, na.rm = T)
    # remove cases where -Inf is introduced
    x[x == -Inf] = 0  ###### May cause problems? Check.
    write.table(x, out_file, quote = F, row.names = F, sep = "\t")
    
    # cycle through the condition pairs in the file and plot each pair
    for (i in 1:dim(repplot)[1]) {
        cat(">>  PLOTTING REPLICATE PLOT ", i, "\n")
        
        # check if the replicate combination exists in the plots
        rep1_1 <- paste(repplot$condition1[i], repplot$rep1_1[i], sep = "_")
        rep1_2 <- paste(repplot$condition1[i], repplot$rep1_2[i], sep = "_")
        rep2_1 <- paste(repplot$condition2[i], repplot$rep2_1[i], sep = "_")
        rep2_2 <- paste(repplot$condition2[i], repplot$rep2_2[i], sep = "_")
        reps <- c(rep1_1, rep1_2, rep2_1, rep2_2)
        if (!any(!(reps %in% names(x)))) {
            
            # prep 1st replicate comparison for plot
            rep1 <- log2(x[, paste(repplot$condition1[i], repplot$rep1_1[i], sep = "_")]/x[, paste(repplot$condition2[i], repplot$rep2_1[i], 
                sep = "_")])
            # prep 2nd replicate comparison for plot
            rep2 <- log2(x[, paste(repplot$condition1[i], repplot$rep1_2[i], sep = "_")]/x[, paste(repplot$condition2[i], repplot$rep2_2[i], 
                sep = "_")])
            
            # remove NA pairs
            idx <- which(is.na(rep1) | is.na(rep2) | is.infinite(rep1) | is.infinite(rep2))
            rep1 <- rep1[-idx]
            rep2 <- rep2[-idx]
            reps.cor <- cor(rep1, rep2, use = "pairwise.complete.obs", method = "pearson")
            
            # set up a square plot centered at 0
            x.lim <- ceiling(max(abs(c(rep1, rep2)), na.rm = T))
            y.lim <- c(-x.lim, x.lim)
            x.lim <- c(-x.lim, x.lim)
            
            # name axes labels
            y.label = paste0(repplot$condition1[i], " vs. ", repplot$condition2[i], "  (", repplot$rep1_1[i], "/", repplot$rep2_1[i], 
                ")")
            x.label = paste0(repplot$condition1[i], " vs. ", repplot$condition2[i], "  (", repplot$rep1_2[i], "/", repplot$rep2_2[i], 
                ")")
            # make plot name
            plot.name <- paste(repplot$condition1[i], " vs. ", repplot$condition2[i], " : R = ", round(reps.cor, 3), sep = "")
            # pdf( paste( dirname(out_file), '/', gsub(' ','_',plot.name) ,'_', repplot$rep1_1[i], '_', repplot$rep1_2[i], '.pdf', sep='')
            # ) plot(rep1, rep2, main=plot.name, xlab=repplot$rep1_1[i], ylab=repplot$rep1_2[i], xlim=x.lim, ylim=y.lim, pch='.') dev.off()
            tmp <- data.frame(rep1, rep2, stringsAsFactors = F)
            p <- ggplot(tmp, aes(x = rep1, y = rep2)) + geom_point() + xlim(x.lim[1], x.lim[2]) + ylim(x.lim[1], x.lim[2]) + ggtitle(plot.name) + 
                labs(x = x.label, y = y.label)
            ggsave(filename = paste(dirname(out_file), "/", gsub(" ", "_", plot.name), "_", repplot$rep1_1[i], "_", repplot$rep1_2[i], 
                ".pdf", sep = ""), plot = p, width = 10, height = 10)
            
            
        } else {
            warning("REPLICATE PLOT ", i, " NOT MADE -- MISSING DATA FROM \n", paste("\t", reps[!(reps %in% names(x))], "\n", collapse = ""))
        }
    }
    
}


# main <- function(opt){ if(opt$command %in% ALLOWED_COMMANDS){ cat(sprintf('>> EXECUTING:\t%s\n',opt$command)) loadLibs()
# if(opt$command == 'concat'){ MQutil.concat(filenames=opt$files, output = opt$output) }else if(opt$command ==
# 'convert-silac'){ MQutil.SILACToLong(filename = opt$files, output = opt$output) }else if(opt$command == 'keys'){
# MQutil.getKeys(filename = opt$files, output = opt$output) }else if(opt$command == 'convert-sites'){
# MQutil.ProteinToSiteConversion (maxq_file = opt$files, output_file = opt$output, ref_proteome_file = opt$proteome, mod_type =
# opt$mod_type) }else if(opt$command == 'annotate'){ MQutil.annotate(input_file = opt$files, output_file = opt$output,
# species=opt$species, uniprot_dir=opt$uniprot_dir) }else if(opt$command == 'results-wide'){ MQutil.resultsWide(input_file =
# opt$files, output_file = opt$output ) }else if(opt$command == 'mapback-sites'){ MQutil.mapSitesBack(input_file = opt$files,
# output_file = opt$output , mapping_file = opt$mapping) }else if(opt$command == 'heatmap'){ MQutil.plotHeatmap(input_file =
# opt$files, output_file = opt$output, labels = opt$labels, lfc_lower = opt$lfc_lower, lfc_upper=opt$lfc_upper, FDR =
# opt$q_value) }else if(opt$command == 'simplify'){ MQutil.simplify(input_file = opt$files, output_file = opt$output) }else
# if(opt$command == 'saint-format'){ MQutil.MaxQToSaint(data_file = opt$files, keys_file = opt$keys, ref_proteome_file =
# opt$proteome, quant_variable = opt$identifier_column) }else if(opt$command == 'data-plots'){ MQutil.dataPlots(input_file =
# opt$files, output_file = opt$output) }else if(opt$command == 'spectral-counts'){ MQutil.spectralCounts(input_file =
# opt$files, keys_file = opt$keys, output_file = opt$output) }else if(opt$command == 'mist'){ MQutil.MISTformat(input_file =
# opt$files, keys_file = opt$keys, output_file = opt$output, species=opt$species, uniprot_dir=opt$uniprot_dir) }else
# if(opt$command == 'mistint'){ MQutil.MISTINTformat(input_file = opt$files, keys_file = opt$keys, output_file = opt$output,
# species=opt$species, uniprot_dir=opt$uniprot_dir) }else if(opt$command == 'samplequant'){ MQutil.sampleQuant(sq_file =
# opt$files, contrast_file = opt$contrast_file, results_file = opt$results_file) }else if(opt$command == 'replicateplots'){
# MQutil.replicatePlots(input_file = opt$files, keys_file = opt$keys, replicate_file = opt$results_file, out_file = opt$output)
# } }else{ cat(sprintf('COMMAND NOT ALLOWED:\t%s\n',opt$command)) cat(sprintf('ALLOWED
# COMMANDS:\t%s\n',paste(ALLOWED_COMMANDS,collapse=','))) } } # opt$command = 'data-plots' # opt$files =
# '~/Projects/HPCKrogan/Data/TBLMSE/data/swissprot/TBLMSE-cox-ub-swissprot-modK-mss-normalized.txt' # opt$output =
# '~/Projects/HPCKrogan/Data/TBLMSE/data/swissprot/TBLMSE-cox-ub-swissprot-modK-mss-normalized.pdf' # opt$command =
# 'saint-format' # opt$files =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-human-exvivo/AP-MS/M2/data/FLU-HUMAN-ALL-APMS-M2-data.txt' #
# opt$keys =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-human-exvivo/AP-MS/M2/data/FLU-HUMAN-ALL-APMS-M2-keys.txt' #
# opt$proteome = '~/Projects/HPCKrogan/Data/Uniprot/homo-sapiens-swissprot.fasta' # opt$command = 'heatmap' # opt$input =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Flu-human-exvivo/H5N1/ph//results//20141203//FLU-HUMAN-H5N1-PH-results.txt' #
# opt$output =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Flu-human-exvivo/H5N1/ph//results//20141203//FLU-HUMAN-H5N1-PH-results.pdf' #
# opt$labels=NULL, names='Protein', cluster_cols=F, display='log2FC' # opt$command = 'results-wide' # opt$input =
# '~/Projects/HPCKrogan/Data/HIV-proteomics/results/20141124-ub-proteins/HIV-UB-SILAC-KROGAN-results.txt' # opt$output =
# '~/Projects/HPCKrogan/Data/HIV-proteomics/results/20141124-ub-proteins/HIV-UB-SILAC-KROGAN-results.txt' # opt$command =
# 'annotate' # opt$input =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-human-exvivo/H5N1/ph//results//20141214//FLU-HUMAN-H5N1-PH-results-mapped.txt'
# # opt$output =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-human-exvivo/H5N1/ph//results//20141214//FLU-HUMAN-H5N1-PH-results-mapped-ann.txt'
# # opt$command = 'convert-silac' # opt$files = '~/Code/RMSQ/tests/abundance/ENTERO-LG-ABU-data.txt' # opt$output =
# '~/Code/RMSQ/tests/abundance/ENTERO-LG-ABU-data-long.txt' # opt$command = 'concat' # opt$files =
# '~/Projects/HPCKrogan/Data/Mtb/Files/073113*evidence.txt' # opt$output =
# '~/Projects/HPCKrogan/Data/HIV-proteomics/Meena/abundance/HIV_vs_MOCK_PROTEIN_evidence_split.txt' # opt$command =
# 'convert-sites' # opt$files = '~/Projects/HPCKrogan/Data/HIV-proteomics/data//UB-silac//HIV-UB-SILAC-KROGAN-data.txt' #
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/data//UB-silac//HIV-UB-SILAC-KROGAN-data-modK.txt' # opt$mod_type =
# 'ub' # opt$proteome = '~/Projects/HPCKrogan/Data/Uniprot/homo-sapiens-swissprot.fasta' # opt$command = 'convert-sites' #
# opt$files =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/Pilot/Mouse-H5N1-Ph//data//FLU-MOUSE-H5N1-PH-PILOT-data.txt'
# # opt$output =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/Pilot/Mouse-H5N1-Ph//data//FLU-MOUSE-H5N1-PH-PILOT-data-modSTY.txt'
# # opt$mod_type = 'ph' # opt$proteome = '~/Projects/HPCKrogan/Data/Uniprot/mus-musculus-uniprot.fasta' # opt$command =
# 'mapback-sites' # opt$files =
# '~/Projects/HPCKrogan/Data/TBLMSE/results/20141124-sites-eqmed-noint/TBLMSE-cox-ub-swissprot-modK-results.txt' # opt$mapping
# = '~/Projects/HPCKrogan/Data/TBLMSE/data/swissprot/TBLMSE-cox-ub-swissprot-modK-mapping.txt' # opt$output =
# '~/Projects/HPCKrogan/Data/TBLMSE/results/20141124-sites-eqmed-noint/TBLMSE-cox-ub-swissprot-modK-results.txt' # opt$command
# = 'heatmap' # opt$files =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/results/20150114/FLU-MOUSE-H1N1-UB-results-ann.txt'
# # opt$output =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/results/20150114/FLU-MOUSE-H1N1-UB-results-ann.pdf'
# # opt$command = 'simplify' # opt$files =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/results/20150114-sites/FLU-MOUSE-H1N1-UB-results.txt'
# # opt$output =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/results/20150114-sites/FLU-MOUSE-H1N1-UB-results-simplified.txt'
# # opt$command = 'heatmap' # opt$files = '~/Code/RMSQ/tests/heatmap/031615-lm-1-6-ph-mss-results.txt' # opt$output =
# '~/Code/RMSQ/tests/heatmap/031615-lm-1-6-ph-mss-results.pdf' # opt$labels = '*' # opt$lfc_lower = -2 # opt$lfc_upper = 2 #
# opt$q_value = 0.05 # opt$command = 'spectral-counts' # opt$files =
# '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/data//FLU-MOUSE-H1N1-UB-data.txt' # opt$keys
# = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/data/FLU-MOUSE-H1N1-UB-keys.txt' #
# opt$output = '~/Desktop/test.txt' # opt=c() # opt$command = 'samplequant' # opt$files = '~/Box
# Sync/projects/CoxJeff/results/20160801_MSStat_results_UB/msstats/TBLMSE-ub-SE-results-mss-sampleQuant.txt' # opt$results_file
# = '~/Box Sync/projects/CoxJeff/results/20160801_MSStat_results_UB/msstats/TBLMSE-ub-SE-results-wide-sites-ann.txt' #
# opt$contrast_file = '~/Box Sync/projects/CoxJeff/data/jj/ub/TBLMSE-cox-ub-silac-SE-contrasts.txt' # # # opt$files = '~/Box
# Sync/projects/FluomicsProteomics/other/Marazzi/results/MSstats3-results-mss-sampleQuant.txt' # opt$contrast_file = '~/Box
# Sync/projects/FluomicsProteomics/other/Marazzi/data/contrasts.txt' # opt$results_file = '~/Box
# Sync/projects/FluomicsProteomics/other/Marazzi/results/MSstats3-results-wide-annotated.txt' main(opt)
