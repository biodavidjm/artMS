#' @import RColorBrewer
#' @import pheatmap
#' @import ggplot2

is.uniprotAc <- function(identifier) {
    grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", identifier)
}

filterMaxqData <- function(data) {
    data_selected <- data[grep("CON__|REV__", data$Proteins, invert = T), ]
    idx <- which(data_selected$Proteins == "")
    if (length(idx) > 0) 
        data_selected <- data_selected[-idx, ]
    return(data_selected)
}

explodeMaxQProteinGroups <- function(data) {
    return(data)
}

removeMaxQProteinGroups <- function(data) {
    data_selected <- data[grep(";", data$Proteins, invert = T), ]
    return(data_selected)
}

#' MELTING AND CASTING

naMax <- function(x) {
    return(max(x, na.rm = T))
}

castMaxQToWide <- function(d_long) {
    data_w <- dcast.data.table(Proteins + Sequence + Charge ~ RawFile + IsotopeLabelType, data = d_long, value.var = "Intensity", 
        fun.aggregate = sum, fill = NA)
    # setnames(data_w,2,'Sequence')
    return(data_w)
}

castMaxQToWidePTM <- function(d_long) {
    data_w <- dcast.data.table(Proteins + Modified.sequence + Charge ~ RawFile + IsotopeLabelType, data = d_long, value.var = "Intensity", 
        fun.aggregate = sum, fill = NA)
    setnames(data_w, 2, "Sequence")
    return(data_w)
}

#' eg. getIsotopeLabel('VE20130426_04_dbdbk_Light')
getIsotopeLabel <- function(str) {
    last_idx <- regexpr("\\_[^\\_]*$", str)[1]
    return(substr(str, last_idx + 1, nchar(str)))
}

## eg. getFileNameWithoutLabel('VE20130426_04_dbdbk_Light')
getFileNameWithoutLabel <- function(str) {
    last_idx <- regexpr("\\_[^\\_]*$", str)[1]
    return(substr(str, 0, last_idx - 1))
}

meltMaxQToLong <- function(data_w, na.rm = F) {
    data_l <- reshape2::melt(data_w, id.vars = c("Proteins", "Sequence", "Charge"), na.rm = na.rm)
    setnames(data_l, old = 4:5, new = c("RawFile", "Intensity"))
    data_l[, `:=`("IsotopeLabelType", NA)]
    data_l$IsotopeLabelType <- apply(data_l, 1, function(x) getIsotopeLabel(x["RawFile"]))
    data_l$RawFile <- apply(data_l, 1, function(x) getFileNameWithoutLabel(x["RawFile"]))
    return(data_l)
}

fillMissingMaxQEntries <- function(data_w, perRun = F) {
    mins <- apply(data_w[, 4:ncol(data_w), with = F], 2, function(x) min(x[x > 0], na.rm = T))
    for (j in (4:ncol(data_w))) if (perRun) {
        set(data_w, which(is.na(data_w[[j]])), j, mins[j - 3])
    } else {
        set(data_w, which(is.na(data_w[[j]])), j, min(mins))
    }
    return(data_w)
}

cleanMissingMaxQEntries <- function(data_l) {
    data_l[data_l$Intensity <= 0 | is.infinite(data_l$Intensity) | is.nan(data_l$Intensity), ]$Intensity = NA
    return(data_l)
}

flattenMaxQTechRepeats <- function(data_l) {
    
}

flattenKeysTechRepeats <- function(keys) {
    keys_agg <- aggregate(. ~ BioReplicate, data = keys, FUN = function(x) unique(x)[1])
    keys_agg$Run <- keys_agg$BioReplicate
    return(keys_agg)
}


#' @title Merge MaxQuant Data with Keys
#' @description This will merge a MaxQuant evidence file and keys file.
#' @param data Data from the evidence file in a data.table .
#' @param keys The information from the keys file in a data.table .
#' @param by What column(s) to merge the data and keys data on. Default = 'RawFile'
#' @keywords MSStats MaxQuant mss
#' mergeMaxQDataWithKeys()
#' @export
mergeMaxQDataWithKeys <- function(data, keys, by = c("RawFile")) {
    # Check if the number of RawFiles is the same.
    unique_data <- unique(data$RawFile)
    unique_keys <- unique(keys$RawFile)
    
    if (length(unique_keys) != length(unique_data)) {
        keys_not_found <- setdiff(unique_keys, unique_data)
        data_not_found <- setdiff(unique_data, unique_keys)
        cat(sprintf("keys found: %s \t keys not in data file:\n%s\n", length(unique_keys) - length(keys_not_found), paste(keys_not_found, 
            collapse = "\t")))
        cat(sprintf("data found: %s \t data not in keys file:\n%s\n", length(unique_data) - length(data_not_found), paste(data_not_found, 
            collapse = "\t")))
    } else {
        cat("\nCheck point: the number of RawFiles in both keys and evidences file is identical\n\n")
    }
    
    ## select only required attributes from MQ format
    data <- merge(data, keys, by = by)
    return(data)
}

normalizePerCondition <- function(d, NORMALIZATION_METHOD = "scale") {
    unique_conditions <- unique(d$Condition)
    d_tmp <- c()
    
    for (u in unique_conditions) {
        print(sprintf("normalizing\t%s", u))
        ss <- data_w[d$Condition == u, ]
        tmp <- normalizeSingle(ss, NORMALIZATION_METHOD)
        d_tmp <- rbind(d_tmp, tmp)
    }
    d_tmp
}

na.replace <- function(v, value = 0) {
    v[is.na(v)] = value
    return(v)
}

myNormalizeMedianValues <- function(x) {
    narrays <- NCOL(x)
    if (narrays == 1) 
        return(x)
    cmed <- log(apply(x, 2, median, na.rm = TRUE))
    cmed <- exp(cmed - mean(cmed, na.rm = T))
    t(t(x)/cmed)
}

normalizeSingle <- function(data_w, NORMALIZATION_METHOD = "scale") {
    
    data_part <- as.matrix(data_w[, 4:ncol(data_w), with = F])
    if (NORMALIZATION_METHOD == "scale") {
        res <- myNormalizeMedianValues(data_part)
    } else if (NORMALIZATION_METHOD == "quantile") {
        res <- normalizeBetweenArrays(data_part, method = NORMALIZATION_METHOD)
    }
    data_part_n <- normalizeBetweenArrays(data_part, method = NORMALIZATION_METHOD)
    res_f <- data.table(cbind(data_w[, 1:3, with = F], res))
    return(res_f)
}

dataToMSSFormat <- function(d) {
    tmp <- data.frame(ProteinName = d$Proteins, PeptideSequence = d$Sequence, PrecursorCharge = NA, FragmentIon = NA, ProductCharge = d$Charge, 
        IsotopeLabelType = d$IsotopeLabelType, Condition = d$Condition, BioReplicate = d$BioReplicate, Run = d$Run, Intensity = d$Intensity)
    tmp
}

samplePeptideBarplot <- function(data_f, config) {
    # set up data into ggplot compatible format
    data_f <- data.table(data_f, labels = paste(data_f$RawFile, data_f$Condition, data_f$BioReplicate))
    data_f <- data_f[with(data_f, order(labels, decreasing = T)), ]
    
    # plot the peptide counts for all the samples TOGETHER
    p <- ggplot(data = data_f, aes(x = labels))
    p <- p + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1, family = "mono")) + ggtitle("Unique peptides per run\n after filtering") + 
        coord_flip()
    ggsave(filename = gsub(".txt", "-peptidecounts.pdf", config$files$output), plot = p, width = 8, height = 10)
    
    # plot the peptide counts for all the samples PER BAIT
    p <- ggplot(data = data_f, aes(x = BioReplicate))
    p <- p + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1, family = "mono")) + ggtitle("Unique peptides per run\n after filtering") + 
        facet_wrap(~Condition, scales = "free") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(filename = gsub(".txt", "-peptidecounts-perBait.pdf", config$files$output), plot = p, width = 8, height = 10)
    
}

sampleCorrelationHeatmap <- function(data_w, keys, config) {
    mat = log2(data_w[, 4:ncol(data_w), with = F])
    mat[is.na(mat)] = 0
    mat_cor = cor(mat, method = "pearson", use = "everything")
    ordered_keys = keys[with(keys, order(RawFile)), ]  ## we want to make informarive row names so order by RawFile because that's how data_w is ordered
    mat_names = paste(ordered_keys$Condition, ordered_keys$BioReplicate, ordered_keys$Run)
    colnames(mat_cor) = mat_names
    rownames(mat_cor) = mat_names
    colors_pos = colorRampPalette(brewer.pal("Blues", n = 5))(10)
    colors_neg = rev(colorRampPalette(brewer.pal("Reds", n = 5))(10))
    colors_tot = c(colors_neg, colors_pos)
    pheatmap(mat = mat_cor, cellwidth = 10, cellheight = 10, scale = "none", filename = gsub(".txt", "-heatmap.pdf", config$files$output), 
        color = colors_tot, breaks = seq(from = -1, to = 1, by = 0.1), fontfamily = "mono")
    
}

plotHeat <- function(mss_F, out_file, labelOrder = NULL, names = "Protein", cluster_cols = F, display = "log2FC") {
    heat_data = data.frame(mss_F, names = names)
    # heat_data = mss_F[,c('uniprot_id','Label','log2FC')]
    
    ## create matrix from log2FC or p-value as user defined
    if (display == "log2FC") {
        # Issues with extreme_val later if we have Inf/-Inf values.
        if (sum(is.infinite(heat_data$log2FC)) > 0) {
            idx <- is.infinite(heat_data$log2FC)
            heat_data$log2FC[idx] <- NA
        }
        heat_data_w = dcast(names ~ Label, data = heat_data, value.var = "log2FC")
    } else if (display == "pvalue") {
        heat_data$adj.pvalue = -log10(heat_data$adj.pvalue + 10^-16)
        heat_data_w = dcast(names ~ Label, data = heat_data, value.var = "adj.pvalue")
    }
    
    ## try gene_names = uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
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
        cat("\t Saving heatmap\n")
        pheatmap(heat_data_w, scale = "none", cellheight = 10, cellwidth = 10, filename = out_file, color = colors_tot, breaks = seq(from = -extreme_val, 
            to = extreme_val, by = bin_size), cluster_cols = cluster_cols, fontfamily = "mono")
    } else {
        heat_data_w = heat_data_w[, labelOrder]
        pheatmap(heat_data_w, scale = "none", cellheight = 10, cellwidth = 10, filename = out_file, color = colors_tot, breaks = seq(from = -extreme_val, 
            to = extreme_val, by = bin_size), cluster_cols = cluster_cols, fontfamily = "mono")
    }
    
    heat_data_w
}

significantHits <- function(mss_results, labels = "*", LFC = c(-2, 2), FDR = 0.05) {
    ## get subset based on labels
    selected_results = mss_results[grep(labels, mss_results$Label), ]
    cat(sprintf("\tAVAILABLE LABELS FOR HEATMAP:\t%s\n", paste(unique(mss_results$Label), collapse = ",")))
    cat(sprintf("\tSELECTED LABELS FOR HEATMAP:\t%s\n", paste(unique(selected_results$Label), collapse = ",")))
    significant_proteins = selected_results[(!is.na(selected_results$log2FC) & selected_results$adj.pvalue <= FDR & (selected_results$log2FC >= 
        LFC[2] | selected_results$log2FC <= LFC[1])), "Protein"]
    significant_results = selected_results[selected_results$Protein %in% significant_proteins, ]
    return(significant_results)
}

logScale <- function(data, format = "wide", base = 2) {
    if (format == "wide") {
        data[, 4:ncol(data)] = log(data[, 4:ncol(data)], base = base)
    } else {
        data$Intensity = log(data$Intensity, base = base)
    }
    return(data)
}

################################### doesnt work so far with new code

peptideDistribution <- function(data_l, output_file, PDF = T) {
    ## look at peptide distribution of all proteins
    if (PDF) 
        pdf(output_file, width = 10, height = 7)
    p = ggplot(data = data_l, aes(x = Raw.file, y = Intensity))
    print(p + geom_boxplot() + theme(axis.text.x = element_text(angle = -90)))
    if (PDF) 
        dev.off()
}

peptideIntensityPerFile <- function(ref_peptides, output_file, PDF = T) {
    if (PDF) 
        pdf(output_file, width = 10, height = 7)
    p = ggplot(data = ref_peptides, aes(x = Raw.file, y = Intensity, group = protein_id))
    print(p + geom_line(colour = "grey") + stat_summary(aes(group = 1), geom = "point", fun.y = median, shape = 17, size = 3, colour = "darkred") + 
        ylab("Intensity") + theme(axis.text.x = element_text(angle = -90)))
    if (PDF) 
        dev.off()
}



prettyPrintHeatmapLabels <- function(uniprot_acs, uniprot_ids, gene_names) {
    # uniprot_ids_trunc = gsub('([A-Z,0-9]+)_([A-Z,0-9]+)','\\1',uniprot_ids) longest_id = max(nchar(uniprot_ids_trunc))
    # tmp_frame = data.frame(t=uniprot_ids_trunc, s=longest_id-nchar(uniprot_ids_trunc)+1, g=gene_names,a=uniprot_acs,
    # stringsAsFactors=F) tmp_frame[is.na(tmp_frame$t),]$t=tmp_frame[is.na(tmp_frame$t),]$a result = apply(tmp_frame, 1,
    # function(x)paste0(x[1],paste(rep(' ',x[2]),collapse=''),x[3])) result = apply(tmp_frame, 1, function(x)paste0(x[4],'
    # ',x[3],collapse=''))
    result = paste(uniprot_acs, uniprot_ids, gene_names, sep = " ")
    return(result)
}

normalizeToReference <- function(data_l_ref, ref_protein, PDF = T, output_file) {
    
    data_l_ref$Intensity = log2(data_l_ref$Intensity)
    
    ## SAMPLES WITH REF PROTEINS BEFORE
    peptideDistribution(data_l_ref, gsub(".txt", "-protein-dist-b.pdf", output_file))
    
    ## compute the average value for background reference peptides over all samples (background+bait) normalize against background
    ## reference proteins
    ref_peptides = data_l_ref[grep(ref_protein, data_l_ref$Proteins), ]
    ref_peptides = data.frame(ref_peptides, protein_id = as.character(paste0(ref_peptides$Sequence, ref_peptides$Charge)), stringsAsFactors = F)
    ref_peptides_counts = aggregate(Raw.file ~ protein_id, data = ref_peptides, FUN = function(x) length(unique(x)))
    unique_files = length(unique(data_l_ref$Raw.file))
    complete_observations = ref_peptides_counts[ref_peptides_counts$Raw.file == unique_files, "protein_id"]
    ref_peptides_complete = ref_peptides[ref_peptides$protein_id %in% complete_observations, ]
    
    ## look at their individual signal over replicates
    peptideIntensityPerFile(ref_peptides_complete[!is.na(ref_peptides_complete$Intensity) & is.finite(ref_peptides_complete$Intensity), 
        ], gsub(".txt", "-peptide-signal-b.pdf", output_file))
    
    ref_peptides_avg_per_rep = aggregate(Intensity ~ Raw.file, data = ref_peptides_complete, FUN = median)
    ref_peptides_avg_all = median(ref_peptides_avg_per_rep$Intensity)
    ref_peptides_avg_per_rep = data.frame(Raw.file = ref_peptides_avg_per_rep$Raw.file, correction = ref_peptides_avg_all - ref_peptides_avg_per_rep$Intensity)
    ref_peptides_complete = merge(ref_peptides_complete, ref_peptides_avg_per_rep, by = "Raw.file", all.x = T)
    ref_peptides_complete = data.frame(ref_peptides_complete, IntensityCorrected = ref_peptides_complete$Intensity + ref_peptides_complete$correction)
    
    ref_peptides_complete$Intensity = ref_peptides_complete$IntensityCorrected
    peptideIntensityPerFile(ref_peptides_complete[!is.na(ref_peptides_complete$Intensity) & is.finite(ref_peptides_complete$Intensity), 
        ], gsub(".txt", "-peptide-signal-a.pdf", output_file))
    
    data_l_ref = merge(data_l_ref, ref_peptides_avg_per_rep, by = "Raw.file", all.x = T)
    data_l_ref$Intensity = data_l_ref$Intensity + data_l_ref$correction
    data_l_ref = data_l_ref[, -(ncol(data_l_ref))]
    
    ## SAMPLES WITH REF PROTEINS AFTER
    peptideDistribution(data_l_ref, gsub(".txt", "-protein-dist-a.pdf", output_file))
    data_l_ref$Intensity = 2^(data_l_ref$Intensity)
    
    return(data_l_ref)
}

simplifyAggregate <- function(str, sep = ",", numeric = F) {
    str_vec = unlist(str_split(str, pattern = sep))
    if (numeric) {
        str_vec = sort(unique(as.numeric(str_vec)))
    } else {
        str_vec = sort(unique(str_vec))
    }
    
    str_new = str_join(as.character(str_vec), collapse = ";")
    return(str_new)
}

simplifyOutput <- function(input) {
    input$Protein = apply(input, 1, function(x) simplifyAggregate(unname(x["Protein"]), sep = ";"))
    if (any(grepl("mod_sites", colnames(input)))) {
        input$mod_sites = apply(input, 1, function(x) simplifyAggregate(unname(x["mod_sites"])))
    }
    if (any(grepl("uniprot_ac", colnames(input)))) {
        input$uniprot_ac = apply(input, 1, function(x) simplifyAggregate(unname(x["uniprot_ac"])))
    }
    if (any(grepl("entrezgene", colnames(input)))) {
        input$entrezgene = apply(input, 1, function(x) simplifyAggregate(unname(x["entrezgene"]), numeric = T))
    }
    if (any(grepl("uniprot_genename", colnames(input)))) {
        input$uniprot_genename = apply(input, 1, function(x) simplifyAggregate(unname(x["uniprot_genename"])))
    }
    if (any(grepl("description", colnames(input)))) {
        input$description = apply(input, 1, function(x) simplifyAggregate(unname(x["description"])))
    }
    input[, `:=`(sample_1, gsub("([A-Z,0-9,a-z,_,\\s]+)\\-{1}([A-Z,0-9,a-z,_,\\s]+)", "\\1", input$Label))]
    input[, `:=`(sample_2, gsub("([A-Z,0-9,a-z,_,\\s]+)\\-{1}([A-Z,0-9,a-z,_,\\s]+)", "\\2", input$Label))]
    input[, `:=`(Label, NULL)]
    if (any(grepl("group_id", colnames(input)))) {
        input[, `:=`(group_id, NULL)]
    }
    return(input)
}
