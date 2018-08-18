





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
