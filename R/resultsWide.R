#' @import reshape2


#' @title Convert default MSStat results file to Wide format.
#' @description Converts the normal MSStats output file into 'wide' format where each row represents a protein's results, and each column represents the comparison made by MSStats. The fold change and p-value of each comparison will be it's own column.
#' @param input_file The filepath to the MaxQuant searched results file (txt tab delimited file). This can either be in either long or wide format.
#' @param save_file The filepath to the intended output file (txt tab delimited file).
#' @param return_results Whether to return the results at the end of the function (default =F).
#' @keywords wide MSStats results
#' resultsWide()
resultsWide = function(input_file, save_file = F, return_results = F) {
    input = read.delim(input_file, stringsAsFactors = F)
    input_l = melt(data = input[, c("Protein", "Label", "log2FC", "adj.pvalue")], id.vars = c("Protein", "Label"))
    
    ## then cast to get combinations of LFCV/PVAl and Label as columns
    input_w = dcast(Protein ~ Label + variable, data = input_l, value.var = c("value"))
    if (save_file) {
        output_file = gsub(".txt", "-wide.txt", input_file)
        write.table(input_w, file = output_file, eol = "\n", sep = "\t", quote = F, row.names = F, col.names = T)
    }
    if (return_results) {
        return(input_w)
    }
}
