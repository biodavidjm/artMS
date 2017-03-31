#' @import data.table
#' @title Check if an input is a file or a data object
#' @description This function is used in order to make it so a user can submit either a path to a data file or a data object in data.frame or data.table form.
#' @param input_file The filepath/object to be checked.
#' @param is.evidence Whether or not the file to be read in is an evidence file. This will assign proper classes to the evidence file when being read in.
#' @keywords file, evidence, input
#' checkIfFile()
#' @export
checkIfFile <- function(input_file, is.evidence=FALSE){
    # check if already a data.frame or data.table
  if(is.data.table(input_file) | is.data.frame(input_file)){
    x <-  data.table(input_file)
  }else if(tryCatch(file.exists(input_file), error=function(e) return(FALSE))){   # check if file path is legit
    if(is.evidence){
      x <- read_evidence_file(input_file)
    }else{
      x <- fread(input_file, integer64 = 'double')
    }
  }else{
    stop("There's something wrong with the file/object you submitted:\n\t", input_file,"\nPlease check that the file directory is correct, or that the data is in `data.frame` or `data.table` format.")
  }
  return(x)
}
