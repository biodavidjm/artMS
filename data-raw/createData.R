# CREATE SAMPLE DATA FILES FOR THE PACKAGE

setwd("~/github/biodavidjm/artMS/")

# GENERATE RANDOM FILE
randomDF <- data.frame(replicate(10,sample(0:1,100,rep=TRUE)))
save(randomDF, file = 'data/randomDF.rdata')

# PH FILES
# MaxQuant Evidence file
ph_evidence <- read.delim('data-raw/ph/evidence.txt', stringsAsFactors = F)
save(ph_keys, file='data/ph_keys.rdata')

# Keys file (experimental design)
ph_keys <- read.delim("data-raw/ph/keys.txt", stringsAsFactors = F)
save(ph_evidence, file='data/ph_evidence.rdata')
