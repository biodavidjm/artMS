# Create data files

# Generate Random File
randomDF <- data.frame(replicate(10,sample(0:1,100,rep=TRUE)))

save(randomDF, file = '~/github/biodavidjm/artMS/data/randomDF.Rdata')

