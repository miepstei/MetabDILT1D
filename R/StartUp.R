a <- read.csv("/Users/me394/Dropbox/Academic/EBPOD/Metabolomics/DILT1DFileLocations.csv", header = TRUE, stringsAsFactors = FALSE)
fileList <- as.list(a$FILE_LOCATION)
names(fileList) <- a$FILE_VARIABLE_NAME
unlink(a)

