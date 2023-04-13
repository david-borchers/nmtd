library(rdryad)
library(readxl)

filepath <- dryad_download(dois = "10.5061/dryad.msbcc2fv7")

Karoodata <- read_excel(filepath[[1]], sheet = "TTD karoo bird data")

head(Karoodata)
save(Karoodata,file="./data/Karoodata.rda")
