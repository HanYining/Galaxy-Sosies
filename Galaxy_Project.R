#Reading in the Data in R
library(data.table)
system.time(newData <- fread("galaxie_data.csv", header = TRUE))
newData <- replace(newData, newData == "-9999", NA)
newData <- data.frame(newData)


data2Bins <- .bincode(newData$specz, breaks = quantile(newData$specz, seq(0,1, by = 0.1)), 
                      include.lowest = TRUE)
newData <- cbind(newData, data2Bins)
head(newData)
str(newData)

#Removes column if more than 25% of 0's from each unique bin

deleteCol <- c()
for(i in 1:unique(newData$data2Bins)) {
  for(j in 1:ncol(newData)) {
    if(length(which(newData[,j] == 0)) > 0.25*nrow(newData)) {
      deleteCol <- append(deleteCol,j)
    }
  }
}
data2 <- newData[,-deleteCol]

#Separating Flux Data and Photometric Data
fluxData <- data2[,c(2,4,53:63)]
photoData <- data2[,c(2,63, 8:52)]

#Quick summary of data2
head(data2)
str(data2)

#PCA on Photometric Data by bin that spits out 

lapply(photoData[,c(3:47)], FUN = class)
par(mfrow = c(2,2), mar=c(4,4,2,1))
by(photoData, INDICES = photoData$data2Bins, FUN = function(z) {
  pca <- prcomp(z[,unlist(lapply(z, FUN = class)) == "numeric"], scale = TRUE, center = TRUE)
  #plot(pca$x[,1:2], pch=16, main=z[1,"Photometric"])
  summary(pca)
  z
})



