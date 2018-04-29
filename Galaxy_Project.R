#Code so far 

#Reading in the Data in R
system.time(newData <- fread("galaxie_data.csv", header = TRUE))
newData <- replace(newData, newData == "-9999", NA)
newData <- data.frame(newData)

#Removes column if more than 25% of 0's
deleteCol <- c()
for(i in 1:ncol(newData)) {
  if(length(which(newData[,i] == 0)) > 0.25*nrow(newData)) {
    deleteCol <- append(deleteCol,i) 
  }
}
data2 <- newData[,-deleteCol]

#Quick summary of data2
head(data2)
str(data2)

#Binning redshifts in 10 bins
data2Bins <- .bincode(data2$specz, breaks = quantile(data2$specz, seq(0,1, by = 0.1)), 
                      include.lowest = TRUE)
data3 <- cbind(data2, data2Bins)
head(data3)

#PCA in flux data (not normalized)
fluxData <- data3[,53:63]
fluxPCA <- prcomp(na.omit(fluxData), center = TRUE, scale = TRUE)
print(fluxPCA)
plot(fluxPCA, type = "l")
summary(fluxPCA)

#PCA on photometric data
library(dplyr)
library(ggfortify)
photoData <- data3[,c(63,8:52)]

#PCA on bin 1
bin1 <- filter(photoData, photoData$data2Bins == 1)
bin1PCA <- prcomp(na.omit(bin1[,2:45]), center = TRUE, scale = TRUE)
print(bin1PCA)
#plot(bin1PCA, type = "l")
summary(bin1PCA)
autoplot(bin1PCA, loadings = TRUE)

#PCA on bin 2
bin2 <- filter(photoData, photoData$data2Bins == 2)
bin2PCA <- prcomp(na.omit(bin2[,2:45]), center = TRUE, scale = TRUE)
print(bin2PCA)
#plot(bin1PCA, type = "l")
summary(bin2PCA)
autoplot(bin2PCA, loadings = TRUE)

#PCA on bin 3
bin3 <- filter(photoData, photoData$data2Bins == 3)
bin3PCA <- prcomp(na.omit(bin3[,2:45]), center = TRUE, scale = TRUE)
print(bin3PCA)
#plot(bin1PCA, type = "l")
summary(bin3PCA)
autoplot(bin3PCA, loadings = TRUE)

#PCA on bin 4
bin4 <- filter(photoData, photoData$data2Bins == 4)
bin4PCA <- prcomp(na.omit(bin4[,2:45]), center = TRUE, scale = TRUE)
print(bin4PCA)
#plot(bin1PCA, type = "l")
summary(bin4PCA)
autoplot(bin4PCA, loadings = TRUE)


#PCA on bin 5
bin5 <- filter(photoData, photoData$data2Bins == 5)
bin5PCA <- prcomp(na.omit(bin5[,2:45]), center = TRUE, scale = TRUE)
print(bin5PCA)
#plot(bin1PCA, type = "l")
summary(bin5PCA)
autoplot(bin5PCA, loadings = TRUE)



