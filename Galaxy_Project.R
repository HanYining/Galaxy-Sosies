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

#PCA in flux data
fluxData <- data3[,53:63]
fluxPCA <- prcomp(na.omit(fluxData), center = TRUE, scale = TRUE)
print(fluxPCA)
plot(fluxPCA, type = "l")

#First 4 PCA explains 54% of the variance
summary(fluxPCA)



#PCA post normalization (finding the highest peak fluxes)





