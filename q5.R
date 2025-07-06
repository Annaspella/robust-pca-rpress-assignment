library(robustbase)
library(rrcov)
library(vegan)
library(R.matlab)

# Analysis by reading the csv data set: table.csv
Data <- read.csv(file.choose(), header=TRUE)
head(Data)
dim(Data)
Data.new = Data[seq(1, nrow(Data), 3),]
head(Data.new)
dim(Data.new)


par(mfrow=c(1,2))

Data.std <- scale(Data.new, center=TRUE, scale=TRUE)
rob1 <- PcaHubert(Data.std, k=6)
#screeplot(rob, type="l", xlim=c(1,8))
#L <- loadings(rob)
#tolEllipsePlot(L[,2:3], classic=TRUE)
pca.scoreplot(rob1, i=1, j=2, id.n=TRUE)

clas1 <- PcaClassic(Data.std, k=2)
#screeplot(clas, type="l", xlim=c(1,8))
#LL <- clas$rotation
#tolEllipsePlot(LL[,2:3], classic=TRUE)
pca.scoreplot(clas1, i=1, j=2, id.n=TRUE)


plot(rob1)
plot(clas1)

# Analysis by reading the .mat file: GaitData.mat
data <- readMat(file.choose())
df <- data.frame(matrix(unlist(data), nrow=48, byrow=FALSE),stringsAsFactors=FALSE)
df.new <- df[seq(1,nrow(df),3),]

par(mfrow=c(1,1))

df.std <- scale(df.new, center=TRUE, scale=TRUE)
rob2 <- PcaHubert(df.std, k=6)
#screeplot(rob, type="l", xlim=c(1,8))
#L <- loadings(rob)
#tolEllipsePlot(L[,2:3], classic=TRUE)
pca.scoreplot(rob2, i=2, j=3, id.n=TRUE)

clas2 <- PcaClassic(df.std, k=4)
#screeplot(clas, type="l", xlim=c(1,8))
#LL <- clas$rotation
#tolEllipsePlot(LL[,2:3], classic=TRUE)
pca.scoreplot(clas2, i=2, j=3, id.n=TRUE)

plot(rob2, pch=16)
plot(clas2, pch=16)


