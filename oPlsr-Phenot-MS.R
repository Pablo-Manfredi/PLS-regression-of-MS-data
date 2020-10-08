#!/usr/bin/env Rscript
#If these packages are not installed yet
#Run the following hashed commands
##########################   START PACKAGE LOADING : UNHASH FOLLOWING COMMANDS THE FIRST TIME
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("pls")
library(pls)
#install.packages("data.table")
library(data.table)
##########################   END PACKAGE LOADING


##########################   START Load VIP functions from Bjørn-Helge Mevik at https://mevik.net/work/software/VIP.R
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))}

## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))}
##########################   END Load VIP functions from Bjørn-Helge Mevik at https://mevik.net/work/software/VIP.R


##############  Get the whole data from file & format as dataframe
phenotype = 'Persist' 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
`FileData` <- fread(file="Table4R.tsv", sep="\t", stringsAsFactors = TRUE, header = TRUE)
Mydatframe <- data.frame (phenotype=(I(scale(`FileData`[,2], center = FALSE, scale = FALSE))))
LengthOfData <- ncol(`FileData`)
MSdata<- data.frame (`FileData`[,3:LengthOfData])
# Standardizing the MS data proven far better predictive power here
Mydatframe$MS <- as.matrix ((I(scale(MSdata, center = TRUE, scale = TRUE))))
# No Standardisation of phenotypes for plots readability. Might be required for some datasets.
Mydatframe$Persist <- as.vector ((I(scale(`FileData`[,2], center = FALSE, scale = FALSE))))
rownames(Mydatframe)<- as.matrix(`FileData`[,1])


##############  The first pls regression with the whole dataset
ncompUsedHere = 2
Proteomic.plsR <- plsr(Persist ~ MS, ncomp = ncompUsedHere+2, data = Mydatframe, validation = "LOO", method = "oscorespls")
summary(Proteomic.plsR)
par(family = 'mono') #just a font adjustment
plot(RMSEP(Proteomic.plsR), legendpos = "topright", main=paste("Sum of Errors for", " Whole Data Set Models",sep=''))
plot(Proteomic.plsR, ncomp = ncompUsedHere, asp = 1, line = TRUE, main=paste("Leave One Out Cross Validation", " Cycle 0",sep=':'))#, labels="names")
plot(Proteomic.plsR, ncomp = ncompUsedHere, asp = 1, line = TRUE, labels="names", main=paste("Leave one out cross Validation", " cycle 0",sep=''))


##############  start of the FEATURE REDUCTION
VIPcutoff = 1
RowsData <- nrow(Mydatframe)
ColumnsData <- ncol(Mydatframe)
VIPvals <- as.data.frame(VIP (Proteomic.plsR))
# Print a file with VIPs from the first oplsr
write.table(t(VIPvals), file='Persist.VIP', sep="\t")
# Creat the list of variables for the next oplsr cycle
VIPvarAbove1 <- subset(VIPvals[ncompUsedHere,], select=(VIPvals[ncompUsedHere,]>VIPcutoff ))
Matching <-  ((colnames(Mydatframe$MS)) %in% (names(VIPvarAbove1)))
#Creat A new data frame with the desired subset of MS data
MyLOOPdatframe <- Mydatframe[]
MyLOOPdatframe$loopMS <- as.matrix(Mydatframe$MS[,select=Matching]) 


##############  Iterating cycles of FEATURE REDUCTION
x <- 0;
StatsCollector <- matrix(ncol = 3)
StatsCollector <- rbind(StatsCollector, c((x), Proteomic.plsR$validation$adj[2], dim(Proteomic.plsR$projection)[1]))
while(x < 15)# can be increased, usually before 10 cycles the data is exhausted 
{
  LOOP.Proteomic.plsR <- plsr(Persist ~ loopMS, ncomp = ncompUsedHere+2, data = MyLOOPdatframe, validation = "LOO", method = "oscorespls")
  print(paste(phenotype, " feature reduction cycle ", x+1, sep=''))
  summary(LOOP.Proteomic.plsR)
  StatsCollector <- rbind(StatsCollector, c((x+1), LOOP.Proteomic.plsR$validation$adj[2], dim(LOOP.Proteomic.plsR$projection)[1]))
  plot(RMSEP(LOOP.Proteomic.plsR), legendpos = "topright", main=paste("Sum of Errors for", " reduction cycle ", x+1, sep=''))
  plot(LOOP.Proteomic.plsR, ncomp = 2, asp = 1, line = TRUE, main=paste(phenotype, " feature reduction cycle ", x+1, sep=''))#, labels="names")
  plot(LOOP.Proteomic.plsR, ncomp = 2, asp = 1, line = TRUE, labels="names", main=paste(phenotype, " feature reduction cycle ", x+1, sep=''))
  TABFileOUT <- 'VIP_run_'
  xfak <- x+1;
  TABFileOUT <-paste(TABFileOUT, xfak, sep = "")
  TABFileOUT <-paste(TABFileOUT, VIPcutoff, sep = ".abo")
  TABFileOUT <-paste(TABFileOUT, ncompUsedHere, sep = ".comp")
  extension <- '.tab';
  TABFileOUT <-paste(TABFileOUT, extension, sep = "")
  TABFileOUT <-paste(TABFileOUT, "Values", sep = ".")
  write.table(MyLOOPdatframe$loopMS, file=TABFileOUT, sep="\t", row.names = row.names(MyLOOPdatframe),col.names=NA)
  VIPvalsLOOP <- as.data.frame(VIP (LOOP.Proteomic.plsR))
  VIPvarAbove1 <- subset(VIPvalsLOOP[ncompUsedHere,], select=(VIPvalsLOOP[ncompUsedHere,]>VIPcutoff ))
  Matching <-  ((colnames(MyLOOPdatframe$loopMS)) %in% (names(VIPvarAbove1)))
  #Creat A new data frame with the desired subset of MS data
  MyLOOPdatframe$loopMS <- data.matrix(MyLOOPdatframe$loopMS[,select=Matching]) 
  
  # Print The VIP file
  TABFileOUT <- 'VIP_run_'
  x <- x+1;
  TABFileOUT <-paste(TABFileOUT, x, sep = "")
  TABFileOUT <-paste(TABFileOUT, VIPcutoff, sep = ".abo")
  TABFileOUT <-paste(TABFileOUT, ncompUsedHere, sep = ".comp")
  extension <- '.tab';
  TABFileOUT <-paste(TABFileOUT, extension, sep = "")
  write.table(t(VIPvalsLOOP), file=TABFileOUT, sep="\t",col.names=NA)
}

##############  Plot the statistics of the Feature reduction (adapted from Robert W. Baer)
## aesthetics
par(mar=c(5, 4, 4, 6) + 0.1)
## Plot first set of data and draw its axis
plot(StatsCollector[,1], StatsCollector[,2], pch=16, axes=FALSE, type="b",col="black", main="Feature Reduction", xlab="", ylab="")
axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
mtext("CV Error score",side=2,line=3.5)
box()
## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(StatsCollector[,1], StatsCollector[,3], pch=15,  xlab="", ylab="", axes=FALSE, type="b", col="red")
## a little farther out (line=4) to make room for labels
mtext("Number of variables",side=4,col="red",line=4) 
axis(4, col="red",col.axis="red",las=1)

## Draw the time axis
axis(1,pretty(StatsCollector[,1]))
mtext("Cycles",side=1,col="black",line=2.5)  

## Add Legend
legend("topright",legend=c("Error","Variables"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

