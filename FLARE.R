# =============================================================================

###### !!!!!!!!!!!!! run DATASETs.R first !!!!!!!!!!!!! #####

# =============================================================================
# Health Prognostics in Multi-sensor Systems (FLARE Methodology)
#
# This script implements the entire analysis pipeline presented in the paper:
# "Health Prognostics in Multi-sensor Systems Based on Multivariate Functional
# Data Analysis"
#
# The script will:
# 1.  INSTALL PACKAGES
# 2.  LOAD LIBRARIES
# 3.  SMOOTHING
# 4.  GENERATE MULTIVARIATE FUNCTIONAL DATA
# 5.  PERFORM MULTIVARIATE FUCNTIONAL PRINCIPAL COMPONENT ANALYSIS (MFPCA)
# 6.  CLASSIFY TRAIN DATA GROUPS
# 7.  CLASSIFY TEST DATA GROUPS (YOUDEN INDEX)
# 8.  CALCULATE DISTANCES FOR DIFFERENT GROUPS
# 9.  CALCULATE MULTIVARIATE EUCLIDIAN DISTANCE FOR LOW AND HIGH SCORE GROUPS
# 10. CALCULATE RMSEs
# =============================================================================

# -----------------------------------------------------------------------------
# 1. INSTALL PACKAGES
# -----------------------------------------------------------------------------

# Load required packages
install.packages("funData")
install.packages("ggplot2")
install.packages("fda")
install.packages("devtools")
install.packages("philentropy")
install.packages("multimode")
install.packages("PredictionR")
install.packages("LaplacesDemon")
install.packages("diptest")
install.packages("mousetrap")
install.packages("cutoff")
install.packages("bbmle")

# -----------------------------------------------------------------------------
# 2. LOAD LIBRARIES
# -----------------------------------------------------------------------------

# Load required libraries
library(funData)
library(ggplot2)
library(fda)
library(MFPCA)
library(devtools)
library(philentropy)
library(multimode)
library(PredictionR)
library(LaplacesDemon)
library(diptest)
library(mousetrap)
library(cutoff)
library(bbmle)

# -----------------------------------------------------------------------------
# 3. SMOOTHING
# -----------------------------------------------------------------------------

tFD001_T24_train <-  t(FD001_T24_train)
tFD001_T30_train <-  t(FD001_T30_train)
tFD001_T50_train <-  t(FD001_T50_train)
tFD001_T24_train <-  t(FD001_T24_train)
tFD001_P30_train <-  t(FD001_P30_train)
tFD001_ps30_train <- t(FD001_ps30_train)
tFD001_phi_train <-  t(FD001_phi_train)
tFD001_BPR_train <-  t(FD001_BPR_train)
tFD001_W31_train <-  t(FD001_W31_train)
tFD001_W32_train <-  t(FD001_W32_train)

# Check argvals for registration
FD001argvalssmooth

# Visualize discrete data
plot(FD001argvalssmooth,tFD001_T24_train, 
     ylab="T30 Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001argvalssmooth,tFD001_T30_train, 
     ylab="T30 Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001argvalssmooth,tFD001_T50_train,
     ylab="T50 Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001argvalssmooth,tFD001_P30_train,
     ylab="P30 Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001argvalssmooth,tFD001_ps30_train,
     ylab="ps30 Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001argvalssmooth,tFD001_phi_train,
     ylab="phi Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001argvalssmooth,tFD001_BPR_train,
     ylab="BPR Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001argvalssmooth,tFD001_W31_train,
     ylab="W31 Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001argvalssmooth,tFD001_W32_train,
     ylab="W32 Sensor Value", 
     xlab="All 100 Engines observations registered between [0-1] interval")
#####################################

dim(FD001argvalssmooth)
is.numeric(FD001argvalssmooth)

#Generate Bspline basis
bsplinebasis2<- create.bspline.basis(c(0,1), 8)
bsplinebasis2
plot(bsplinebasis2)

# Generate MEANS
FD001meanT24<-mean(na.omit(as.vector(tFD001_T24_train)))
FD001meanT30<-mean(na.omit(as.vector(tFD001_T30_train)))
FD001meanT50<-mean(na.omit(as.vector(tFD001_T50_train)))
FD001meanP30<-mean(na.omit(as.vector(tFD001_P30_train)))
FD001meanps30<-mean(na.omit(as.vector(tFD001_ps30_train)))
FD001meanphi<-mean(na.omit(as.vector(tFD001_phi_train)))
FD001meanBPR<-mean(na.omit(as.vector(tFD001_BPR_train)))
FD001meanW31<-mean(na.omit(as.vector(tFD001_W31_train)))
FD001meanW32<-mean(na.omit(as.vector(tFD001_W32_train)))

# Smooth data
FD001smoothallT24 <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_T24_train[,i])))
  FD001smoothallT24[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_T24_train[,i]))/FD001meanT24,bsplinebasis2)[["fd"]][["coefs"]]
}
FD001smoothallT30 <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_T30_train[,i])))
  FD001smoothallT30[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_T30_train[,i]))/FD001meanT30,bsplinebasis2)[["fd"]][["coefs"]]
}
FD001smoothallT50 <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_T50_train[,i])))
  FD001smoothallT50[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_T50_train[,i]))/FD001meanT50,bsplinebasis2)[["fd"]][["coefs"]]
}
FD001smoothallP30 <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_P30_train[,i])))
  FD001smoothallP30[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_P30_train[,i]))/FD001meanP30,bsplinebasis2)[["fd"]][["coefs"]]
}
FD001smoothallps30 <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_ps30_train[,i])))
  FD001smoothallps30[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_ps30_train[,i]))/FD001meanps30,bsplinebasis2)[["fd"]][["coefs"]]
}
FD001smoothallphi <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_phi_train[,i])))
  FD001smoothallphi[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_phi_train[,i]))/FD001meanphi,bsplinebasis2)[["fd"]][["coefs"]]
}
FD001smoothallBPR <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_BPR_train[,i])))
  FD001smoothallBPR[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_BPR_train[,i]))/FD001meanBPR,bsplinebasis2)[["fd"]][["coefs"]]
}
FD001smoothallW31 <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_W31_train[,i])))
  FD001smoothallW31[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_W31_train[,i]))/FD001meanW31,bsplinebasis2)[["fd"]][["coefs"]]
}
FD001smoothallW32 <- matrix(data = NA, nrow=8, ncol=length(FD001_train))
for (i in 1:length(FD001_train)) {
  ENGINEallArgvals <-seq(0,1, length.out= length(na.omit(tFD001_W32_train[,i])))
  FD001smoothallW32[,i]=smooth.basis(ENGINEallArgvals,as.vector(na.omit(tFD001_W32_train[,i]))/FD001meanW32,bsplinebasis2)[["fd"]][["coefs"]]
}

####################################
newargvals <- seq(0,1,length=8)
#################################

FD001fdsmoothallT24<-Data2fd(argvals = newargvals, y=FD001smoothallT24)
FD001fdsmoothallT24$coefs<-FD001smoothallT24
FD001fdsmoothallT24$basis$nbasis<-8
FD001fdsmoothallT24$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallT24$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
FD001fdsmoothallT24

plot(FD001fdsmoothallT24)

FD001fdsmoothallT30<-Data2fd(argvals = newargvals, y=FD001smoothallT30)
FD001fdsmoothallT30$coefs<-FD001smoothallT30
FD001fdsmoothallT30$basis$nbasis<-8
FD001fdsmoothallT30$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallT30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD001fdsmoothallT50<-Data2fd(argvals = newargvals, y=FD001smoothallT50)
FD001fdsmoothallT50$coefs<-FD001smoothallT50
FD001fdsmoothallT50$basis$nbasis<-8
FD001fdsmoothallT50$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallT50$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD001fdsmoothallP30<-Data2fd(argvals = newargvals, y=FD001smoothallP30)
FD001fdsmoothallP30$coefs<-FD001smoothallP30
FD001fdsmoothallP30$basis$nbasis<-8
FD001fdsmoothallP30$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallP30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD001fdsmoothallps30<-Data2fd(argvals = newargvals, y=FD001smoothallps30)
FD001fdsmoothallps30$coefs<-FD001smoothallps30
FD001fdsmoothallps30$basis$nbasis<-8
FD001fdsmoothallps30$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallps30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD001fdsmoothallphi<-Data2fd(argvals = newargvals, y=FD001smoothallphi)
FD001fdsmoothallphi$coefs<-FD001smoothallphi
FD001fdsmoothallphi$basis$nbasis<-8
FD001fdsmoothallphi$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallphi$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD001fdsmoothallBPR<-Data2fd(argvals = newargvals, y=FD001smoothallBPR)
FD001fdsmoothallBPR$coefs<-FD001smoothallBPR
FD001fdsmoothallBPR$basis$nbasis<-8
FD001fdsmoothallBPR$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallBPR$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD001fdsmoothallW31<-Data2fd(argvals = newargvals, y=FD001smoothallW31)
FD001fdsmoothallW31$coefs<-FD001smoothallW31
FD001fdsmoothallW31$basis$nbasis<-8
FD001fdsmoothallW31$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallW31$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD001fdsmoothallW32<-Data2fd(argvals = newargvals, y=FD001smoothallW32)
FD001fdsmoothallW32$coefs<-FD001smoothallW32
FD001fdsmoothallW32$basis$nbasis<-8
FD001fdsmoothallW32$basis$params<- c(0.2,0.4,0.6,0.8)
FD001fdsmoothallW32$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

#############################

FD002fdsmoothallT24<-Data2fd(argvals = newargvals2, y=FD002smoothallT24)
FD002fdsmoothallT24$coefs<-FD002smoothallT24
FD002fdsmoothallT24$basis$nbasis<-30
FD002fdsmoothallT24$basis$params<- c(0.03703704,0.07407407,0.11111111,0.14814815,
                                     0.18518519,0.22222222,0.25925926,0.29629630,
                                     0.33333333,0.37037037,0.40740741,0.4444444,
                                     0.48148148,0.51851852,0.55555556,0.59259259,
                                     0.62962963,0.66666667,0.70370370,0.74074074,
                                     0.77777778,0.81481481,0.85185185,0.88888889,
                                     0.92592593,0.96296296)
FD002fdsmoothallT24$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8","bspl4.9","bspl4.10",
                                   "bspl4.11","bspl4.12","bspl4.13","bspl4.14","bspl4.15","bspl4.16","bspl4.17","bspl4.18","bspl4.19","bspl4.20",
                                   "bspl4.21","bspl4.22","bspl4.23","bspl4.24","bspl4.25","bspl4.26","bspl4.27","bspl4.28","bspl4.29","bspl4.30")
plot(FD002fdsmoothallT24)

FD002fdsmoothallT24<-Data2fd(argvals = newargvals, y=FD002smoothallT24)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallT24$coefs<-FD002smoothallT24
FD002fdsmoothallT24$basis$nbasis<-8
FD002fdsmoothallT24$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallT24$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
plot(FD002fdsmoothallT24)

FD002fdsmoothallT30<-Data2fd(argvals = newargvals, y=FD002smoothallT30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallT30$coefs<-FD002smoothallT30
FD002fdsmoothallT30$basis$nbasis<-8
FD002fdsmoothallT30$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallT30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD002fdsmoothallT50<-Data2fd(argvals = newargvals, y=FD002smoothallT50)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallT50$coefs<-FD002smoothallT50
FD002fdsmoothallT50$basis$nbasis<-8
FD002fdsmoothallT50$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallT50$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD002fdsmoothallP30<-Data2fd(argvals = newargvals, y=FD002smoothallP30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallP30$coefs<-FD002smoothallP30
FD002fdsmoothallP30$basis$nbasis<-8
FD002fdsmoothallP30$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallP30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD002fdsmoothallps30<-Data2fd(argvals = newargvals, y=FD002smoothallps30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallps30$coefs<-FD002smoothallps30
FD002fdsmoothallps30$basis$nbasis<-8
FD002fdsmoothallps30$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallps30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD002fdsmoothallphi<-Data2fd(argvals = newargvals, y=FD002smoothallphi)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallphi$coefs<-FD002smoothallphi
FD002fdsmoothallphi$basis$nbasis<-8
FD002fdsmoothallphi$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallphi$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD002fdsmoothallBPR<-Data2fd(argvals = newargvals, y=FD002smoothallBPR)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallBPR$coefs<-FD002smoothallBPR
FD002fdsmoothallBPR$basis$nbasis<-8
FD002fdsmoothallBPR$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallBPR$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD002fdsmoothallW31<-Data2fd(argvals = newargvals, y=FD002smoothallW31)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallW31$coefs<-FD002smoothallW31
FD002fdsmoothallW31$basis$nbasis<-8
FD002fdsmoothallW31$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallW31$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD002fdsmoothallW32<-Data2fd(argvals = newargvals, y=FD002smoothallW32)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD002fdsmoothallW32$coefs<-FD002smoothallW32
FD002fdsmoothallW32$basis$nbasis<-8
FD002fdsmoothallW32$basis$params<- c(0.2,0.4,0.6,0.8)
FD002fdsmoothallW32$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
#############################

FD003fdsmoothallT24<-Data2fd(argvals = newargvals, y=FD003smoothallT24)
FD003fdsmoothallT24$coefs<-FD003smoothallT24
FD003fdsmoothallT24$basis$nbasis<-8
FD003fdsmoothallT24$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallT24$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
FD003fdsmoothallT24

plot(FD003fdsmoothallT24, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="T24",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)

FD003fdsmoothallT30<-Data2fd(argvals = newargvals, y=FD003smoothallT30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD003fdsmoothallT30$coefs<-FD003smoothallT30
FD003fdsmoothallT30$basis$nbasis<-8
FD003fdsmoothallT30$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallT30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

plot(FD003fdsmoothallT30, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="T30",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)

FD003fdsmoothallT50<-Data2fd(argvals = newargvals, y=FD003smoothallT50)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD003fdsmoothallT50$coefs<-FD003smoothallT50
FD003fdsmoothallT50$basis$nbasis<-8
FD003fdsmoothallT50$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallT50$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
plot(FD003fdsmoothallT50, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="T50",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)

FD003fdsmoothallP30<-Data2fd(argvals = newargvals, y=FD003smoothallP30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD003fdsmoothallP30$coefs<-FD003smoothallP30
FD003fdsmoothallP30$basis$nbasis<-8
FD003fdsmoothallP30$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallP30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
plot(FD003fdsmoothallP30, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="P30",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)


FD003fdsmoothallps30<-Data2fd(argvals = newargvals, y=FD003smoothallps30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD003fdsmoothallps30$coefs<-FD003smoothallps30
FD003fdsmoothallps30$basis$nbasis<-8
FD003fdsmoothallps30$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallps30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
plot(FD003fdsmoothallps30, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="ps30",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)

FD003fdsmoothallphi<-Data2fd(argvals = newargvals, y=FD003smoothallphi)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD003fdsmoothallphi$coefs<-FD003smoothallphi
FD003fdsmoothallphi$basis$nbasis<-8
FD003fdsmoothallphi$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallphi$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
plot(FD003fdsmoothallphi, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="phi",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)

FD003fdsmoothallBPR<-Data2fd(argvals = newargvals, y=FD003smoothallBPR)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD003fdsmoothallBPR$coefs<-FD003smoothallBPR
FD003fdsmoothallBPR$basis$nbasis<-8
FD003fdsmoothallBPR$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallBPR$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
plot(FD003fdsmoothallBPR, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="BPR",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)

FD003fdsmoothallW31<-Data2fd(argvals = newargvals, y=FD003smoothallW31)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD003fdsmoothallW31$coefs<-FD003smoothallW31
FD003fdsmoothallW31$basis$nbasis<-8
FD003fdsmoothallW31$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallW31$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
plot(FD003fdsmoothallW31, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="W31",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)

FD003fdsmoothallW32<-Data2fd(argvals = newargvals, y=FD003smoothallW32)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD003fdsmoothallW32$coefs<-FD003smoothallW32
FD003fdsmoothallW32$basis$nbasis<-8
FD003fdsmoothallW32$basis$params<- c(0.2,0.4,0.6,0.8)
FD003fdsmoothallW32$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
plot(FD003fdsmoothallW32, xlab="Cycle Time, registered interval", ylab="Sensor Value", lty=1,
     main="W32",
     cex.lab=2.4, cex.axis=2.2,cex.main=3, cex.sub=2)

#############################

FD004fdsmoothallT24<-Data2fd(argvals = newargvals, y=FD004smoothallT24)
FD004fdsmoothallT24$coefs<-FD004smoothallT24
FD004fdsmoothallT24$basis$nbasis<-8
FD004fdsmoothallT24$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallT24$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
FD004fdsmoothallT24

plot(FD004fdsmoothallT24)

FD004fdsmoothallT30<-Data2fd(argvals = newargvals, y=FD004smoothallT30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD004fdsmoothallT30$coefs<-FD004smoothallT30
FD004fdsmoothallT30$basis$nbasis<-8
FD004fdsmoothallT30$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallT30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD004fdsmoothallT50<-Data2fd(argvals = newargvals, y=FD004smoothallT50)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD004fdsmoothallT50$coefs<-FD004smoothallT50
FD004fdsmoothallT50$basis$nbasis<-8
FD004fdsmoothallT50$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallT50$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD004fdsmoothallP30<-Data2fd(argvals = newargvals, y=FD004smoothallP30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD004fdsmoothallP30$coefs<-FD004smoothallP30
FD004fdsmoothallP30$basis$nbasis<-8
FD004fdsmoothallP30$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallP30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD004fdsmoothallps30<-Data2fd(argvals = newargvals, y=FD004smoothallps30)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD004fdsmoothallps30$coefs<-FD004smoothallps30
FD004fdsmoothallps30$basis$nbasis<-8
FD004fdsmoothallps30$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallps30$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD004fdsmoothallphi<-Data2fd(argvals = newargvals, y=FD004smoothallphi)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD004fdsmoothallphi$coefs<-FD004smoothallphi
FD004fdsmoothallphi$basis$nbasis<-8
FD004fdsmoothallphi$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallphi$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD004fdsmoothallBPR<-Data2fd(argvals = newargvals, y=FD004smoothallBPR)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD004fdsmoothallBPR$coefs<-FD004smoothallBPR
FD004fdsmoothallBPR$basis$nbasis<-8
FD004fdsmoothallBPR$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallBPR$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD004fdsmoothallW31<-Data2fd(argvals = newargvals, y=FD004smoothallW31)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD004fdsmoothallW31$coefs<-FD004smoothallW31
FD004fdsmoothallW31$basis$nbasis<-8
FD004fdsmoothallW31$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallW31$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")

FD004fdsmoothallW32<-Data2fd(argvals = newargvals, y=FD004smoothallW32)
###çok önemli smoothingi fd ye çevirdikten sonra düzeltme gerekiyor!!!!
FD004fdsmoothallW32$coefs<-FD004smoothallW32
FD004fdsmoothallW32$basis$nbasis<-8
FD004fdsmoothallW32$basis$params<- c(0.2,0.4,0.6,0.8)
FD004fdsmoothallW32$basis$names<-c("bspl4.1","bspl4.2","bspl4.3","bspl4.4","bspl4.5","bspl4.6","bspl4.7","bspl4.8")
#############################

#Generate Smooth Plots (fd objects)
plot(FD001fdsmoothallT24,
     ylab="T24 / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001fdsmoothallT30,
     ylab="T30 / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001fdsmoothallT50,
     ylab="T50 / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001fdsmoothallP30,
     ylab="P30 / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001fdsmoothallps30,
     ylab="ps30 / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001fdsmoothallphi,
     ylab="phi / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001fdsmoothallBPR,
     ylab="BPR / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001fdsmoothallW31,
     ylab="W31 / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
plot(FD001fdsmoothallW32,
     ylab="W32 / Smoothed Spline functions",
     xlab="All 100 Engines observations registered between [0-1] interval")
###############################

###############################
newargvals3 <- seq(0,1,length=200) ##tryig because  fd has 12 coefs
###############################

#Generate funData objects
FD001funDatasmoothallT24<-fd2funData(fdobj=FD001fdsmoothallT24,argvals =newargvals3)
FD001funDatasmoothallT30<-fd2funData(FD001fdsmoothallT30,argvals = newargvals3)
FD001funDatasmoothallT50<-fd2funData(FD001fdsmoothallT50,argvals = newargvals3)
FD001funDatasmoothallP30<-fd2funData(FD001fdsmoothallP30,argvals = newargvals3)
FD001funDatasmoothallps30<-fd2funData(FD001fdsmoothallps30,argvals = newargvals3)
FD001funDatasmoothallphi<-fd2funData(FD001fdsmoothallphi,argvals = newargvals3)
FD001funDatasmoothallBPR<-fd2funData(FD001fdsmoothallBPR,argvals = newargvals3)
FD001funDatasmoothallW31<-fd2funData(FD001fdsmoothallW31,argvals = newargvals3)
FD001funDatasmoothallW32<-fd2funData(FD001fdsmoothallW32,argvals = newargvals3)

autoplot(FD001funDatasmoothallT24)
class(FD001funDatasmoothallT24)

###copare is ok !!
plot(FD001fdsmoothallT24, main = "fd object")
plot(FD001funDatasmoothallT24, main = "funData object")

##check all
autoplot(FD001funDatasmoothallT30)
autoplot(FD001funDatasmoothallT50)
autoplot(FD001funDatasmoothallP30)
autoplot(FD001funDatasmoothallps30)
autoplot(FD001funDatasmoothallphi)
autoplot(FD001funDatasmoothallBPR)
autoplot(FD001funDatasmoothallW31)
autoplot(FD001funDatasmoothallW32)


# -----------------------------------------------------------------------------
# 4. GENERATE MULTIVARIATE FUNCTIONAL DATA
# -----------------------------------------------------------------------------

#Generate multiFunData objects for MFPCA
FD001multifunallsensors<-multiFunData(FD001funDatasmoothallT24,FD001funDatasmoothallT30,FD001funDatasmoothallT50,
                                      FD001funDatasmoothallP30,FD001funDatasmoothallps30,FD001funDatasmoothallphi,
                                      FD001funDatasmoothallBPR,FD001funDatasmoothallW31,FD001funDatasmoothallW32)
FD001multifunallsensors
dimSupp(FD001multifunallsensors)


# -----------------------------------------------------------------------------
# 5. PERFORM MULTIVARIATE FUCNTIONAL PRINCIPAL COMPONENT ANALYSIS (MFPCA)
# -----------------------------------------------------------------------------

#Perform MFPCA
FD001BsplineRegisteredsmoothedMFPCA <- MFPCA(FD001multifunallsensors, M = 3, 
                                        uniExpansions = BSplineExpansions, fit = TRUE)
FD001BsplineRegisteredsmoothedMFPCA

#Generate meanFunction plots
 univMeanFunctions <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$meanFunction)
 univMeanFunctions
 univMeanFunctions[[1]]
 univMeanFunctions[[2]]
 univMeanFunctions[[3]]
 univMeanFunctions[[4]]
 univMeanFunctions[[5]]
 univMeanFunctions[[6]]
 univMeanFunctions[[7]]
 univMeanFunctions[[8]]
 univMeanFunctions[[9]]

#Generate SCREEPLOTS
screeplot(FD001BsplineRegisteredsmoothedMFPCA, main = "Screeplot - lines")

screeplot(FD001BsplineRegisteredsmoothedMFPCA, type = "barplot", main = "Screeplot - barplot")

####################################################################

T24g <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[1]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
T24g

T30g <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[2]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
T30g

T50g <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[3]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
T50g

P30g <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[4]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
P30g

ps30g <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[5]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
ps30g

phig <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[6]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
phig

BPRg <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[7]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
BPRg

W31g <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[8]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
W31g

W32g <- autoplot(FD001BsplineRegisteredsmoothedMFPCA$functions[[9]]) + geom_hline(yintercept = 0, col = "grey", lwd = 1.5) + geom_line(aes(colour = obs), lwd = 1.25)+ labs(x = "Cycle Time") + scale_color_manual(values = c("#B80C0C", "#0C0CB8", "#129412"), name = "Princ. Comp.", labels= c("1st", "2nd", "3rd")) 
W32g


get_legend<-function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

gridExtra::grid.arrange(T24g + theme(legend.position = "none"), 
                        T30g + theme(legend.position = "none"), 
                        T50g + theme(legend.position = "none"), 
                        P30g + theme(legend.position = "none"), 
                        ps30g + theme(legend.position = "none"), 
                        phig + theme(legend.position = "none"),
                        BPRg + theme(legend.position = "none"), 
                        W31g + theme(legend.position = "none"), 
                        W32g + theme(legend.position = "none"),
                        nrow= 3, ncol = 3, widths = c(1, 1, 1))

#######################
#Generate dataframe
FD001df <- rbind(data.frame(PC = 1, val = FD001BsplineRegisteredsmoothedMFPCA$scores[, 1]), 
            data.frame(PC = 2, val = FD001BsplineRegisteredsmoothedMFPCA$scores[, 2]), 
            data.frame(PC = 3, val = FD001BsplineRegisteredsmoothedMFPCA$scores[, 3]))

#Check Scores
FD001BsplineRegisteredsmoothedMFPCA$scores

#Generate Score Plots
scoreplot(FD001BsplineRegisteredsmoothedMFPCA)

# -----------------------------------------------------------------------------
# 6. CLASSIFY TRAIN DATA GROUPS
# -----------------------------------------------------------------------------

#Check score densities
FD001dfpc1 <- FD001df[,2]
FD001dfpc1_100 <- FD001dfpc1[1:length(FD001_train)]
FD001dfpc1
FD001transformeddfpc1<- (1/(FD001dfpc1_100+1000))
FD001transformeddfpc1

hist(FD001dfpc1)
hist(FD001transformeddfpc1)

##### histogram of first PC score & density find using MLE ####
hist(FD001dfpc1, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "1/(FPC1+1000)",
     main = "")
lines(density(FD001dfpc1), # density plot
      lwd = 2, # thickness of line
      col = "red")

##### histogram of first PC score & density find using MLE ####
hist(FD001transformeddfpc1, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "1/(FPC1+1000)",
     main = "")
lines(density(FD001transformeddfpc1), # density plot
      lwd = 2, # thickness of line
      col = "red")

modetest(FD001transformeddfpc1)
####  In the above analysis, the p-value 
####  is virtually 0 — supporting the alternative 
####  hypothesis of more than one mode present in the 
####  data at the 5% level of significance.

FD001antimode<-locmodes(FD001transformeddfpc1,mod0=2,display=TRUE)
FD001antimode<-FD001antimode$locations[2]
FD001antimode

FD001lowscoredata<-FD001transformeddfpc1[which(FD001transformeddfpc1>FD001antimode )]
FD001bigscoredata<-FD001transformeddfpc1[which(FD001transformeddfpc1<FD001antimode )]
length(FD001bigscoredata)

hist(FD001transformeddfpc1, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "1st FPC",
     main = "",
)
lines(density(FD001bigscoredata), # density plot
      lwd = 2, # thickness of line
      col = "red")
lines(density(FD001lowscoredata), # density plot
      lwd = 2, # thickness of line
      col = "red")
	  
	  
# -----------------------------------------------------------------------------
# 7. CLASSIFY TEST DATA GROUPS (YOUDEN INDEX)
# -----------------------------------------------------------------------------



DATAT24<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_T24.csv", header = TRUE, row.names = 1)

x<-DATAT24$Initial.value.of.Sensor.1[DATAT24$Classification..MFPC.Scores.=="Low"]
y<-DATAT24$Initial.value.of.Sensor.1[DATAT24$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(G-F)
which((G-F)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((G-F)==J)
c=gridxy[k[1]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

T24seperation<-c

T24_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/1-T24train_test_all.csv", header = TRUE,row.names = 1)
T24_all_train_test
T24_all_train_test<-as.matrix(T24_all_train_test)
 
T24TESTinit<-matrix(NA, nrow = 100 , ncol = 2)
T24TESTinit[,2]<-as.vector(T24_all_train_test[101:200,2])
T24TESTinit[,1]<-seq(1:100)

T24TESTinitlow<-T24TESTinit[,1][which(T24TESTinit[,2]>T24seperation)]  
T24TESTinitbig<-T24TESTinit[,1][which(T24TESTinit[,2]<T24seperation)]  

#############
#############
#############
#############

DATAT30<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_T30.csv", header = TRUE, row.names = 1)

x<-DATAT30$Initial.value.of.Sensor.1[DATAT30$Classification..MFPC.Scores.=="Low"]
y<-DATAT30$Initial.value.of.Sensor.1[DATAT30$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(G-F)
which((G-F)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((G-F)==J)
c=gridxy[k[20]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

T30seperation<-c

T30_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/2-T30train_test_all.csv", header = TRUE,row.names = 1)
T30_all_train_test
T30_all_train_test<-as.matrix(T30_all_train_test)
 
T30TESTinit<-matrix(NA, nrow = 100 , ncol = 2)
T30TESTinit[,2]<-as.vector(T30_all_train_test[101:200,2])
T30TESTinit[,1]<-seq(1:100)

T30TESTinitlow<-T30TESTinit[,1][which(T30TESTinit[,2]>T30seperation)]  
T30TESTinitbig<-T30TESTinit[,1][which(T30TESTinit[,2]<T30seperation)]  

#############
#############
#############
#############

DATAT50<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_T50.csv", header = TRUE, row.names = 1)

x<-DATAT50$Initial.value.of.Sensor.1[DATAT50$Classification..MFPC.Scores.=="Low"]
y<-DATAT50$Initial.value.of.Sensor.1[DATAT50$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(G-F)
which((G-F)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((G-F)==J)
c=gridxy[k[10]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

T50seperation<-c

T50_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/3-T50train_test_all.csv", header = TRUE,row.names = 1)
T50_all_train_test
T50_all_train_test<-as.matrix(T50_all_train_test)
 
T50TESTinit<-matrix(NA, nrow = 100 , ncol = 2)
T50TESTinit[,2]<-as.vector(T50_all_train_test[101:200,2])
T50TESTinit[,1]<-seq(1:100)

T50TESTinitlow<-T50TESTinit[,1][which(T50TESTinit[,2]>T50seperation)]  
T50TESTinitbig<-T50TESTinit[,1][which(T50TESTinit[,2]<T50seperation)]  

#############
#############
#############
#############

DATAP30<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_P30.csv", header = TRUE, row.names = 1)

x<-DATAP30$Initial.value.of.Sensor.1[DATAP30$Classification..MFPC.Scores.=="Low"]
y<-DATAP30$Initial.value.of.Sensor.1[DATAP30$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(F-G)
which((F-G)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((F-G)==J)
c=gridxy[k[1]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

P30seperation<-c

P30_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/4-P30train_test_all.csv", header = TRUE,row.names = 1)
P30_all_train_test
P30_all_train_test<-as.matrix(P30_all_train_test)
 
P30TESTinit<-matrix(NA, nrow = 100 , ncol = 2)
P30TESTinit[,2]<-as.vector(P30_all_train_test[101:200,2])
P30TESTinit[,1]<-seq(1:100)

P30TESTinitlow<-P30TESTinit[,1][which(P30TESTinit[,2]<P30seperation)]  
P30TESTinitbig<-P30TESTinit[,1][which(P30TESTinit[,2]>P30seperation)]  

#############
#############
#############
#############

DATAps30<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_ps30.csv", header = TRUE, row.names = 1)

x<-DATAps30$Initial.value.of.Sensor.1[DATAps30$Classification..MFPC.Scores.=="Low"]
y<-DATAps30$Initial.value.of.Sensor.1[DATAps30$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(G-F)
which((G-F)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((G-F)==J)
c=gridxy[k[35]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

ps30seperation<-c

ps30_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/5-ps30train_test_all.csv", header = TRUE,row.names = 1)
ps30_all_train_test
ps30_all_train_test<-as.matrix(ps30_all_train_test)
 
ps30TESTinit<-matrix(NA, nrow = 100 , ncol = 2)
ps30TESTinit[,2]<-as.vector(ps30_all_train_test[101:200,2])
ps30TESTinit[,1]<-seq(1:100)

ps30TESTinitlow<-ps30TESTinit[,1][which(ps30TESTinit[,2]>ps30seperation)]  
ps30TESTinitbig<-ps30TESTinit[,1][which(ps30TESTinit[,2]<ps30seperation)]  


#############

DATAphi<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_phi.csv", header = TRUE, row.names = 1)

x<-DATAphi$Initial.value.of.Sensor.1[DATAphi$Classification..MFPC.Scores.=="Low"]
y<-DATAphi$Initial.value.of.Sensor.1[DATAphi$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(F-G)
which((F-G)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((F-G)==J)
c=gridxy[k[5]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

phiseperation<-c

phi_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/6-phitrain_test_all.csv", header = TRUE,row.names = 1)
phi_all_train_test
phi_all_train_test<-as.matrix(phi_all_train_test)

phiTESTinit<-matrix(NA, nrow = 100 , ncol = 2)
phiTESTinit[,2]<-as.vector(phi_all_train_test[101:200,2])
phiTESTinit[,1]<-seq(1:100)

phiTESTinitlow<-phiTESTinit[,1][which(phiTESTinit[,2]<phiseperation)]
phiTESTinitbig<-phiTESTinit[,1][which(phiTESTinit[,2]>phiseperation)] 

#############

DATABPR<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_BPR.csv", header = TRUE, row.names = 1)

x<-DATABPR$Initial.value.of.Sensor.1[DATABPR$Classification..MFPC.Scores.=="Low"]
y<-DATABPR$Initial.value.of.Sensor.1[DATABPR$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(G-F)
which((G-F)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((G-F)==J)
c=gridxy[k[6]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

BPRseperation<-c

BPR_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/7-BPRtrain_test_all.csv", header = TRUE,row.names = 1)
BPR_all_train_test
BPR_all_train_test<-as.matrix(BPR_all_train_test)
 
BPRTESTinit<-matrix(NA, nrow = 100 , ncol = 2)
BPRTESTinit[,2]<-as.vector(BPR_all_train_test[101:200,2])
BPRTESTinit[,1]<-seq(1:100)

BPRTESTinitlow<-BPRTESTinit[,1][which(BPRTESTinit[,2]>BPRseperation)]  
BPRTESTinitbig<-BPRTESTinit[,1][which(BPRTESTinit[,2]<BPRseperation)]  

#############

DATAW31<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_W31.csv", header = TRUE, row.names = 1)

x<-DATAW31$Initial.value.of.Sensor.1[DATAW31$Classification..MFPC.Scores.=="Low"]
y<-DATAW31$Initial.value.of.Sensor.1[DATAW31$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(F-G)
which((F-G)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((F-G)==J)
c=gridxy[k[5]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

W31seperation<-c

W31_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/8-W31train_test_all.csv", header = TRUE,row.names = 1)
W31_all_train_test
W31_all_train_test<-as.matrix(W31_all_train_test)
 
W31TESTinit<-matrix(NA, nrow = 100 , ncol = 2)
W31TESTinit[,2]<-as.vector(W31_all_train_test[101:200,2])
W31TESTinit[,1]<-seq(1:100)

W31TESTinitlow<-W31TESTinit[,1][which(W31TESTinit[,2]<W31seperation)] 
W31TESTinitbig<-W31TESTinit[,1][which(W31TESTinit[,2]>W31seperation)] 

#############
#############
#############
#############

DATAW32<-read.csv("C:/DrCy/erre/data_for_registration/initials/train initials_new_W32.csv", header = TRUE, row.names = 1)

x<-DATAW32$Initial.value.of.Sensor.1[DATAW32$Classification..MFPC.Scores.=="Low"]
y<-DATAW32$Initial.value.of.Sensor.1[DATAW32$Classification..MFPC.Scores.=="Big"]

#Compute the empirical distribution function for both separately in a grid:
gridxy=seq(min(x,y),max(x,y),length.out=1000)
F=numeric()
G=numeric()
for (i in 1:1000){
  F[i]=mean(x<=gridxy[i])
  G[i]=mean(y<=gridxy[i])
}

#The Youden index is the maximum difference between these two functions: J=max()
plot(gridxy,F,type="l",ylim=c(0,1))
lines(gridxy,G,col="red")
J=max(F-G)
which((F-G)==J)

#This last sentence gives you the optimal cutoff point. It may be ties (not probable but not impossible). In such a case, we can just take any of the cases. Let's define then:
k=which((F-G)==J)
c=gridxy[k[5]]
c 

#So those curves whose initial point is below c will be classified as red (see **) and black otherwise
abline(v=c,col="blue",lty=3)

p1=hist(x)
p2=hist(y)

plot(p1, col=rgb(0,0,1,1/4), xlim=c(min(gridxy),max(gridxy)))  # first histogram
plot(p2, col=rgb(1,0,0,1/4), xlim=c(min(gridxy),max(gridxy)), add=T)
abline(v=c,col="blue",lty=3)

W32seperation<-c

W32_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/9-W32train_test_all.csv", header = TRUE,row.names = 1)
W32_all_train_test
W32_all_train_test<-as.matrix(W32_all_train_test)
 
W32TESTinit<-matrix(NA, nrow = 100 , ncol = 2)
W32TESTinit[,2]<-as.vector(W32_all_train_test[101:200,2])
W32TESTinit[,1]<-seq(1:100)

W32TESTinitlow<-W32TESTinit[,1][which(W32TESTinit[,2]<W32seperation)]  
W32TESTinitbig<-W32TESTinit[,1][which(W32TESTinit[,2]>W32seperation)]  


T24TESTinitlow
T30TESTinitlow
T50TESTinitlow
P30TESTinitlow
ps30TESTinitlow
phiTESTinitlow
BPRTESTinitlow
W31TESTinitlow
W32TESTinitlow

T24TESTinitbig
T30TESTinitbig
T50TESTinitbig
P30TESTinitbig
ps30TESTinitbig
phiTESTinitbig
BPRTESTinitbig
W31TESTinitbig
W32TESTinitbig

allclassified<- matrix(NA,  ncol= 9, nrow = 100)
###############

for (i in 1:100) 
{
  
  if(T24TESTinit[i,2]>T24seperation)
  {allclassified[i,1]<- "low"}
  else{allclassified[i,1]<- "big"
  }
  if(T30TESTinit[i,2]>T30seperation)
  {allclassified[i,2]<- "low"}
  else{allclassified[i,2]<- "big"
  }
  if(T50TESTinit[i,2]>T50seperation)
  {allclassified[i,3]<- "low"}
  else{allclassified[i,3]<- "big"
  }
  if(P30TESTinit[i,2]<P30seperation)
  {allclassified[i,4]<- "low"}
  else{allclassified[i,4]<- "big"
  }
  if(ps30TESTinit[i,2]>ps30seperation)
  {allclassified[i,5]<- "low"}
  else{allclassified[i,5]<- "big"
  }
  if(phiTESTinit[i,2]<phiseperation)
  {allclassified[i,6]<- "low"}
  else{allclassified[i,6]<- "big"
  }
  if(BPRTESTinit[i,2]>BPRseperation)
  {allclassified[i,7]<- "low"}
  else{allclassified[i,7]<- "big"
  }
  if(W31TESTinit[i,2]<W31seperation)
  {allclassified[i,8]<- "low"}
  else{allclassified[i,8]<- "big"
  }
  if(W32TESTinit[i,2]<W32seperation)
  {allclassified[i,9]<- "low"}
  else{allclassified[i,9]<- "big"
  }
}


TESTclassified<- matrix(NA,  ncol= 2, nrow = 100)
TESTclassified[,1]<-seq(1:100)
for (i in 1:100) {
  if(length(allclassified[i,][which(allclassified[i,]=="low")])>4)
  {TESTclassified[i,2]<-"low"}
  else{TESTclassified[i,2]<-"big"}
}

TestLowClassEng<-as.numeric(TESTclassified[,1][which(TESTclassified[,2]=="low")])
TestBigClassEng<-as.numeric(TESTclassified[,1][which(TESTclassified[,2]=="big")])

######
TESTclassifieddeneme<- matrix(NA,  ncol= 2, nrow = 100)
TESTclassifieddeneme[,1]<-seq(1:100)
for (i in 1:100) {
  if(length(allclassified[i,][which(allclassified[i,]=="low")])>5)
  {TESTclassifieddeneme[i,2]<-"low"}
  else if(length(allclassified[i,][which(allclassified[i,]=="big")])>5)
  {TESTclassifieddeneme[i,2]<-"big"}
  else {TESTclassifieddeneme[i,2]<-"5 yada 4"}
}
length(TestLowClassEng)
length(TestBigClassEng)


##############TEST DATA CLASSIFICATION DONE HERE (YOUDEN INDEX)  ################


# -----------------------------------------------------------------------------
# 8. CALCULATE DISTANCES FOR DIFFERENT GROUPS
# -----------------------------------------------------------------------------

###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION#########################


T24_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/1-T24train_test_all.csv", header = TRUE,row.names = 1)
T24_all_train_test
T24_all_train_test<-as.matrix(T24_all_train_test)


#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  T24_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(T24_all_train_test[100+i,])))
  list_test_matrix[[i]]<-T24_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
T24TrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
T24TrainScoresmatrix<-matrix(NA,100,3)
T24TrainScoresmatrix[,1]<-seq(1:100)
T24TrainScoresmatrix[,2]<-T24TrainScores[,1]
T24TrainScoresmatrix[,3]<-T24TrainScores[,2]

TrainLowScores<-as.numeric(T24TrainScoresmatrix[,1][which(T24TrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(T24TrainScoresmatrix[,1][which(T24TrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_T24_SCORE<-list()
for (i in TestLowClassEng)
{ 
  T24_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(T24_all_train_test[100+i,])))
  list_test_matrix_T24_SCORE[[i]]<-T24_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  T24_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(T24_all_train_test[100+i,])))
  list_test_matrix_T24_SCORE[[i]]<-T24_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  T24_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_T24_SCORE[[i]])[1,]))
  for (j in 1:66) {
    T24_test_matrix[j,]<-T24_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_T24_SCORE[[i]])[1,])]
    list_test_matrix_T24_SCORE[[i]][j,]<-T24_test_matrix[j,]
  }
  list_test_matrix_T24_SCORE[[i]][66,]<-T24_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  T24_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_T24_SCORE[[i]])[1,]))
  for (j in 1:36) {
    T24_test_matrix[j,]<-T24_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_T24_SCORE[[i]])[1,])]
    list_test_matrix_T24_SCORE[[i]][j,]<-T24_test_matrix[j,]
  }
  list_test_matrix_T24_SCORE[[i]][36,]<-T24_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_T24_SCORE[[4]]
test3<-list_test_matrix_T24_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_T24_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_T24 <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_T24_SCORE[[i]][66,])))
  noNA_test_matrix_T24 <-na.omit(list_test_matrix_T24_SCORE[[i]])
  list_test_matrix_SCORE_T24_noNA[[i]]<-noNA_test_matrix_T24
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_T24 <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_T24_SCORE[[i]][36,])))
  noNA_test_matrix_T24 <-na.omit(list_test_matrix_T24_SCORE[[i]])
  list_test_matrix_SCORE_T24_noNA[[i]]<-noNA_test_matrix_T24
}
#####list with transpose of each matrix

#list_test_matrix_t_SCORE_T24_noNA<-list()
#for (i in 1:100)
#{ 
#  list_test_matrix_t_SCORE_T24_noNA[[i]]<-t(list_test_matrix_SCORE_T24_noNA[[i]])
#}

list_test_matrix_t_SCORE_T24_noNA<-list()
for (i in TestLowClassEng)
{ 
  list_test_matrix_t_SCORE_T24_noNA[[i]]<-t(list_test_matrix_SCORE_T24_noNA[[i]])
}
for (i in TestBigClassEng)
{ 
  list_test_matrix_t_SCORE_T24_noNA[[i]]<-t(list_test_matrix_SCORE_T24_noNA[[i]])
}


######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_T24_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_T24_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T24_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_T24_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_T24_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_T24_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_T24_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_T24_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_T24_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_T24_Big_scores[[j]]<-list_test_smoothing_T24_Big_scores
}


list_test_all_smooth_T24_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_T24_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T24_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_T24_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_T24_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_T24_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_T24_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_T24_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_T24_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_T24_Low_scores[[j]]<-list_test_smoothing_T24_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=7

#lows cutted
plot(list_test_all_smooth_T24_Low_scores[[testengine]][[1]],xlim=c(0,330), ylim=c(641.7,644))
for (i in 1:nrow(list_test_matrix_SCORE_T24_noNA[[testengine]]))
{
  lines(list_test_all_smooth_T24_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_T24_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T24_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_T24_Big_scores[[testengine]][[1]],xlim=c(0,330), ylim=c(641.7,644))
for (i in 1:nrow(list_test_matrix_SCORE_T24_noNA[[testengine]]))
{
  lines(list_test_all_smooth_T24_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_T24_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T24_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_T24_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T24_noNA[[testengine]])]],xlim=c(0,370), ylim=c(641.7,644))
for (i in as.vector(list_test_matrix_SCORE_T24_noNA[[testengine]][,1])) {
  lines(listeT24[[i]])
} 
lines(list_test_all_smooth_T24_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T24_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_T24_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T24_noNA[[testengine]])]],xlim=c(0,370), ylim=c(641.7,644))
for (i in as.vector(list_test_matrix_SCORE_T24_noNA[[testengine]][,1])) {
  lines(listeT24[[i]])
} 
lines(list_test_all_smooth_T24_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T24_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_T24_noNA[[j]])

FUNDATAT24LIST<-list()
for (j in TestLowClassEng) {
  FUNDATAT24LIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T24_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_T24_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_T24_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAT24LIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATAT24LIST[[j]]<-FUNDATAT24LIST_LIST
}

FUNDATAT24LIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATAT24LIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T24_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_T24_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_T24_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAT24LIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATAT24LIST_BIG[[j]]<-FUNDATAT24LIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_T24_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_T24_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(T24_all_train_test)

vectorT24<-as.vector(T24_all_train_test[101:200,2:363])
vectorT24<-na.omit(vectorT24)
mean(na.omit(vectorT24))
meanT24<-mean(na.omit(vectorT24))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_T24_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixT24reg<-matrix(NA,nrow(list_test_matrix_SCORE_T24_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_T24_noNA[[j]])) {
    distancematrixT24reg[i,]<-FUNDATAT24LIST_BIG[[j]][[i]]@X/meanT24
  }
  list_for_distance_T24_reg_big_scores[[j]]<-distancematrixT24reg
}

list_for_distance_T24_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixT24reg<-matrix(NA,nrow(list_test_matrix_SCORE_T24_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_T24_noNA[[j]])) {
    distancematrixT24reg[i,]<-FUNDATAT24LIST[[j]][[i]]@X/meanT24
  }
  list_for_distance_T24_reg_low_scores[[j]]<-distancematrixT24reg
}

############# DISTANCE CALCUALTION YAPILDI
library(philentropy)
list_distance_T24_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_T24_reg_low_scores[[i]], method = "euclidean" )
  list_distance_T24_reg_low_scores[[i]]<-DISTANCE
}

list_distance_T24_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_T24_reg_big_scores[[i]], method = "euclidean" )
  list_distance_T24_reg_big_scores[[i]]<-DISTANCE
}


###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T24 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################


###################T30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################T30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################T30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################T30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
T30_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/2-T30train_test_all.csv", header = TRUE,row.names = 1)
T30_all_train_test
T30_all_train_test<-as.matrix(T30_all_train_test)


#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  T30_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(T30_all_train_test[100+i,])))
  list_test_matrix[[i]]<-T30_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
T30TrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
T30TrainScoresmatrix<-matrix(NA,100,3)
T30TrainScoresmatrix[,1]<-seq(1:100)
T30TrainScoresmatrix[,2]<-T30TrainScores[,1]
T30TrainScoresmatrix[,3]<-T30TrainScores[,2]

TrainLowScores<-as.numeric(T30TrainScoresmatrix[,1][which(T30TrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(T30TrainScoresmatrix[,1][which(T30TrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_T30_SCORE<-list()
for (i in TestLowClassEng)
{ 
  T30_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(T30_all_train_test[100+i,])))
  list_test_matrix_T30_SCORE[[i]]<-T30_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  T30_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(T30_all_train_test[100+i,])))
  list_test_matrix_T30_SCORE[[i]]<-T30_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  T30_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_T30_SCORE[[i]])[1,]))
  for (j in 1:66) {
    T30_test_matrix[j,]<-T30_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_T30_SCORE[[i]])[1,])]
    list_test_matrix_T30_SCORE[[i]][j,]<-T30_test_matrix[j,]
  }
  list_test_matrix_T30_SCORE[[i]][66,]<-T30_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  T30_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_T30_SCORE[[i]])[1,]))
  for (j in 1:36) {
    T30_test_matrix[j,]<-T30_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_T30_SCORE[[i]])[1,])]
    list_test_matrix_T30_SCORE[[i]][j,]<-T30_test_matrix[j,]
  }
  list_test_matrix_T30_SCORE[[i]][36,]<-T30_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_T30_SCORE[[4]]
test3<-list_test_matrix_T30_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_T30_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_T30 <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_T30_SCORE[[i]][66,])))
  noNA_test_matrix_T30 <-na.omit(list_test_matrix_T30_SCORE[[i]])
  list_test_matrix_SCORE_T30_noNA[[i]]<-noNA_test_matrix_T30
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_T30 <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_T30_SCORE[[i]][36,])))
  noNA_test_matrix_T30 <-na.omit(list_test_matrix_T30_SCORE[[i]])
  list_test_matrix_SCORE_T30_noNA[[i]]<-noNA_test_matrix_T30
}
#####list with transpose of each matrix

#list_test_matrix_t_SCORE_T30_noNA<-list()
#for (i in 1:100)
#{ 
#  list_test_matrix_t_SCORE_T30_noNA[[i]]<-t(list_test_matrix_SCORE_T30_noNA[[i]])
#}

list_test_matrix_t_SCORE_T30_noNA<-list()
for (i in TestLowClassEng)
{ 
  list_test_matrix_t_SCORE_T30_noNA[[i]]<-t(list_test_matrix_SCORE_T30_noNA[[i]])
}
for (i in TestBigClassEng)
{ 
  list_test_matrix_t_SCORE_T30_noNA[[i]]<-t(list_test_matrix_SCORE_T30_noNA[[i]])
}


######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_T30_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_T30_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T30_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_T30_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_T30_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_T30_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_T30_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_T30_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_T30_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_T30_Big_scores[[j]]<-list_test_smoothing_T30_Big_scores
}


list_test_all_smooth_T30_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_T30_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T30_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_T30_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_T30_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_T30_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_T30_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_T30_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_T30_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_T30_Low_scores[[j]]<-list_test_smoothing_T30_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=7

#lows cutted
plot(list_test_all_smooth_T30_Low_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(1580,1610))
for (i in 1:nrow(list_test_matrix_SCORE_T30_noNA[[testengine]]))
{
  lines(list_test_all_smooth_T30_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_T30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T30_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_T30_Big_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(1580,1610))
for (i in 1:nrow(list_test_matrix_SCORE_T30_noNA[[testengine]]))
{
  lines(list_test_all_smooth_T30_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_T30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T30_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_T30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T30_noNA[[testengine]])]],xlim=c(0,370), ylim=c(1580,1610))
for (i in as.vector(list_test_matrix_SCORE_T30_noNA[[testengine]][,1])) {
  lines(listeT30[[i]])
} 
lines(list_test_all_smooth_T30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T30_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_T30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T30_noNA[[testengine]])]],xlim=c(0,370), ylim=c(1580,1610))
for (i in as.vector(list_test_matrix_SCORE_T30_noNA[[testengine]][,1])) {
  lines(listeT30[[i]])
} 
lines(list_test_all_smooth_T30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T30_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_T30_noNA[[j]])

FUNDATAT30LIST<-list()
for (j in TestLowClassEng) {
  FUNDATAT30LIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T30_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_T30_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_T30_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAT30LIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATAT30LIST[[j]]<-FUNDATAT30LIST_LIST
}

FUNDATAT30LIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATAT30LIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T30_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_T30_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_T30_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAT30LIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATAT30LIST_BIG[[j]]<-FUNDATAT30LIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_T30_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_T30_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(T30_all_train_test)

vectorT30<-as.vector(T30_all_train_test[101:200,2:363])
vectorT30<-na.omit(vectorT30)
mean(na.omit(vectorT30))
meanT30<-mean(na.omit(vectorT30))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_T30_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixT30reg<-matrix(NA,nrow(list_test_matrix_SCORE_T30_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_T30_noNA[[j]])) {
    distancematrixT30reg[i,]<-FUNDATAT30LIST_BIG[[j]][[i]]@X/meanT30
  }
  list_for_distance_T30_reg_big_scores[[j]]<-distancematrixT30reg
}

list_for_distance_T30_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixT30reg<-matrix(NA,nrow(list_test_matrix_SCORE_T30_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_T30_noNA[[j]])) {
    distancematrixT30reg[i,]<-FUNDATAT30LIST[[j]][[i]]@X/meanT30
  }
  list_for_distance_T30_reg_low_scores[[j]]<-distancematrixT30reg
}

############# DISTANCE CALCUALTION YAPILDI
list_distance_T30_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_T30_reg_low_scores[[i]], method = "euclidean" )
  list_distance_T30_reg_low_scores[[i]]<-DISTANCE
}

list_distance_T30_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_T30_reg_big_scores[[i]], method = "euclidean" )
  list_distance_T30_reg_big_scores[[i]]<-DISTANCE
}

###################T30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################



###################T50 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################T50 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################T50 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################T50 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################

T50_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/3-T50train_test_all.csv", header = TRUE,row.names = 1)
T50_all_train_test
T50_all_train_test<-as.matrix(T50_all_train_test)

#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  T50_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(T50_all_train_test[100+i,])))
  list_test_matrix[[i]]<-T50_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
T50TrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
T50TrainScoresmatrix<-matrix(NA,100,3)
T50TrainScoresmatrix[,1]<-seq(1:100)
T50TrainScoresmatrix[,2]<-T50TrainScores[,1]
T50TrainScoresmatrix[,3]<-T50TrainScores[,2]

TrainLowScores<-as.numeric(T50TrainScoresmatrix[,1][which(T50TrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(T50TrainScoresmatrix[,1][which(T50TrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_T50_SCORE<-list()
for (i in TestLowClassEng)
{ 
  T50_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(T50_all_train_test[100+i,])))
  list_test_matrix_T50_SCORE[[i]]<-T50_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  T50_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(T50_all_train_test[100+i,])))
  list_test_matrix_T50_SCORE[[i]]<-T50_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  T50_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_T50_SCORE[[i]])[1,]))
  for (j in 1:66) {
    T50_test_matrix[j,]<-T50_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_T50_SCORE[[i]])[1,])]
    list_test_matrix_T50_SCORE[[i]][j,]<-T50_test_matrix[j,]
  }
  list_test_matrix_T50_SCORE[[i]][66,]<-T50_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  T50_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_T50_SCORE[[i]])[1,]))
  for (j in 1:36) {
    T50_test_matrix[j,]<-T50_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_T50_SCORE[[i]])[1,])]
    list_test_matrix_T50_SCORE[[i]][j,]<-T50_test_matrix[j,]
  }
  list_test_matrix_T50_SCORE[[i]][36,]<-T50_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_T50_SCORE[[4]]
test3<-list_test_matrix_T50_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_T50_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_T50 <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_T50_SCORE[[i]][66,])))
  noNA_test_matrix_T50 <-na.omit(list_test_matrix_T50_SCORE[[i]])
  list_test_matrix_SCORE_T50_noNA[[i]]<-noNA_test_matrix_T50
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_T50 <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_T50_SCORE[[i]][36,])))
  noNA_test_matrix_T50 <-na.omit(list_test_matrix_T50_SCORE[[i]])
  list_test_matrix_SCORE_T50_noNA[[i]]<-noNA_test_matrix_T50
}
#####list with transpose of each matrix

#list_test_matrix_t_SCORE_T50_noNA<-list()
#for (i in 1:100)
#{ 
#  list_test_matrix_t_SCORE_T50_noNA[[i]]<-t(list_test_matrix_SCORE_T50_noNA[[i]])
#}

list_test_matrix_t_SCORE_T50_noNA<-list()
for (i in TestLowClassEng)
{ 
  list_test_matrix_t_SCORE_T50_noNA[[i]]<-t(list_test_matrix_SCORE_T50_noNA[[i]])
}
for (i in TestBigClassEng)
{ 
  list_test_matrix_t_SCORE_T50_noNA[[i]]<-t(list_test_matrix_SCORE_T50_noNA[[i]])
}




######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_T50_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_T50_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T50_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_T50_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_T50_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_T50_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_T50_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_T50_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_T50_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_T50_Big_scores[[j]]<-list_test_smoothing_T50_Big_scores
}


list_test_all_smooth_T50_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_T50_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T50_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_T50_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_T50_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_T50_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_T50_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_T50_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_T50_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_T50_Low_scores[[j]]<-list_test_smoothing_T50_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=7

#lows cutted
plot(list_test_all_smooth_T50_Low_scores[[testengine]][[1]],xlim=c(0,330), ylim=c(1390,1435))
for (i in 1:nrow(list_test_matrix_SCORE_T50_noNA[[testengine]]))
{
  lines(list_test_all_smooth_T50_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_T50_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T50_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_T50_Big_scores[[testengine]][[1]],xlim=c(0,330), ylim=c(1390,1435))
for (i in 1:nrow(list_test_matrix_SCORE_T50_noNA[[testengine]]))
{
  lines(list_test_all_smooth_T50_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_T50_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T50_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_T50_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T50_noNA[[testengine]])]],xlim=c(0,330), ylim=c(1390,1435))
for (i in as.vector(list_test_matrix_SCORE_T50_noNA[[testengine]][,1])) {
  lines(listeT50[[i]])
} 
lines(list_test_all_smooth_T50_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_T50_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_T50_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T50_noNA[[testengine]])]],xlim=c(0,330), ylim=c(1390,1435))
for (i in as.vector(list_test_matrix_SCORE_T50_noNA[[testengine]][,1])) {
  lines(listeT50[[i]])
} 
lines(list_test_all_smooth_T50_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_T50_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_T50_noNA[[j]])

FUNDATAT50LIST<-list()
for (j in TestLowClassEng) {
  FUNDATAT50LIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T50_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_T50_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_T50_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAT50LIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATAT50LIST[[j]]<-FUNDATAT50LIST_LIST
}

FUNDATAT50LIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATAT50LIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_T50_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_T50_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_T50_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAT50LIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATAT50LIST_BIG[[j]]<-FUNDATAT50LIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_T50_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_T50_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(T50_all_train_test)

vectorT50<-as.vector(T50_all_train_test[101:200,2:363])
vectorT50<-na.omit(vectorT50)
mean(na.omit(vectorT50))
meanT50<-mean(na.omit(vectorT50))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_T50_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixT50reg<-matrix(NA,nrow(list_test_matrix_SCORE_T50_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_T50_noNA[[j]])) {
    distancematrixT50reg[i,]<-FUNDATAT50LIST_BIG[[j]][[i]]@X/meanT50
  }
  list_for_distance_T50_reg_big_scores[[j]]<-distancematrixT50reg
}

list_for_distance_T50_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixT50reg<-matrix(NA,nrow(list_test_matrix_SCORE_T50_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_T50_noNA[[j]])) {
    distancematrixT50reg[i,]<-FUNDATAT50LIST[[j]][[i]]@X/meanT50
  }
  list_for_distance_T50_reg_low_scores[[j]]<-distancematrixT50reg
}

############# DISTANCE CALCUALTION YAPILDI
list_distance_T50_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_T50_reg_low_scores[[i]], method = "euclidean" )
  list_distance_T50_reg_low_scores[[i]]<-DISTANCE
}

list_distance_T50_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_T50_reg_big_scores[[i]], method = "euclidean" )
  list_distance_T50_reg_big_scores[[i]]<-DISTANCE
}


###################T50 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T50 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T50 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################T50 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################


###################P30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################P30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################P30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################P30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################


P30_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/4-P30train_test_all.csv", header = TRUE,row.names = 1)
P30_all_train_test
P30_all_train_test<-as.matrix(P30_all_train_test)



#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  P30_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(P30_all_train_test[100+i,])))
  list_test_matrix[[i]]<-P30_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
P30TrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
P30TrainScoresmatrix<-matrix(NA,100,3)
P30TrainScoresmatrix[,1]<-seq(1:100)
P30TrainScoresmatrix[,2]<-P30TrainScores[,1]
P30TrainScoresmatrix[,3]<-P30TrainScores[,2]

TrainLowScores<-as.numeric(P30TrainScoresmatrix[,1][which(P30TrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(P30TrainScoresmatrix[,1][which(P30TrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_P30_SCORE<-list()
for (i in TestLowClassEng)
{ 
  P30_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(P30_all_train_test[100+i,])))
  list_test_matrix_P30_SCORE[[i]]<-P30_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  P30_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(P30_all_train_test[100+i,])))
  list_test_matrix_P30_SCORE[[i]]<-P30_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  P30_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_P30_SCORE[[i]])[1,]))
  for (j in 1:66) {
    P30_test_matrix[j,]<-P30_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_P30_SCORE[[i]])[1,])]
    list_test_matrix_P30_SCORE[[i]][j,]<-P30_test_matrix[j,]
  }
  list_test_matrix_P30_SCORE[[i]][66,]<-P30_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  P30_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_P30_SCORE[[i]])[1,]))
  for (j in 1:36) {
    P30_test_matrix[j,]<-P30_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_P30_SCORE[[i]])[1,])]
    list_test_matrix_P30_SCORE[[i]][j,]<-P30_test_matrix[j,]
  }
  list_test_matrix_P30_SCORE[[i]][36,]<-P30_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_P30_SCORE[[4]]
test3<-list_test_matrix_P30_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_P30_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_P30 <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_P30_SCORE[[i]][66,])))
  noNA_test_matrix_P30 <-na.omit(list_test_matrix_P30_SCORE[[i]])
  list_test_matrix_SCORE_P30_noNA[[i]]<-noNA_test_matrix_P30
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_P30 <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_P30_SCORE[[i]][36,])))
  noNA_test_matrix_P30 <-na.omit(list_test_matrix_P30_SCORE[[i]])
  list_test_matrix_SCORE_P30_noNA[[i]]<-noNA_test_matrix_P30
}
#####list with transpose of each matrix

#list_test_matrix_t_SCORE_P30_noNA<-list()
#for (i in 1:100)
#{ 
#  list_test_matrix_t_SCORE_P30_noNA[[i]]<-t(list_test_matrix_SCORE_P30_noNA[[i]])
#}

list_test_matrix_t_SCORE_P30_noNA<-list()
for (i in TestLowClassEng)
{ 
  list_test_matrix_t_SCORE_P30_noNA[[i]]<-t(list_test_matrix_SCORE_P30_noNA[[i]])
}
for (i in TestBigClassEng)
{ 
  list_test_matrix_t_SCORE_P30_noNA[[i]]<-t(list_test_matrix_SCORE_P30_noNA[[i]])
}




######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_P30_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_P30_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_P30_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_P30_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_P30_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_P30_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_P30_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_P30_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_P30_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_P30_Big_scores[[j]]<-list_test_smoothing_P30_Big_scores
}


list_test_all_smooth_P30_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_P30_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_P30_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_P30_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_P30_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_P30_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_P30_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_P30_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_P30_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_P30_Low_scores[[j]]<-list_test_smoothing_P30_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=7

#lows cutted
plot(list_test_all_smooth_P30_Low_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(550.5,555.2))
for (i in 1:nrow(list_test_matrix_SCORE_P30_noNA[[testengine]]))
{
  lines(list_test_all_smooth_P30_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_P30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_P30_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_P30_Big_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(550.5,555.2))
for (i in 1:nrow(list_test_matrix_SCORE_P30_noNA[[testengine]]))
{
  lines(list_test_all_smooth_P30_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_P30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_P30_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_P30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_P30_noNA[[testengine]])]],xlim=c(0,370), ylim=c(550.5,555.2))
for (i in as.vector(list_test_matrix_SCORE_P30_noNA[[testengine]][,1])) {
  lines(listeP30[[i]])
} 
lines(list_test_all_smooth_P30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_P30_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_P30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_P30_noNA[[testengine]])]],xlim=c(0,370), ylim=c(550.5,555.2))
for (i in as.vector(list_test_matrix_SCORE_P30_noNA[[testengine]][,1])) {
  lines(listeP30[[i]])
} 
lines(list_test_all_smooth_P30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_P30_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_P30_noNA[[j]])

FUNDATAP30LIST<-list()
for (j in TestLowClassEng) {
  FUNDATAP30LIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_P30_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_P30_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_P30_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAP30LIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATAP30LIST[[j]]<-FUNDATAP30LIST_LIST
}

FUNDATAP30LIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATAP30LIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_P30_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_P30_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_P30_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAP30LIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATAP30LIST_BIG[[j]]<-FUNDATAP30LIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_P30_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_P30_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(P30_all_train_test)

vectorP30<-as.vector(P30_all_train_test[101:200,2:363])
vectorP30<-na.omit(vectorP30)
mean(na.omit(vectorP30))
meanP30<-mean(na.omit(vectorP30))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_P30_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixP30reg<-matrix(NA,nrow(list_test_matrix_SCORE_P30_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_P30_noNA[[j]])) {
    distancematrixP30reg[i,]<-FUNDATAP30LIST_BIG[[j]][[i]]@X/meanP30
  }
  list_for_distance_P30_reg_big_scores[[j]]<-distancematrixP30reg
}

list_for_distance_P30_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixP30reg<-matrix(NA,nrow(list_test_matrix_SCORE_P30_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_P30_noNA[[j]])) {
    distancematrixP30reg[i,]<-FUNDATAP30LIST[[j]][[i]]@X/meanP30
  }
  list_for_distance_P30_reg_low_scores[[j]]<-distancematrixP30reg
}

############# DISTANCE CALCUALTION YAPILDI
list_distance_P30_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_P30_reg_low_scores[[i]], method = "euclidean" )
  list_distance_P30_reg_low_scores[[i]]<-DISTANCE
}

list_distance_P30_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_P30_reg_big_scores[[i]], method = "euclidean" )
  list_distance_P30_reg_big_scores[[i]]<-DISTANCE
}


###################P30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################P30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################P30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################P30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################

###################ps30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################ps30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################ps30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################ps30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################

# 
ps30_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/5-ps30train_test_all.csv", header = TRUE,row.names = 1)
ps30_all_train_test
ps30_all_train_test<-as.matrix(ps30_all_train_test)

#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  ps30_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(ps30_all_train_test[100+i,])))
  list_test_matrix[[i]]<-ps30_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
ps30TrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
ps30TrainScoresmatrix<-matrix(NA,100,3)
ps30TrainScoresmatrix[,1]<-seq(1:100)
ps30TrainScoresmatrix[,2]<-ps30TrainScores[,1]
ps30TrainScoresmatrix[,3]<-ps30TrainScores[,2]

TrainLowScores<-as.numeric(ps30TrainScoresmatrix[,1][which(ps30TrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(ps30TrainScoresmatrix[,1][which(ps30TrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_ps30_SCORE<-list()
for (i in TestLowClassEng)
{ 
  ps30_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(ps30_all_train_test[100+i,])))
  list_test_matrix_ps30_SCORE[[i]]<-ps30_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  ps30_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(ps30_all_train_test[100+i,])))
  list_test_matrix_ps30_SCORE[[i]]<-ps30_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  ps30_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_ps30_SCORE[[i]])[1,]))
  for (j in 1:66) {
    ps30_test_matrix[j,]<-ps30_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_ps30_SCORE[[i]])[1,])]
    list_test_matrix_ps30_SCORE[[i]][j,]<-ps30_test_matrix[j,]
  }
  list_test_matrix_ps30_SCORE[[i]][66,]<-ps30_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  ps30_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_ps30_SCORE[[i]])[1,]))
  for (j in 1:36) {
    ps30_test_matrix[j,]<-ps30_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_ps30_SCORE[[i]])[1,])]
    list_test_matrix_ps30_SCORE[[i]][j,]<-ps30_test_matrix[j,]
  }
  list_test_matrix_ps30_SCORE[[i]][36,]<-ps30_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_ps30_SCORE[[4]]
test3<-list_test_matrix_ps30_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_ps30_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_ps30 <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_ps30_SCORE[[i]][66,])))
  noNA_test_matrix_ps30 <-na.omit(list_test_matrix_ps30_SCORE[[i]])
  list_test_matrix_SCORE_ps30_noNA[[i]]<-noNA_test_matrix_ps30
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_ps30 <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_ps30_SCORE[[i]][36,])))
  noNA_test_matrix_ps30 <-na.omit(list_test_matrix_ps30_SCORE[[i]])
  list_test_matrix_SCORE_ps30_noNA[[i]]<-noNA_test_matrix_ps30
}
#####list with transpose of each matrix

#list_test_matrix_t_SCORE_ps30_noNA<-list()
#for (i in 1:100)
#{ 
#  list_test_matrix_t_SCORE_ps30_noNA[[i]]<-t(list_test_matrix_SCORE_ps30_noNA[[i]])
#}

list_test_matrix_t_SCORE_ps30_noNA<-list()
for (i in TestLowClassEng)
{ 
  list_test_matrix_t_SCORE_ps30_noNA[[i]]<-t(list_test_matrix_SCORE_ps30_noNA[[i]])
}
for (i in TestBigClassEng)
{ 
  list_test_matrix_t_SCORE_ps30_noNA[[i]]<-t(list_test_matrix_SCORE_ps30_noNA[[i]])
}




######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_ps30_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_ps30_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_ps30_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_ps30_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_ps30_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_ps30_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_ps30_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_ps30_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_ps30_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_ps30_Big_scores[[j]]<-list_test_smoothing_ps30_Big_scores
}


list_test_all_smooth_ps30_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_ps30_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_ps30_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_ps30_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_ps30_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_ps30_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_ps30_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_ps30_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_ps30_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_ps30_Low_scores[[j]]<-list_test_smoothing_ps30_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=8

#lows cutted
plot(list_test_all_smooth_ps30_Low_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(47.05,48.4))
for (i in 1:nrow(list_test_matrix_SCORE_ps30_noNA[[testengine]]))
{
  lines(list_test_all_smooth_ps30_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_ps30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_ps30_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_ps30_Big_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(47.05,48.4))
for (i in 1:nrow(list_test_matrix_SCORE_ps30_noNA[[testengine]]))
{
  lines(list_test_all_smooth_ps30_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_ps30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_ps30_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_ps30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_ps30_noNA[[testengine]])]],xlim=c(0,370), ylim=c(47.05,48.4))
for (i in as.vector(list_test_matrix_SCORE_ps30_noNA[[testengine]][,1])) {
  lines(listeps30[[i]])
} 
lines(list_test_all_smooth_ps30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_ps30_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_ps30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_ps30_noNA[[testengine]])]],xlim=c(0,370), ylim=c(47.05,48.4))
for (i in as.vector(list_test_matrix_SCORE_ps30_noNA[[testengine]][,1])) {
  lines(listeps30[[i]])
} 
lines(list_test_all_smooth_ps30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_ps30_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_ps30_noNA[[j]])

FUNDATAps30LIST<-list()
for (j in TestLowClassEng) {
  FUNDATAps30LIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_ps30_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_ps30_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_ps30_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAps30LIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATAps30LIST[[j]]<-FUNDATAps30LIST_LIST
}

FUNDATAps30LIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATAps30LIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_ps30_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_ps30_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_ps30_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAps30LIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATAps30LIST_BIG[[j]]<-FUNDATAps30LIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_ps30_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_ps30_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(ps30_all_train_test)

vectorps30<-as.vector(ps30_all_train_test[101:200,2:363])
vectorps30<-na.omit(vectorps30)
mean(na.omit(vectorps30))
meanps30<-mean(na.omit(vectorps30))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_ps30_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixps30reg<-matrix(NA,nrow(list_test_matrix_SCORE_ps30_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_ps30_noNA[[j]])) {
    distancematrixps30reg[i,]<-FUNDATAps30LIST_BIG[[j]][[i]]@X/meanps30
  }
  list_for_distance_ps30_reg_big_scores[[j]]<-distancematrixps30reg
}

list_for_distance_ps30_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixps30reg<-matrix(NA,nrow(list_test_matrix_SCORE_ps30_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_ps30_noNA[[j]])) {
    distancematrixps30reg[i,]<-FUNDATAps30LIST[[j]][[i]]@X/meanps30
  }
  list_for_distance_ps30_reg_low_scores[[j]]<-distancematrixps30reg
}

############# DISTANCE CALCUALTION YAPILDI
list_distance_ps30_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_ps30_reg_low_scores[[i]], method = "euclidean" )
  list_distance_ps30_reg_low_scores[[i]]<-DISTANCE
}

list_distance_ps30_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_ps30_reg_big_scores[[i]], method = "euclidean" )
  list_distance_ps30_reg_big_scores[[i]]<-DISTANCE
}



###################ps30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################ps30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################ps30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################ps30 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################

###################phi DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################phi DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################phi DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################phi DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################

# 
phi_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/6-phitrain_test_all.csv", header = TRUE,row.names = 1)
phi_all_train_test
phi_all_train_test<-as.matrix(phi_all_train_test)

#### list for 100 test engine - all 101 x number of observation matrix

#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  phi_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(phi_all_train_test[100+i,])))
  list_test_matrix[[i]]<-phi_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
phiTrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
phiTrainScoresmatrix<-matrix(NA,100,3)
phiTrainScoresmatrix[,1]<-seq(1:100)
phiTrainScoresmatrix[,2]<-phiTrainScores[,1]
phiTrainScoresmatrix[,3]<-phiTrainScores[,2]

TrainLowScores<-as.numeric(phiTrainScoresmatrix[,1][which(phiTrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(phiTrainScoresmatrix[,1][which(phiTrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_phi_SCORE<-list()
for (i in TestLowClassEng)
{ 
  phi_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(phi_all_train_test[100+i,])))
  list_test_matrix_phi_SCORE[[i]]<-phi_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  phi_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(phi_all_train_test[100+i,])))
  list_test_matrix_phi_SCORE[[i]]<-phi_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  phi_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_phi_SCORE[[i]])[1,]))
  for (j in 1:66) {
    phi_test_matrix[j,]<-phi_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_phi_SCORE[[i]])[1,])]
    list_test_matrix_phi_SCORE[[i]][j,]<-phi_test_matrix[j,]
  }
  list_test_matrix_phi_SCORE[[i]][66,]<-phi_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  phi_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_phi_SCORE[[i]])[1,]))
  for (j in 1:36) {
    phi_test_matrix[j,]<-phi_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_phi_SCORE[[i]])[1,])]
    list_test_matrix_phi_SCORE[[i]][j,]<-phi_test_matrix[j,]
  }
  list_test_matrix_phi_SCORE[[i]][36,]<-phi_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_phi_SCORE[[4]]
test3<-list_test_matrix_phi_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_phi_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_phi <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_phi_SCORE[[i]][66,])))
  noNA_test_matrix_phi <-na.omit(list_test_matrix_phi_SCORE[[i]])
  list_test_matrix_SCORE_phi_noNA[[i]]<-noNA_test_matrix_phi
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_phi <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_phi_SCORE[[i]][36,])))
  noNA_test_matrix_phi <-na.omit(list_test_matrix_phi_SCORE[[i]])
  list_test_matrix_SCORE_phi_noNA[[i]]<-noNA_test_matrix_phi
}
#####list with transpose of each matrix

list_test_matrix_t_SCORE_phi_noNA<-list()
for (i in 1:100)
{ 
  list_test_matrix_t_SCORE_phi_noNA[[i]]<-t(list_test_matrix_SCORE_phi_noNA[[i]])
}



######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_phi_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_phi_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_phi_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_phi_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_phi_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_phi_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_phi_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_phi_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_phi_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_phi_Big_scores[[j]]<-list_test_smoothing_phi_Big_scores
}


list_test_all_smooth_phi_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_phi_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_phi_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_phi_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_phi_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_phi_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_phi_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_phi_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_phi_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_phi_Low_scores[[j]]<-list_test_smoothing_phi_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=8

#lows cutted
plot(list_test_all_smooth_phi_Low_scores[[testengine]][[1]],xlim=c(0,330), ylim=c(519.1,522.8))
for (i in 1:nrow(list_test_matrix_SCORE_phi_noNA[[testengine]]))
{
  lines(list_test_all_smooth_phi_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_phi_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_phi_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_phi_Big_scores[[testengine]][[1]],xlim=c(0,330), ylim=c(519.1,522.8))
for (i in 1:nrow(list_test_matrix_SCORE_phi_noNA[[testengine]]))
{
  lines(list_test_all_smooth_phi_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_phi_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_phi_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_phi_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_phi_noNA[[testengine]])]],xlim=c(0,330), ylim=c(519.1,522.8))
for (i in as.vector(list_test_matrix_SCORE_phi_noNA[[testengine]][,1])) {
  lines(listephi[[i]])
} 
lines(list_test_all_smooth_phi_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_phi_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_phi_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_phi_noNA[[testengine]])]],xlim=c(0,330), ylim=c(519.1,522.8))
for (i in as.vector(list_test_matrix_SCORE_phi_noNA[[testengine]][,1])) {
  lines(listephi[[i]])
} 
lines(list_test_all_smooth_phi_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_phi_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_phi_noNA[[j]])

FUNDATAphiLIST<-list()
for (j in TestLowClassEng) {
  FUNDATAphiLIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_phi_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_phi_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_phi_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAphiLIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATAphiLIST[[j]]<-FUNDATAphiLIST_LIST
}

FUNDATAphiLIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATAphiLIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_phi_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_phi_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_phi_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAphiLIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATAphiLIST_BIG[[j]]<-FUNDATAphiLIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_phi_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_phi_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(phi_all_train_test)

vectorphi<-as.vector(phi_all_train_test[101:200,2:363])
vectorphi<-na.omit(vectorphi)
mean(na.omit(vectorphi))
meanphi<-mean(na.omit(vectorphi))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_phi_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixphireg<-matrix(NA,nrow(list_test_matrix_SCORE_phi_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_phi_noNA[[j]])) {
    distancematrixphireg[i,]<-FUNDATAphiLIST_BIG[[j]][[i]]@X/meanphi
  }
  list_for_distance_phi_reg_big_scores[[j]]<-distancematrixphireg
}

list_for_distance_phi_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixphireg<-matrix(NA,nrow(list_test_matrix_SCORE_phi_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_phi_noNA[[j]])) {
    distancematrixphireg[i,]<-FUNDATAphiLIST[[j]][[i]]@X/meanphi
  }
  list_for_distance_phi_reg_low_scores[[j]]<-distancematrixphireg
}

############# DISTANCE CALCUALTION YAPILDI
list_distance_phi_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_phi_reg_low_scores[[i]], method = "euclidean" )
  list_distance_phi_reg_low_scores[[i]]<-DISTANCE
}

list_distance_phi_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_phi_reg_big_scores[[i]], method = "euclidean" )
  list_distance_phi_reg_big_scores[[i]]<-DISTANCE
}


###################phi DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################phi DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################phi DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################phi DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################

###################BPR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################BPR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################BPR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################BPR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################

BPR_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/7-BPRtrain_test_all.csv", header = TRUE,row.names = 1)
BPR_all_train_test
BPR_all_train_test<-as.matrix(BPR_all_train_test)


#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  BPR_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(BPR_all_train_test[100+i,])))
  list_test_matrix[[i]]<-BPR_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
BPRTrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
BPRTrainScoresmatrix<-matrix(NA,100,3)
BPRTrainScoresmatrix[,1]<-seq(1:100)
BPRTrainScoresmatrix[,2]<-BPRTrainScores[,1]
BPRTrainScoresmatrix[,3]<-BPRTrainScores[,2]

TrainLowScores<-as.numeric(BPRTrainScoresmatrix[,1][which(BPRTrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(BPRTrainScoresmatrix[,1][which(BPRTrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_BPR_SCORE<-list()
for (i in TestLowClassEng)
{ 
  BPR_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(BPR_all_train_test[100+i,])))
  list_test_matrix_BPR_SCORE[[i]]<-BPR_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  BPR_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(BPR_all_train_test[100+i,])))
  list_test_matrix_BPR_SCORE[[i]]<-BPR_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  BPR_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_BPR_SCORE[[i]])[1,]))
  for (j in 1:66) {
    BPR_test_matrix[j,]<-BPR_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_BPR_SCORE[[i]])[1,])]
    list_test_matrix_BPR_SCORE[[i]][j,]<-BPR_test_matrix[j,]
  }
  list_test_matrix_BPR_SCORE[[i]][66,]<-BPR_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  BPR_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_BPR_SCORE[[i]])[1,]))
  for (j in 1:36) {
    BPR_test_matrix[j,]<-BPR_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_BPR_SCORE[[i]])[1,])]
    list_test_matrix_BPR_SCORE[[i]][j,]<-BPR_test_matrix[j,]
  }
  list_test_matrix_BPR_SCORE[[i]][36,]<-BPR_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_BPR_SCORE[[4]]
test3<-list_test_matrix_BPR_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_BPR_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_BPR <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_BPR_SCORE[[i]][66,])))
  noNA_test_matrix_BPR <-na.omit(list_test_matrix_BPR_SCORE[[i]])
  list_test_matrix_SCORE_BPR_noNA[[i]]<-noNA_test_matrix_BPR
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_BPR <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_BPR_SCORE[[i]][36,])))
  noNA_test_matrix_BPR <-na.omit(list_test_matrix_BPR_SCORE[[i]])
  list_test_matrix_SCORE_BPR_noNA[[i]]<-noNA_test_matrix_BPR
}
#####list with transpose of each matrix

#list_test_matrix_t_SCORE_BPR_noNA<-list()
#for (i in 1:100)
#{ 
#  list_test_matrix_t_SCORE_BPR_noNA[[i]]<-t(list_test_matrix_SCORE_BPR_noNA[[i]])
#}

list_test_matrix_t_SCORE_BPR_noNA<-list()
for (i in TestLowClassEng)
{ 
  list_test_matrix_t_SCORE_BPR_noNA[[i]]<-t(list_test_matrix_SCORE_BPR_noNA[[i]])
}
for (i in TestBigClassEng)
{ 
  list_test_matrix_t_SCORE_BPR_noNA[[i]]<-t(list_test_matrix_SCORE_BPR_noNA[[i]])
}



######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_BPR_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_BPR_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_BPR_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_BPR_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_BPR_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_BPR_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_BPR_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_BPR_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_BPR_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_BPR_Big_scores[[j]]<-list_test_smoothing_BPR_Big_scores
}


list_test_all_smooth_BPR_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_BPR_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_BPR_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_BPR_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_BPR_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_BPR_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_BPR_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_BPR_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_BPR_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_BPR_Low_scores[[j]]<-list_test_smoothing_BPR_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=8

#lows cutted
plot(list_test_all_smooth_BPR_Low_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(8.37,8.55))
for (i in 1:nrow(list_test_matrix_SCORE_BPR_noNA[[testengine]]))
{
  lines(list_test_all_smooth_BPR_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_BPR_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_BPR_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_BPR_Big_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(8.37,8.55))
for (i in 1:nrow(list_test_matrix_SCORE_BPR_noNA[[testengine]]))
{
  lines(list_test_all_smooth_BPR_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_BPR_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_BPR_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_BPR_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_BPR_noNA[[testengine]])]],xlim=c(0,370), ylim=c(8.37,8.55))
for (i in as.vector(list_test_matrix_SCORE_BPR_noNA[[testengine]][,1])) {
  lines(listeBPR[[i]])
} 
lines(list_test_all_smooth_BPR_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_BPR_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_BPR_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_BPR_noNA[[testengine]])]],xlim=c(0,370), ylim=c(8.37,8.55))
for (i in as.vector(list_test_matrix_SCORE_BPR_noNA[[testengine]][,1])) {
  lines(listeBPR[[i]])
} 
lines(list_test_all_smooth_BPR_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_BPR_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_BPR_noNA[[j]])

FUNDATABPRLIST<-list()
for (j in TestLowClassEng) {
  FUNDATABPRLIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_BPR_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_BPR_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_BPR_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATABPRLIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATABPRLIST[[j]]<-FUNDATABPRLIST_LIST
}

FUNDATABPRLIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATABPRLIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_BPR_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_BPR_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_BPR_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATABPRLIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATABPRLIST_BIG[[j]]<-FUNDATABPRLIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_BPR_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_BPR_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(BPR_all_train_test)

vectorBPR<-as.vector(BPR_all_train_test[101:200,2:363])
vectorBPR<-na.omit(vectorBPR)
mean(na.omit(vectorBPR))
meanBPR<-mean(na.omit(vectorBPR))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_BPR_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixBPRreg<-matrix(NA,nrow(list_test_matrix_SCORE_BPR_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_BPR_noNA[[j]])) {
    distancematrixBPRreg[i,]<-FUNDATABPRLIST_BIG[[j]][[i]]@X/meanBPR
  }
  list_for_distance_BPR_reg_big_scores[[j]]<-distancematrixBPRreg
}

list_for_distance_BPR_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixBPRreg<-matrix(NA,nrow(list_test_matrix_SCORE_BPR_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_BPR_noNA[[j]])) {
    distancematrixBPRreg[i,]<-FUNDATABPRLIST[[j]][[i]]@X/meanBPR
  }
  list_for_distance_BPR_reg_low_scores[[j]]<-distancematrixBPRreg
}

############# DISTANCE CALCUALTION YAPILDI
list_distance_BPR_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_BPR_reg_low_scores[[i]], method = "euclidean" )
  list_distance_BPR_reg_low_scores[[i]]<-DISTANCE
}

list_distance_BPR_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_BPR_reg_big_scores[[i]], method = "euclidean" )
  list_distance_BPR_reg_big_scores[[i]]<-DISTANCE
}


###################BPR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################BPR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################BPR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################BPR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################

###################W31 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################W31 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################W31 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################W31 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################


W31_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/8-W31train_test_all.csv", header = TRUE,row.names = 1)
W31_all_train_test
W31_all_train_test<-as.matrix(W31_all_train_test)


#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  W31_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(W31_all_train_test[100+i,])))
  list_test_matrix[[i]]<-W31_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
W31TrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
W31TrainScoresmatrix<-matrix(NA,100,3)
W31TrainScoresmatrix[,1]<-seq(1:100)
W31TrainScoresmatrix[,2]<-W31TrainScores[,1]
W31TrainScoresmatrix[,3]<-W31TrainScores[,2]

TrainLowScores<-as.numeric(W31TrainScoresmatrix[,1][which(W31TrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(W31TrainScoresmatrix[,1][which(W31TrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_W31_SCORE<-list()
for (i in TestLowClassEng)
{ 
  W31_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(W31_all_train_test[100+i,])))
  list_test_matrix_W31_SCORE[[i]]<-W31_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  W31_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(W31_all_train_test[100+i,])))
  list_test_matrix_W31_SCORE[[i]]<-W31_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  W31_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_W31_SCORE[[i]])[1,]))
  for (j in 1:66) {
    W31_test_matrix[j,]<-W31_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_W31_SCORE[[i]])[1,])]
    list_test_matrix_W31_SCORE[[i]][j,]<-W31_test_matrix[j,]
  }
  list_test_matrix_W31_SCORE[[i]][66,]<-W31_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  W31_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_W31_SCORE[[i]])[1,]))
  for (j in 1:36) {
    W31_test_matrix[j,]<-W31_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_W31_SCORE[[i]])[1,])]
    list_test_matrix_W31_SCORE[[i]][j,]<-W31_test_matrix[j,]
  }
  list_test_matrix_W31_SCORE[[i]][36,]<-W31_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_W31_SCORE[[4]]
test3<-list_test_matrix_W31_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_W31_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_W31 <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_W31_SCORE[[i]][66,])))
  noNA_test_matrix_W31 <-na.omit(list_test_matrix_W31_SCORE[[i]])
  list_test_matrix_SCORE_W31_noNA[[i]]<-noNA_test_matrix_W31
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_W31 <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_W31_SCORE[[i]][36,])))
  noNA_test_matrix_W31 <-na.omit(list_test_matrix_W31_SCORE[[i]])
  list_test_matrix_SCORE_W31_noNA[[i]]<-noNA_test_matrix_W31
}
#####list with transpose of each matrix

#list_test_matrix_t_SCORE_W31_noNA<-list()
#for (i in 1:100)
#{ 
#  list_test_matrix_t_SCORE_W31_noNA[[i]]<-t(list_test_matrix_SCORE_W31_noNA[[i]])
#}

list_test_matrix_t_SCORE_W31_noNA<-list()
for (i in TestLowClassEng)
{ 
  list_test_matrix_t_SCORE_W31_noNA[[i]]<-t(list_test_matrix_SCORE_W31_noNA[[i]])
}
for (i in TestBigClassEng)
{ 
  list_test_matrix_t_SCORE_W31_noNA[[i]]<-t(list_test_matrix_SCORE_W31_noNA[[i]])
}



######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_W31_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_W31_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_W31_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_W31_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_W31_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_W31_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_W31_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_W31_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_W31_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_W31_Big_scores[[j]]<-list_test_smoothing_W31_Big_scores
}


list_test_all_smooth_W31_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_W31_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_W31_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_W31_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_W31_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_W31_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_W31_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_W31_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_W31_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_W31_Low_scores[[j]]<-list_test_smoothing_W31_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=1

#lows cutted
plot(list_test_all_smooth_W31_Low_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(38.2,39.2))
for (i in 1:nrow(list_test_matrix_SCORE_W31_noNA[[testengine]]))
{
  lines(list_test_all_smooth_W31_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_W31_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W31_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_W31_Big_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(38.2,39.2))
for (i in 1:nrow(list_test_matrix_SCORE_W31_noNA[[testengine]]))
{
  lines(list_test_all_smooth_W31_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_W31_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W31_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_W31_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W31_noNA[[testengine]])]],xlim=c(0,370), ylim=c(38.2,39.2))
for (i in as.vector(list_test_matrix_SCORE_W31_noNA[[testengine]][,1])) {
  lines(listeW31[[i]])
} 
lines(list_test_all_smooth_W31_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W31_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_W31_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W31_noNA[[testengine]])]],xlim=c(0,370), ylim=c(38.2,39.2))
for (i in as.vector(list_test_matrix_SCORE_W31_noNA[[testengine]][,1])) {
  lines(listeW31[[i]])
} 
lines(list_test_all_smooth_W31_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W31_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_W31_noNA[[j]])

FUNDATAW31LIST<-list()
for (j in TestLowClassEng) {
  FUNDATAW31LIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_W31_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_W31_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_W31_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAW31LIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATAW31LIST[[j]]<-FUNDATAW31LIST_LIST
}

FUNDATAW31LIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATAW31LIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_W31_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_W31_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_W31_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAW31LIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATAW31LIST_BIG[[j]]<-FUNDATAW31LIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_W31_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_W31_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(W31_all_train_test)

vectorW31<-as.vector(W31_all_train_test[101:200,2:363])
vectorW31<-na.omit(vectorW31)
mean(na.omit(vectorW31))
meanW31<-mean(na.omit(vectorW31))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_W31_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixW31reg<-matrix(NA,nrow(list_test_matrix_SCORE_W31_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_W31_noNA[[j]])) {
    distancematrixW31reg[i,]<-FUNDATAW31LIST_BIG[[j]][[i]]@X/meanW31
  }
  list_for_distance_W31_reg_big_scores[[j]]<-distancematrixW31reg
}

list_for_distance_W31_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixW31reg<-matrix(NA,nrow(list_test_matrix_SCORE_W31_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_W31_noNA[[j]])) {
    distancematrixW31reg[i,]<-FUNDATAW31LIST[[j]][[i]]@X/meanW31
  }
  list_for_distance_W31_reg_low_scores[[j]]<-distancematrixW31reg
}

############# DISTANCE CALCUALTION YAPILDI
list_distance_W31_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_W31_reg_low_scores[[i]], method = "euclidean" )
  list_distance_W31_reg_low_scores[[i]]<-DISTANCE
}

list_distance_W31_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_W31_reg_big_scores[[i]], method = "euclidean" )
  list_distance_W31_reg_big_scores[[i]]<-DISTANCE
}



###################W31 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################W31 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################W31 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################W31 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################

###################W32 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################W32 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################W32 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################W32 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################


W32_all_train_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/9-W32train_test_all.csv", header = TRUE,row.names = 1)
W32_all_train_test
W32_all_train_test<-as.matrix(W32_all_train_test)


#### list for 100 test engine - all 101 x number of observation matrix
list_test_matrix<-list()
for (i in 1:100)
{ 
  W32_test_matrix_length <- matrix(data = NA, nrow=101, ncol=length(na.omit(W32_all_train_test[100+i,])))
  list_test_matrix[[i]]<-W32_test_matrix_length
}

####train datasında low ve high scoreları ayırdım.
W32TrainScores<-read.csv("C:/DrCy/erre/data_for_registration/train mfpc scores_with low-big.csv", header = TRUE,row.names = 1)
W32TrainScoresmatrix<-matrix(NA,100,3)
W32TrainScoresmatrix[,1]<-seq(1:100)
W32TrainScoresmatrix[,2]<-W32TrainScores[,1]
W32TrainScoresmatrix[,3]<-W32TrainScores[,2]

TrainLowScores<-as.numeric(W32TrainScoresmatrix[,1][which(W32TrainScoresmatrix[,3]=="Low")])
TrainBigScores<-as.numeric(W32TrainScoresmatrix[,1][which(W32TrainScoresmatrix[,3]=="Big")])

TestLowClassEng
TestBigClassEng
list_test_matrix_W32_SCORE<-list()
for (i in TestLowClassEng)
{ 
  W32_test_matrix_length <- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(na.omit(W32_all_train_test[100+i,])))
  list_test_matrix_W32_SCORE[[i]]<-W32_test_matrix_length
}

for (i in TestBigClassEng)
{ 
  W32_test_matrix_length <- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(na.omit(W32_all_train_test[100+i,])))
  list_test_matrix_W32_SCORE[[i]]<-W32_test_matrix_length
}

#### fill all lists with 100 train + row101 for related  test engine
for (i in TestLowClassEng) { 
  W32_test_matrix<- matrix(data = NA, nrow=(length(TrainLowScores)+1), ncol=length(as.matrix(list_test_matrix_W32_SCORE[[i]])[1,]))
  for (j in 1:66) {
    W32_test_matrix[j,]<-W32_all_train_test[TrainLowScores[j],1:length(as.matrix(list_test_matrix_W32_SCORE[[i]])[1,])]
    list_test_matrix_W32_SCORE[[i]][j,]<-W32_test_matrix[j,]
  }
  list_test_matrix_W32_SCORE[[i]][66,]<-W32_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

for (i in TestBigClassEng) { 
  W32_test_matrix<- matrix(data = NA, nrow=(length(TrainBigScores)+1), ncol=length(as.matrix(list_test_matrix_W32_SCORE[[i]])[1,]))
  for (j in 1:36) {
    W32_test_matrix[j,]<-W32_all_train_test[TrainBigScores[j],1:length(as.matrix(list_test_matrix_W32_SCORE[[i]])[1,])]
    list_test_matrix_W32_SCORE[[i]][j,]<-W32_test_matrix[j,]
  }
  list_test_matrix_W32_SCORE[[i]][36,]<-W32_all_train_test[100+i, 1:length(as.matrix(list_test_matrix[[i]])[1,])]
}

test3<-list_test_matrix_W32_SCORE[[4]]
test3<-list_test_matrix_W32_SCORE[[5]]


length(TestLowClassEng)
#### ignore NAs for all
#####list with ignoring train data with observation less than test data
list_test_matrix_SCORE_W32_noNA<-list()
for (i in TestLowClassEng)
{
  noNA_test_matrix_W32 <- matrix(data = NA, nrow=66, ncol=length(na.omit(list_test_matrix_W32_SCORE[[i]][66,])))
  noNA_test_matrix_W32 <-na.omit(list_test_matrix_W32_SCORE[[i]])
  list_test_matrix_SCORE_W32_noNA[[i]]<-noNA_test_matrix_W32
}
for (i in TestBigClassEng)
{ 
  noNA_test_matrix_W32 <- matrix(data = NA, nrow=36, ncol=length(na.omit(list_test_matrix_W32_SCORE[[i]][36,])))
  noNA_test_matrix_W32 <-na.omit(list_test_matrix_W32_SCORE[[i]])
  list_test_matrix_SCORE_W32_noNA[[i]]<-noNA_test_matrix_W32
}
#####list with transpose of each matrix

#list_test_matrix_t_SCORE_W32_noNA<-list()
#for (i in 1:100)
#{ 
#  list_test_matrix_t_SCORE_W32_noNA[[i]]<-t(list_test_matrix_SCORE_W32_noNA[[i]])
#}

list_test_matrix_t_SCORE_W32_noNA<-list()
for (i in TestLowClassEng)
{ 
  list_test_matrix_t_SCORE_W32_noNA[[i]]<-t(list_test_matrix_SCORE_W32_noNA[[i]])
}
for (i in TestBigClassEng)
{ 
  list_test_matrix_t_SCORE_W32_noNA[[i]]<-t(list_test_matrix_SCORE_W32_noNA[[i]])
}



######smoothing for each engine
no_of_splines=5
j=10
i=2
list_test_all_smooth_W32_Big_scores<-list()
for (j in TestBigClassEng) {
  list_test_smoothing_W32_Big_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_W32_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_W32_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_W32_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_W32_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_W32_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_W32_Big_scores[[i]]<-Smooth
  }
  list_test_all_smooth_W32_Big_scores[[j]]<-list_test_smoothing_W32_Big_scores
}


list_test_all_smooth_W32_Low_scores<-list()
for (j in TestLowClassEng) {
  list_test_smoothing_W32_Low_scores<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[j]]))
  {
    Smooth<- smooth.basis( argvals = seq(1,length(na.omit(list_test_matrix_t_SCORE_W32_noNA[[j]][,1]))-1, length.out= length(na.omit(list_test_matrix_t_SCORE_W32_noNA[[j]][,1]))-1),
                           y= as.vector(list_test_matrix_SCORE_W32_noNA[[j]][i,2:length(na.omit(list_test_matrix_t_SCORE_W32_noNA[[j]][,1]))]), 
                           fdParobj = create.bspline.basis(c(1,length(na.omit(list_test_matrix_t_SCORE_W32_noNA[[j]][,1]))-1),no_of_splines))
    list_test_smoothing_W32_Low_scores[[i]]<-Smooth
  }
  list_test_all_smooth_W32_Low_scores[[j]]<-list_test_smoothing_W32_Low_scores
}



TestLowClassEng
TestBigClassEng
###burası 66 ve 36 lık smoothlar
###plot for cutted at observation
testengine=8

#lows cutted
plot(list_test_all_smooth_W32_Low_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(22.97,23.5))
for (i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[testengine]]))
{
  lines(list_test_all_smooth_W32_Low_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

#bigs cutted
plot(list_test_all_smooth_W32_Big_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(22.97,23.5))
for (i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[testengine]]))
{
  lines(list_test_all_smooth_W32_Big_scores[[testengine]][[i]], col="black")
}
lines(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")



###plot for (full observation trains - cutted test)

#lows full
plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(22.97,23.5))
for (i in as.vector(list_test_matrix_SCORE_W32_noNA[[testengine]][,1])) {
  lines(listeW32[[i]])
} 
lines(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

#bigs full
plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(22.97,23.5))
for (i in as.vector(list_test_matrix_SCORE_W32_noNA[[testengine]][,1])) {
  lines(listeW32[[i]])
} 
lines(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")



j
TestLowClassEng
TestBigClassEng

#################DISTANCE########################
#########################REGISTRATION OF SMOOTHED
#####Distance için 20 noktadan observation alınıyor

i=1
j=8

1:nrow(list_test_matrix_SCORE_W32_noNA[[j]])

FUNDATAW32LIST<-list()
for (j in TestLowClassEng) {
  FUNDATAW32LIST_LIST<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_W32_Low_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_W32_Low_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAW32LIST_LIST[[i]]<-funDatasmooth
  }
  FUNDATAW32LIST[[j]]<-FUNDATAW32LIST_LIST
}

FUNDATAW32LIST_BIG<-list()
for (j in TestBigClassEng) {
  FUNDATAW32LIST_LIST_BIG<-list()
  for (i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[j]])) {
    funDatasmooth<-fd2funData(list_test_all_smooth_W32_Big_scores[[j]][[i]]$fd,
                              argvals = seq(from=1, to=length(list_test_all_smooth_W32_Big_scores[[j]][[i]]$argvals), length.out= 20))
    FUNDATAW32LIST_LIST_BIG[[i]]<-funDatasmooth
  }
  FUNDATAW32LIST_BIG[[j]]<-FUNDATAW32LIST_LIST_BIG
}

#funDatasmooth<-fd2funData(list_test_all_smooth_W32_Big_scores[[j]][[1]]$fd,
#                          argvals = seq(from=1, to=length(list_test_all_smooth_W32_Big_scores[[j]][[1]]$argvals), length.out= 20))
#funDatasmooth@X


#####DISTANCE Calculation
##Bu observationlardan matrix oluşturuluyor.

dim(W32_all_train_test)

vectorW32<-as.vector(W32_all_train_test[101:200,2:363])
vectorW32<-na.omit(vectorW32)
mean(na.omit(vectorW32))
meanW32<-mean(na.omit(vectorW32))

#####BURADA DISTANCE HESABI ICIN MATRIXLER LISTEYE ALINDI
list_for_distance_W32_reg_big_scores<-list()
for (j in TestBigClassEng) {
  distancematrixW32reg<-matrix(NA,nrow(list_test_matrix_SCORE_W32_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[j]])) {
    distancematrixW32reg[i,]<-FUNDATAW32LIST_BIG[[j]][[i]]@X/meanW32
  }
  list_for_distance_W32_reg_big_scores[[j]]<-distancematrixW32reg
}

list_for_distance_W32_reg_low_scores<-list()
for (j in TestLowClassEng) {
  distancematrixW32reg<-matrix(NA,nrow(list_test_matrix_SCORE_W32_noNA[[j]]),20) 
  for(i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[j]])) {
    distancematrixW32reg[i,]<-FUNDATAW32LIST[[j]][[i]]@X/meanW32
  }
  list_for_distance_W32_reg_low_scores[[j]]<-distancematrixW32reg
}

############# DISTANCE CALCUALTION YAPILDI
list_distance_W32_reg_low_scores<-list()
for (i in TestLowClassEng) {
  DISTANCE<-distance(list_for_distance_W32_reg_low_scores[[i]], method = "euclidean" )
  list_distance_W32_reg_low_scores[[i]]<-DISTANCE
}

list_distance_W32_reg_big_scores<-list()
for (i in TestBigClassEng) {
  DISTANCE<-distance(list_for_distance_W32_reg_big_scores[[i]], method = "euclidean" )
  list_distance_W32_reg_big_scores[[i]]<-DISTANCE
}


###################W32 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################W32 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################W32 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################W32 DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################



# -----------------------------------------------------------------------------
# 9. CALCULATE MULTIVARIATE EUCLIDIAN DISTANCE FOR LOW AND HIGH SCORE GROUPS
# -----------------------------------------------------------------------------



###################MULTIVAR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################MULTIVAR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################MULTIVAR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
###################MULTIVAR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - START#########################
#Multivariate Distance calculation
list_all_distance_merged_low_scores<-list()
for (i in TestLowClassEng) {
  
  a<-list_distance_T24_reg_low_scores[[i]]+
    list_distance_T30_reg_low_scores[[i]]+
    list_distance_T50_reg_low_scores[[i]]+
    list_distance_P30_reg_low_scores[[i]]+
    list_distance_ps30_reg_low_scores[[i]]+
    list_distance_phi_reg_low_scores[[i]]+
    list_distance_BPR_reg_low_scores[[i]]+
    list_distance_W31_reg_low_scores[[i]]+
    list_distance_W32_reg_low_scores[[i]]
  
  list_all_distance_merged_low_scores[[i]]<-a
}

list_all_distance_merged_big_scores<-list()
for (i in TestBigClassEng) {
  
  a<-list_distance_T24_reg_big_scores[[i]]+
    list_distance_T30_reg_big_scores[[i]]+
    list_distance_T50_reg_big_scores[[i]]+
    list_distance_P30_reg_big_scores[[i]]+
    list_distance_ps30_reg_big_scores[[i]]+
    list_distance_phi_reg_big_scores[[i]]+
    list_distance_BPR_reg_big_scores[[i]]+
    list_distance_W31_reg_big_scores[[i]]+
    list_distance_W32_reg_big_scores[[i]]
  list_all_distance_merged_big_scores[[i]]<-a
}

############## SORTED DISTANCES FOR ALL ENGINES
list_dist_all_sorted<-list()
for (i in TestLowClassEng) {
  DISTANCE<-matrix(NA,nrow(list_test_matrix_SCORE_W31_noNA[[i]]),2)
  DISTANCE[,1]<-list_test_matrix_SCORE_W31_noNA[[i]][,1]
  DISTANCE[,2]<-list_all_distance_merged_low_scores[[i]][,nrow(list_test_matrix_SCORE_W31_noNA[[i]])]
  DISTANCE<-DISTANCE[order(DISTANCE[,2],decreasing = FALSE),]
  DISTANCE
  list_dist_all_sorted[[i]]<-DISTANCE
}

for (i in TestBigClassEng) {
  DISTANCE<-matrix(NA,nrow(list_test_matrix_SCORE_W31_noNA[[i]]),2)
  DISTANCE[,1]<-list_test_matrix_SCORE_W31_noNA[[i]][,1]
  DISTANCE[,2]<-list_all_distance_merged_big_scores[[i]][,nrow(list_test_matrix_SCORE_W31_noNA[[i]])]
  DISTANCE<-DISTANCE[order(DISTANCE[,2],decreasing = FALSE),]
  DISTANCE
  list_dist_all_sorted[[i]]<-DISTANCE
}



######################PLOTS

TestBigClassEng
TestLowClassEng
no_of_closest=5
testengine=93

#lows
par(mfrow=c(3,3))
par(cex.lab=cex, cex.axis=cex, cex.main=cex)
cex=2
plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(641.7,644), xlab="Cycle Time")
title(main="T24", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeT24[[i]], col="red")
}
lines(list_test_all_smooth_T24_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="blue")
text(70, 643, "Test Engine no:70", col="blue")
text(200, 642.5, "5 closest train engines (low score group)", col="red")

plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(1580,1610), xlab="Cycle Time")
title(main="T30", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeT30[[i]])
}
lines(list_test_all_smooth_T30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(1390,1435), xlab="Cycle Time")
title(main="T50", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeT50[[i]])
}
lines(list_test_all_smooth_T50_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(550.5,555.2), xlab="Cycle Time")
title(main="P30", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeP30[[i]])
}
lines(list_test_all_smooth_P30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(47.05,48.4), xlab="Cycle Time")
title(main="ps30", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeps30[[i]])
}
lines(list_test_all_smooth_ps30_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(519.1,522.8), xlab="Cycle Time")
title(main="phi", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listephi[[i]])
}
lines(list_test_all_smooth_phi_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(8.37,8.55), xlab="Cycle Time")
title(main="BPR", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeBPR[[i]])
}
lines(list_test_all_smooth_BPR_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")


plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(38.2,39.2), xlab="Cycle Time")
title(main="W31", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeW31[[i]])
}
lines(list_test_all_smooth_W31_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")


plot(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(22.97,23.5), xlab="Cycle Time")
title(main="W32", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeW32[[i]])
}
lines(list_test_all_smooth_W32_Low_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")


testengine=82
#Bigs
par(mfrow=c(3,3))
par(cex.lab=cex, cex.axis=cex, cex.main=cex)
cex=2
plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(641.7,644), xlab="Cycle Time")
title(main="T24", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeT24[[i]])
}
lines(list_test_all_smooth_T24_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="blue")
text(60, 642.5, "Test Engine no:81", col="blue")

text(220, 642.3, "5 closest train engines (big score group)", col="black")

plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(1580,1610), xlab="Cycle Time")
title(main="T30", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeT30[[i]])
}
lines(list_test_all_smooth_T30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="blue")

plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(1390,1435), xlab="Cycle Time")
title(main="T50", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeT50[[i]])
}
lines(list_test_all_smooth_T50_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(550.5,555.2), xlab="Cycle Time")
title(main="P30", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeP30[[i]])
}
lines(list_test_all_smooth_P30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(47.05,48.4), xlab="Cycle Time")
title(main="ps30", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeps30[[i]])
}
lines(list_test_all_smooth_ps30_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(519.1,522.8), xlab="Cycle Time")
title(main="phi", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listephi[[i]])
}
lines(list_test_all_smooth_phi_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")

plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(8.37,8.55), xlab="Cycle Time")
title(main="BPR", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeBPR[[i]])
}
lines(list_test_all_smooth_BPR_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")


plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(38.2,39.2), xlab="Cycle Time")
title(main="W31", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeW31[[i]])
}
lines(list_test_all_smooth_W31_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")


plot(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]],xlim=c(0,370), ylim=c(22.97,23.5), xlab="Cycle Time")
title(main="W32", cex.main=3, )
for (i in c(list_dist_all_sorted[[testengine]][2:(no_of_closest+1),1]))
{
  lines(listeW32[[i]])
}
lines(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")



#bigs cutted
#plot(list_test_all_smooth_W32_Big_scores[[testengine]][[1]],xlim=c(0,370), ylim=c(22.97,23.5))
#for (i in 1:nrow(list_test_matrix_SCORE_W32_noNA[[testengine]]))
#{
#  lines(list_test_all_smooth_W32_Big_scores[[testengine]][[i]], col="black")
#}
#lines(list_test_all_smooth_W32_Big_scores[[testengine]][[nrow(list_test_matrix_SCORE_W32_noNA[[testengine]])]], col="red")


###################MULTIVAR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################MULTIVAR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################MULTIVAR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################MULTIVAR DISTANCE CALCULATION AFTER MFPCA CLASSIFICATION - FINISH#########################



# -----------------------------------------------------------------------------
# 10. CALCULATE RMSEs
# -----------------------------------------------------------------------------


###################RUL PREDICTION AFTER MFPCA CLASSIFICATION - START#########################
###################RUL PREDICTION AFTER MFPCA CLASSIFICATION - START#########################
###################RUL PREDICTION AFTER MFPCA CLASSIFICATION - START#########################
###################RUL PREDICTION AFTER MFPCA CLASSIFICATION - START#########################

####TRUE RULimport
RUL_test<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/1-T24test_RUL.csv", header = TRUE)
RUL_test
RUL_test<-as.matrix(RUL_test)


par(mfrow=c(1,1))
par(cex.lab=cex, cex.axis=cex, cex.main=cex)
cex=2
no_of_closest_engine=8
TRUE_RUL_DECREASING<-RUL_test[order(RUL_test[,3],decreasing = FALSE),]
plot(TRUE_RUL_DECREASING[,3])

RUL_test_bigscore<- matrix(NA,ncol=4, nrow= length(TestBigClassEng))
for (i in 1:length(TestBigClassEng)) {
  RUL_test_bigscore[i,]<- RUL_test[TestBigClassEng[i],]
}

TRUE_RUL_DECREASING_bigscore<-RUL_test_bigscore[order(RUL_test_bigscore[,3],decreasing = FALSE),]
plot(TRUE_RUL_DECREASING_bigscore[,3])

RUL_test_lowscore<- matrix(NA,ncol=4, nrow= length(TestLowClassEng))
for (i in 1:length(TestLowClassEng)) {
  RUL_test_lowscore[i,]<- RUL_test[TestLowClassEng[i],]
}

TRUE_RUL_DECREASING_lowscore<-RUL_test_lowscore[order(RUL_test_lowscore[,3],decreasing = FALSE),]
plot(TRUE_RUL_DECREASING_lowscore[,3])


####RUL Prediction

par(mfrow=c(1,1))
par(cex.lab=cex, cex.axis=cex, cex.main=cex)
cex=2


LIFE_TRAIN<-read.csv("C:/DrCy/erre/data_for_registration/test/NEW/1-T24train_LIFE.csv", header = TRUE)
LIFE_TRAIN<-as.matrix(LIFE_TRAIN)
LIFE_TRAIN

no_of_closest_engine<-7
list_closest_X_lowscore<-list()
for (i in TestLowClassEng ) {
  close<-list_dist_all_sorted[[i]][2:(no_of_closest_engine+1),1]
  list_closest_X_lowscore[[i]]<-close
}

#19 eğer 75 lik yaptıysan
TestBigClassEngexcept49<-TestBigClassEng
TestBigClassEngexcept49[27]<-NA
TestBigClassEngexcept49<-na.omit(TestBigClassEngexcept49)

list_closest_X_bigcore<-list()
for (i in TestBigClassEngexcept49 ) {
  close<-list_dist_all_sorted[[i]][2:(no_of_closest_engine+1),1]
  list_closest_X_bigcore[[i]]<-close
}
list_closest_X_bigcore[[49]]<-c(96,69,92,67)



# FOR ENGINE 49 - 96,69,92,67
#1 NO OF ENGINE
#2 OBS Test
#3 RUL TEST
#4 LIFE TEST
#5 PRED LIFE
#6 PRED RUL
#7 MEAN RUL ERROR
#8 PRED LIFE MEDIAN
#9 PRED RUL MEDIAN
#9 MEDIAN RUL ERROR

RUL_PREDICTION_lowscore<- matrix(data = NA, nrow=length(TestLowClassEng), 10)
newheaderstest<- c("Engine No","No of OBS", "True RUL", "TRUE Life",
                   "PRED Life(Mean)","Pred RUL(Mean)","RUL ERROR(Mean)",
                   "PRED Life(Med)","Pred RUL(Med)","RUL ERROR(Med)")
colnames(RUL_PREDICTION_lowscore) <- newheaderstest
RUL_PREDICTION_lowscore[,1]<-RUL_test_lowscore[,1]
RUL_PREDICTION_lowscore[,2]<-RUL_test_lowscore[,2]
RUL_PREDICTION_lowscore[,3]<-RUL_test_lowscore[,3]
RUL_PREDICTION_lowscore[,4]<-RUL_test_lowscore[,4]
i=1
for (i in 1:length(TestLowClassEng)) {
  RUL_PREDICTION_lowscore[i,5]<-mean(LIFE_TRAIN[list_closest_X_lowscore[[TestLowClassEng[i]]],2])
  RUL_PREDICTION_lowscore[i,6]<-RUL_PREDICTION_lowscore[i,5]-RUL_PREDICTION_lowscore[i,2]
  RUL_PREDICTION_lowscore[i,7]<-RUL_PREDICTION_lowscore[i,6]-RUL_PREDICTION_lowscore[i,3]
  RUL_PREDICTION_lowscore[i,8]<-median(LIFE_TRAIN[list_closest_X_lowscore[[TestLowClassEng[i]]],2]) 
  RUL_PREDICTION_lowscore[i,9]<-RUL_PREDICTION_lowscore[i,8]-RUL_PREDICTION_lowscore[i,2]
  RUL_PREDICTION_lowscore[i,10]<-RUL_PREDICTION_lowscore[i,9]-RUL_PREDICTION_lowscore[i,3]
}

RUL_PREDICTION_bigscore<- matrix(data = NA, nrow=length(TestBigClassEng), 10)
newheaderstest<- c("Engine No","No of OBS", "True RUL", "TRUE Life",
                   "PRED Life(Mean)","Pred RUL(Mean)","RUL ERROR(Mean)",
                   "PRED Life(Med)","Pred RUL(Med)","RUL ERROR(Med)")
colnames(RUL_PREDICTION_bigscore) <- newheaderstest
RUL_PREDICTION_bigscore[,1]<-RUL_test_bigscore[,1]
RUL_PREDICTION_bigscore[,2]<-RUL_test_bigscore[,2]
RUL_PREDICTION_bigscore[,3]<-RUL_test_bigscore[,3]
RUL_PREDICTION_bigscore[,4]<-RUL_test_bigscore[,4]

for (i in 1:length(TestBigClassEng)) {
  RUL_PREDICTION_bigscore[i,5]<-mean(LIFE_TRAIN[list_closest_X_bigcore[[TestBigClassEng[i]]],2])
  RUL_PREDICTION_bigscore[i,6]<-RUL_PREDICTION_bigscore[i,5]-RUL_PREDICTION_bigscore[i,2]
  RUL_PREDICTION_bigscore[i,7]<-RUL_PREDICTION_bigscore[i,6]-RUL_PREDICTION_bigscore[i,3]
  RUL_PREDICTION_bigscore[i,8]<-median(LIFE_TRAIN[list_closest_X_bigcore[[TestBigClassEng[i]]],2]) 
  RUL_PREDICTION_bigscore[i,9]<-RUL_PREDICTION_bigscore[i,8]-RUL_PREDICTION_bigscore[i,2]
  RUL_PREDICTION_bigscore[i,10]<-RUL_PREDICTION_bigscore[i,9]-RUL_PREDICTION_bigscore[i,3]
}

TRUE_RUL_DECREASING_lowscore
ordervector_mean<- as.vector(TRUE_RUL_DECREASING_lowscore[,1]) # ORDER vector for increasing RUL
RUL_PREDICTION_lowscore_SORTED <- RUL_PREDICTION_lowscore[match(ordervector_mean, RUL_PREDICTION_lowscore),]##SORT increasing RUL

ordervector_mean<- as.vector(TRUE_RUL_DECREASING_bigscore[,1]) # ORDER vector for increasing RUL
RUL_PREDICTION_bigscore_SORTED <- RUL_PREDICTION_bigscore[match(ordervector_mean, RUL_PREDICTION_bigscore),]##SORT increasing RUL


######RMSE MEAN
ERRORSQUARE <- RUL_PREDICTION_lowscore_SORTED[,7]^2

SUM_ERRORSQUARE<-sum(ERRORSQUARE)
length(ERRORSQUARE)
SUM_ERRORSQUARE/length(ERRORSQUARE)
sqrt(SUM_ERRORSQUARE/length(ERRORSQUARE))


######RMSE MEDIAN
ERRORSQUAREMED <- RUL_PREDICTION_lowscore_SORTED[,10]^2
SUM_ERRORSQUAREMED<-sum(ERRORSQUAREMED)
SUM_ERRORSQUAREMED/length(ERRORSQUAREMED)
sqrt(SUM_ERRORSQUAREMED/length(ERRORSQUARE))

max(RUL_PREDICTION_lowscore[,7])
min(RUL_PREDICTION_lowscore[,7])
max(RUL_PREDICTION_lowscore[,10])
min(RUL_PREDICTION_lowscore[,10])

length(which(RUL_PREDICTION_lowscore[,10]< (-0.5)) )

length(which(RUL_PREDICTION_lowscore[,10]> (0.5)) )

length(which(RUL_PREDICTION_lowscore[,10]==(-0.5)) )
length(which(RUL_PREDICTION_lowscore[,10]==(0.5)) )
length(which(RUL_PREDICTION_lowscore[,10]==(0)) )

length(which(RUL_PREDICTION_lowscore[,7]< (-0.5)) )

length(which(RUL_PREDICTION_lowscore[,7]> (0.5)) )

length(which(RUL_PREDICTION_lowscore[,10]< (-0.5)) )

length(which(RUL_PREDICTION_lowscore[,10]> (0.5)) )

length(which(RUL_PREDICTION_lowscore[,10]==(-0.5)) )
length(which(RUL_PREDICTION_lowscore[,10]==(0.5)) )
length(which(RUL_PREDICTION_lowscore[,10]==(0)) )


######RMSE MEAN
ERRORSQUARE <- RUL_PREDICTION_bigscore_SORTED[,7]^2

SUM_ERRORSQUARE<-sum(ERRORSQUARE)
length(ERRORSQUARE)
SUM_ERRORSQUARE/length(ERRORSQUARE)
sqrt(SUM_ERRORSQUARE/length(ERRORSQUARE))


######RMSE MEDIAN
ERRORSQUAREMED <- RUL_PREDICTION_bigscore_SORTED[,10]^2
SUM_ERRORSQUAREMED<-sum(ERRORSQUAREMED)
SUM_ERRORSQUAREMED/length(ERRORSQUAREMED)
sqrt(SUM_ERRORSQUAREMED/length(ERRORSQUAREMED))

max(RUL_PREDICTION_bigscore[,7])
min(RUL_PREDICTION_bigscore[,7])
max(RUL_PREDICTION_bigscore[,10])
min(RUL_PREDICTION_bigscore[,10])

length(which(RUL_PREDICTION_bigscore[,10]< (-0.5)) )

length(which(RUL_PREDICTION_bigscore[,10]> (0.5)) )

length(which(RUL_PREDICTION_bigscore[,10]==(-0.5)) )
length(which(RUL_PREDICTION_bigscore[,10]==(0.5)) )
length(which(RUL_PREDICTION_bigscore[,10]==(0)) )

length(which(RUL_PREDICTION_bigscore[,7]< (-0.5)) )

length(which(RUL_PREDICTION_bigscore[,7]> (0.5)) )

length(which(RUL_PREDICTION_bigscore[,10]< (-0.5)) )

length(which(RUL_PREDICTION_bigscore[,10]> (0.5)) )

length(which(RUL_PREDICTION_bigscore[,10]==(-0.5)) )
length(which(RUL_PREDICTION_bigscore[,10]==(0.5)) )
length(which(RUL_PREDICTION_bigscore[,10]==(0)) )




RUL_PREDICTION_lowscore
RUL_PREDICTION_bigscore

total<-rbind(RUL_PREDICTION_lowscore,RUL_PREDICTION_bigscore)

ERRORSQUARE <- total[,7]^2

SUM_ERRORSQUARE<-sum(ERRORSQUARE)
SUM_ERRORSQUARE/100
sqrt(SUM_ERRORSQUARE/100)


######RMSE MEDIAN
ERRORSQUAREMED <- total[,10]^2
SUM_ERRORSQUAREMED<-sum(ERRORSQUAREMED)
SUM_ERRORSQUAREMED/100
sqrt(SUM_ERRORSQUAREMED/100)


###################RUL PREDICTION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################RUL PREDICTION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################RUL PREDICTION AFTER MFPCA CLASSIFICATION - FINISH#########################
###################RUL PREDICTION AFTER MFPCA CLASSIFICATION - FINISH#########################