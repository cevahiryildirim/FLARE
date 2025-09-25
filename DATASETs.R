# =============================================================================
# Health Prognostics in Multi-sensor Systems (FLARE Methodology)
#
# This script implements the entire analysis pipeline presented in the paper:
# "Health Prognostics in Multi-sensor Systems Based on Multivariate Functional
# Data Analysis"
#
# The script will:
# 1. INSTALL PACKAGES
# 2. LOAD LIBRARIES
# 3. SETUP AND DATA LOADING
# 4. CONTROL THE DATA
# =============================================================================


# -----------------------------------------------------------------------------
# 1. INSTALL PACKAGES
# -----------------------------------------------------------------------------

# Load required packages
install.packages("CMAPSS")


# -----------------------------------------------------------------------------
# 2. LOAD LIBRARIES
# -----------------------------------------------------------------------------

# Load required libraries
library(CMAPSS)

# -----------------------------------------------------------------------------
# 3. SETUP AND DATA LOADING
# -----------------------------------------------------------------------------

# Load dataset
data("CMAPSS")
str(CMAPSS$train)

# DATASET sizes
CMAPSS$subsets


# Split data by LIFE, OBS, RUL
FD001_TRAIN_LIFE<-CMAPSS$train$N[1:100]
FD002_TRAIN_LIFE<-CMAPSS$train$N[101:360]
FD003_TRAIN_LIFE<-CMAPSS$train$N[361:460]
FD004_TRAIN_LIFE<-CMAPSS$train$N[461:709]

FD001_TEST_OBS<-CMAPSS$test$N[1:100]
FD002_TEST_OBS<-CMAPSS$test$N[101:359]
FD003_TEST_OBS<-CMAPSS$test$N[360:459]
FD004_TEST_OBS<-CMAPSS$test$N[460:707]

FD001_TEST_RUL<-CMAPSS$test$RUL[1:100]
FD002_TEST_RUL<-CMAPSS$test$RUL[101:359]
FD003_TEST_RUL<-CMAPSS$test$RUL[360:459]
FD004_TEST_RUL<-CMAPSS$test$RUL[460:707]

FD001_TEST_LIFE<-FD001_TEST_OBS+FD001_TEST_RUL
FD002_TEST_LIFE<-FD002_TEST_OBS+FD002_TEST_RUL
FD003_TEST_LIFE<-FD003_TEST_OBS+FD003_TEST_RUL
FD004_TEST_LIFE<-FD004_TEST_OBS+FD004_TEST_RUL

#FD001 argvals
FD001_train<-CMAPSS[["train"]][["N"]][1:100]
FD001_test<-CMAPSS[["test"]][["N"]][1:100]
length(FD001_train)
length(FD001_test)
max(FD001_train)
FD001argvalssmooth<-matrix(NA,max(FD001_train),length(FD001_train))
for (i in 1:length(FD001_train)) {
  FD001argvalssmooth[1:FD001_train[i],i]<-seq(0, 1, length.out=FD001_train[i])
}

#FD002 argvals
FD002_train<-CMAPSS[["train"]][["N"]][101:360]
FD002_test<-CMAPSS[["test"]][["N"]][101:359]
length(FD002_train)
length(FD002_test)
max(FD002_train)
FD002argvalssmooth<-matrix(NA,max(FD002_train),length(FD002_train))
for (i in 1:length(FD002_train)) {
  FD002argvalssmooth[1:FD002_train[i],i]<-seq(0, 1, length.out=FD002_train[i])
}

#FD003 argvals
FD003_train<-CMAPSS[["train"]][["N"]][361:460]
FD003_test<-CMAPSS[["test"]][["N"]][360:459]
length(FD003_train)
length(FD003_test)
max(FD003_train)
FD003argvalssmooth<-matrix(NA,max(FD003_train),length(FD003_train))
for (i in 1:length(FD003_train)) {
  FD003argvalssmooth[1:FD003_train[i],i]<-seq(0, 1, length.out=FD003_train[i])
}

#FD004 argvals
FD004_train<-CMAPSS[["train"]][["N"]][461:709]
FD004_test<-CMAPSS[["train"]][["N"]][460:707]
length(FD004_train)
length(FD004_test)
max(FD004_train)
FD004argvalssmooth<-matrix(NA,max(FD004_train),length(FD004_train))
for (i in 1:length(FD004_train)) {
  FD004argvalssmooth[1:FD004_train[i],i]<-seq(0, 1, length.out=FD004_train[i])
}


##All sensors train & test
allsensors_train<-CMAPSS[["train"]][["x"]]
allsensors_test<-CMAPSS[["test"]][["x"]]

##All train sensors
all_T24_train<-allsensors_train[,"lpc.temp"]
all_T30_train<-allsensors_train[,"hpc.temp"]
all_T50_train<-allsensors_train[,"out.temp"]
all_P30_train<-allsensors_train[,"hpc.pres"]
all_ps30_train<-allsensors_train[,"stat.pres"]
all_phi_train<-allsensors_train[,"phi"]
all_BPR_train<-allsensors_train[,"bypass.ratio"]
all_W31_train<-allsensors_train[,"hpt.bleed"]
all_W32_train<-allsensors_train[,"lpt.bleed"]

##All test sensors
all_T24_test<-allsensors_test[,"lpc.temp"]
all_T30_test<-allsensors_test[,"hpc.temp"]
all_T50_test<-allsensors_test[,"out.temp"]
all_P30_test<-allsensors_test[,"hpc.pres"]
all_ps30_test<-allsensors_test[,"stat.pres"]
all_phi_test<-allsensors_test[,"phi"]
all_BPR_test<-allsensors_test[,"bypass.ratio"]
all_W31_test<-allsensors_test[,"hpt.bleed"]
all_W32_test<-allsensors_test[,"lpt.bleed"]


#Train sensors starting points
FD001START_trn<-1
FD002START_trn<-sum(FD001_train)+1
FD003START_trn<-sum(FD001_train,FD002_train)+1
FD004START_trn<-sum(FD001_train,FD002_train,FD003_train)+1

#Test sensors starting points
FD001START_tst<-1
FD002START_tst<-sum(FD001_test)+1
FD003START_tst<-sum(FD001_test,FD002_test)+1
FD004START_tst<-sum(FD001_test,FD002_test,FD003_test)+1


##FD001 T24 matrix generation
##FD001 T24 matrix generation
##FD001 T24 matrix generation
FD001_T24_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_T24_all)

##FD001 T24 all
for (i in 1:length(FD001_train)) {
  FD001_T24_all[i,1:FD001_train[i]]<-all_T24_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_T24_all[i+length(FD001_train),1:FD001_test[i]]<-all_T24_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 T24 train and test separated
FD001_T24_train<-FD001_T24_all[1:length(FD001_train),]
FD001_T24_test<-FD001_T24_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_T24_train)
dim(FD001_T24_test)

###################################################

##FD002 T24 matrix generation
##FD002 T24 matrix generation
##FD002 T24 matrix generation
FD002_T24_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_T24_all)

##FD002 T24 all
for (i in 1:length(FD002_train)) {
  FD002_T24_all[i,1:FD002_train[i]]<-all_T24_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD002_T24_all[i+length(FD002_train),1:FD002_test[i]]<-all_T24_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}

##FD002 T24 train and test separated
FD002_T24_train<-FD002_T24_all[1:length(FD002_train),]
FD002_T24_test<-FD002_T24_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_T24_test)
dim(FD002_T24_train)

###################################################

##Generate FD003 T24 matrix
FD003_T24_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_T24_all)

##FD003 T24 all
for (i in 1:length(FD003_train)) {
  FD003_T24_all[i,1:FD003_train[i]]<-all_T24_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD003_T24_all[i+length(FD003_train),1:FD003_test[i]]<-all_T24_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}

##FD003 T24 train and test separated
FD003_T24_train<-FD003_T24_all[1:length(FD003_train),]
FD003_T24_test<-FD003_T24_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_T24_test)
dim(FD003_T24_train)

###################################################

## Generate FD004 T24 matrix
## Generate FD004 T24 matrix
## Generate FD004 T24 matrix
FD004_T24_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_T24_all)

##FD004 T24 all
for (i in 1:length(FD004_train)) {
  FD004_T24_all[i,1:FD004_train[i]]<-all_T24_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD004_T24_all[i+length(FD004_train),1:FD004_test[i]]<-all_T24_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 T24 train and test separated
FD004_T24_train<-FD004_T24_all[1:length(FD004_train),]
FD004_T24_test<-FD004_T24_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_T24_train)
dim(FD004_T24_test)

###################################################

##FD001 T30 matrix generation
##FD001 T30 matrix generation
##FD001 T30 matrix generation
FD001_T30_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_T30_all)

##FD001 T30 all
for (i in 1:length(FD001_train)) {
  FD001_T30_all[i,1:FD001_train[i]]<-all_T30_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_T30_all[i+length(FD001_train),1:FD001_test[i]]<-all_T30_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 T30 train and test separated
FD001_T30_train<-FD001_T30_all[1:length(FD001_train),]
FD001_T30_test<-FD001_T30_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_T30_train)
dim(FD001_T30_test)

###################################################

##FD002 T30 matrix generation
##FD002 T30 matrix generation
##FD002 T30 matrix generation
FD002_T30_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_T30_all)

##FD002 T30 all
for (i in 1:length(FD002_train)) {
  FD002_T30_all[i,1:FD002_train[i]]<-all_T30_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD002_T30_all[i+length(FD002_train),1:FD002_test[i]]<-all_T30_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}
##FD002 T30 train and test separated
FD002_T30_train<-FD002_T30_all[1:length(FD002_train),]
FD002_T30_test<-FD002_T30_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_T30_train)
dim(FD002_T30_test)

###################################################

##FD003 T30 matrix generation
##FD003 T30 matrix generation
##FD003 T30 matrix generation
FD003_T30_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_T30_all)

##FD003 T30 all
for (i in 1:length(FD003_train)) {
  FD003_T30_all[i,1:FD003_train[i]]<-all_T30_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD003_T30_all[i+length(FD003_train),1:FD003_test[i]]<-all_T30_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}
##FD003 T30 train and test separated
FD003_T30_train<-FD003_T30_all[1:length(FD003_train),]
FD003_T30_test<-FD003_T30_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_T30_train)
dim(FD003_T30_test)

###################################################

##FD004 T30 matrix generation
##FD004 T30 matrix generation
##FD004 T30 matrix generation
FD004_T30_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_T30_all)

##FD004 T30 all
for (i in 1:length(FD004_train)) {
  FD004_T30_all[i,1:FD004_train[i]]<-all_T30_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD004_T30_all[i+length(FD004_train),1:FD004_test[i]]<-all_T30_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 T30 train and test separated
FD004_T30_train<-FD004_T30_all[1:length(FD004_train),]
FD004_T30_test<-FD004_T30_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_T30_train)
dim(FD004_T30_test)

###################################################

##FD001 T50 matrix generation
##FD001 T50 matrix generation
##FD001 T50 matrix generation
FD001_T50_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_T50_all)

##FD001 T50 all
for (i in 1:length(FD001_train)) {
  FD001_T50_all[i,1:FD001_train[i]]<-all_T50_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_T50_all[i+length(FD001_train),1:FD001_test[i]]<-all_T50_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 T50 train and test separated
FD001_T50_train<-FD001_T50_all[1:length(FD001_train),]
FD001_T50_test<-FD001_T50_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_T50_train)
dim(FD001_T50_test)

###################################################

##FD002 T50 matrix generation
##FD002 T50 matrix generation
##FD002 T50 matrix generation
FD002_T50_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_T50_all)

##FD002 T50 all
for (i in 1:length(FD002_train)) {
  FD002_T50_all[i,1:FD002_train[i]]<-all_T50_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD002_T50_all[i+length(FD002_train),1:FD002_test[i]]<-all_T50_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}
##FD002 T50 train and test separated
FD002_T50_train<-FD002_T50_all[1:length(FD002_train),]
FD002_T50_test<-FD002_T50_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_T50_train)
dim(FD002_T50_test)

###################################################

##FD003 T50 matrix generation
##FD003 T50 matrix generation
##FD003 T50 matrix generation
FD003_T50_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_T50_all)

##FD003 T50 all
for (i in 1:length(FD003_train)) {
  FD003_T50_all[i,1:FD003_train[i]]<-all_T50_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD003_T50_all[i+length(FD003_train),1:FD003_test[i]]<-all_T50_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}
##FD003 T50 train and test separated
FD003_T50_train<-FD003_T50_all[1:length(FD003_train),]
FD003_T50_test<-FD003_T50_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_T50_train)
dim(FD003_T50_test)

###################################################

##FD004 T50 matrix generation
##FD004 T50 matrix generation
##FD004 T50 matrix generation
FD004_T50_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_T50_all)

##FD004 T50 all
for (i in 1:length(FD004_train)) {
  FD004_T50_all[i,1:FD004_train[i]]<-all_T50_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD004_test)) {
  FD004_T50_all[(i+length(FD004_train)),(1:FD004_test[i])]<-all_T50_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 T50 train and test separated
FD004_T50_train<-FD004_T50_all[1:length(FD004_train),]
FD004_T50_test<-FD004_T50_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_T50_train)
dim(FD004_T50_test)

###################################################

##FD001 P30 matrix generation
##FD001 P30 matrix generation
##FD001 P30 matrix generation
FD001_P30_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_P30_all)

##FD001 P30 all
for (i in 1:length(FD001_train)) {
  FD001_P30_all[i,1:FD001_train[i]]<-all_P30_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_P30_all[i+length(FD001_train),1:FD001_test[i]]<-all_P30_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 P30 train and test separated
FD001_P30_train<-FD001_P30_all[1:length(FD001_train),]
FD001_P30_test<-FD001_P30_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_P30_train)
dim(FD001_P30_test)

###################################################

##FD002 P30 matrix generation
##FD002 P30 matrix generation
##FD002 P30 matrix generation
FD002_P30_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_P30_all)

##FD002 P30 all
for (i in 1:length(FD002_train)) {
  FD002_P30_all[i,1:FD002_train[i]]<-all_P30_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD002_test)) {
  FD002_P30_all[i+length(FD002_train),1:FD002_test[i]]<-all_P30_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}
##FD002 P30 train and test separated
FD002_P30_train<-FD002_P30_all[1:length(FD002_train),]
FD002_P30_test<-FD002_P30_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_P30_train)
dim(FD002_P30_test)

###################################################

##FD003 P30 matrix generation
##FD003 P30 matrix generation
##FD003 P30 matrix generation
FD003_P30_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_P30_all)

##FD003 P30 all
for (i in 1:length(FD003_train)) {
  FD003_P30_all[i,1:FD003_train[i]]<-all_P30_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD003_test)) {
  FD003_P30_all[i+length(FD003_train),1:FD003_test[i]]<-all_P30_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}
##FD003 P30 train and test separated
FD003_P30_train<-FD003_P30_all[1:length(FD003_train),]
FD003_P30_test<-FD003_P30_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_P30_train)
dim(FD003_P30_test)

###################################################

##FD004 P30 matrix generation
##FD004 P30 matrix generation
##FD004 P30 matrix generation
FD004_P30_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_P30_all)

##FD004 P30 all
for (i in 1:length(FD004_train)) {
  FD004_P30_all[i,1:FD004_train[i]]<-all_P30_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD004_test)) {
  FD004_P30_all[i+length(FD004_train),1:FD004_test[i]]<-all_P30_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 P30 train and test separated
FD004_P30_train<-FD004_P30_all[1:length(FD004_train),]
FD004_P30_test<-FD004_P30_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_P30_train)
dim(FD004_P30_test)

###################################################

##FD001 ps30 matrix generation
##FD001 ps30 matrix generation
##FD001 ps30 matrix generation
FD001_ps30_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_ps30_all)

##FD001 ps30 all
for (i in 1:length(FD001_train)) {
  FD001_ps30_all[i,1:FD001_train[i]]<-all_ps30_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_ps30_all[i+length(FD001_train),1:FD001_test[i]]<-all_ps30_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 ps30 train and test separated
FD001_ps30_train<-FD001_ps30_all[1:length(FD001_train),]
FD001_ps30_test<-FD001_ps30_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_ps30_train)
dim(FD001_ps30_test)

###################################################

##FD002 ps30 matrix generation
FD002_ps30_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_ps30_all)

##FD002 ps30 all
for (i in 1:length(FD002_train)) {
  FD002_ps30_all[i,1:FD002_train[i]]<-all_ps30_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD002_test)) {
  FD002_ps30_all[i+length(FD002_train),1:FD002_test[i]]<-all_ps30_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}
##FD002 ps30 train and test separated
FD002_ps30_train<-FD002_ps30_all[1:length(FD002_train),]
FD002_ps30_test<-FD002_ps30_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_ps30_train)
dim(FD002_ps30_test)

###################################################

##FD003 ps30 matrix generation
##FD003 ps30 matrix generation
##FD003 ps30 matrix generation
FD003_ps30_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_ps30_all)

##FD003 ps30 all
for (i in 1:length(FD003_train)) {
  FD003_ps30_all[i,1:FD003_train[i]]<-all_ps30_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD003_test)) {
  FD003_ps30_all[i+length(FD003_train),1:FD003_test[i]]<-all_ps30_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}
##FD003 ps30 train and test separated
FD003_ps30_train<-FD003_ps30_all[1:length(FD003_train),]
FD003_ps30_test<-FD003_ps30_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_ps30_train)
dim(FD003_ps30_test)

###################################################

##FD004 ps30 matrix generation
FD004_ps30_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_ps30_all)

##FD004 ps30 all
for (i in 1:length(FD004_train)) {
  FD004_ps30_all[i,1:FD004_train[i]]<-all_ps30_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD004_test)) {
  FD004_ps30_all[i+length(FD004_train),1:FD004_test[i]]<-all_ps30_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 ps30 train and test separated
FD004_ps30_train<-FD004_ps30_all[1:length(FD004_train),]
FD004_ps30_test<-FD004_ps30_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_ps30_train)
dim(FD004_ps30_test)

###################################################

##FD001 phi matrix generation
##FD001 phi matrix generation
##FD001 phi matrix generation
FD001_phi_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_phi_all)

##FD001 phi all
for (i in 1:length(FD001_train)) {
  FD001_phi_all[i,1:FD001_train[i]]<-all_phi_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_phi_all[i+length(FD001_train),1:FD001_test[i]]<-all_phi_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 phi train and test separated
FD001_phi_train<-FD001_phi_all[1:length(FD001_train),]
FD001_phi_test<-FD001_phi_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_phi_train)
dim(FD001_phi_test)

###################################################

##FD002 phi matrix generation
##FD002 phi matrix generation
##FD002 phi matrix generation
FD002_phi_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_phi_all)

##FD002 phi all
for (i in 1:length(FD002_train)) {
  FD002_phi_all[i,1:FD002_train[i]]<-all_phi_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD002_test)) {
  FD002_phi_all[i+length(FD002_train),1:FD002_test[i]]<-all_phi_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}
##FD002 phi train and test separated
FD002_phi_train<-FD002_phi_all[1:length(FD002_train),]
FD002_phi_test<-FD002_phi_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_phi_train)
dim(FD002_phi_test)

###################################################

##FD003 phi matrix generation
##FD003 phi matrix generation
##FD003 phi matrix generation
FD003_phi_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_phi_all)

##FD003 phi all
for (i in 1:length(FD003_train)) {
  FD003_phi_all[i,1:FD003_train[i]]<-all_phi_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD003_test)) {
  FD003_phi_all[i+length(FD003_train),1:FD003_test[i]]<-all_phi_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}
##FD003 phi train and test separated
FD003_phi_train<-FD003_phi_all[1:length(FD003_train),]
FD003_phi_test<-FD003_phi_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_phi_train)
dim(FD003_phi_test)

###################################################

##FD004 phi matrix generation
##FD004 phi matrix generation
##FD004 phi matrix generation
FD004_phi_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_phi_all)

##FD004 phi all
for (i in 1:length(FD004_train)) {
  FD004_phi_all[i,1:FD004_train[i]]<-all_phi_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD004_test)) {
  FD004_phi_all[i+length(FD004_train),1:FD004_test[i]]<-all_phi_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 phi train and test separated
FD004_phi_train<-FD004_phi_all[1:length(FD004_train),]
FD004_phi_test<-FD004_phi_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_phi_train)
dim(FD004_phi_test)

###################################################

##FD001 BPR matrix generation
##FD001 BPR matrix generation
##FD001 BPR matrix generation
FD001_BPR_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_BPR_all)

##FD001 BPR all
for (i in 1:length(FD001_train)) {
  FD001_BPR_all[i,1:FD001_train[i]]<-all_BPR_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_BPR_all[i+length(FD001_train),1:FD001_test[i]]<-all_BPR_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 BPR train and test separated
FD001_BPR_train<-FD001_BPR_all[1:length(FD001_train),]
FD001_BPR_test<-FD001_BPR_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_BPR_train)
dim(FD001_BPR_test)

###################################################

##FD002 BPR matrix generation
##FD002 BPR matrix generation
##FD002 BPR matrix generation
FD002_BPR_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_BPR_all)

##FD002 BPR all
for (i in 1:length(FD002_train)) {
  FD002_BPR_all[i,1:FD002_train[i]]<-all_BPR_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD002_test)) {
  FD002_BPR_all[i+length(FD002_train),1:FD002_test[i]]<-all_BPR_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}
##FD002 BPR train and test separated
FD002_BPR_train<-FD002_BPR_all[1:length(FD002_train),]
FD002_BPR_test<-FD002_BPR_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_BPR_train)
dim(FD002_BPR_test)

###################################################

##FD003 BPR matrix generation
##FD003 BPR matrix generation
##FD003 BPR matrix generation
FD003_BPR_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_BPR_all)

##FD003 BPR all
for (i in 1:length(FD003_train)) {
  FD003_BPR_all[i,1:FD003_train[i]]<-all_BPR_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD003_test)) {
  FD003_BPR_all[i+length(FD003_train),1:FD003_test[i]]<-all_BPR_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}
##FD003 BPR train and test separated
FD003_BPR_train<-FD003_BPR_all[1:length(FD003_train),]
FD003_BPR_test<-FD003_BPR_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_BPR_train)
dim(FD003_BPR_test)

###################################################

##FD004 BPR matrix generation
##FD004 BPR matrix generation
##FD004 BPR matrix generation
FD004_BPR_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_BPR_all)

##FD004 BPR all
for (i in 1:length(FD004_train)) {
  FD004_BPR_all[i,1:FD004_train[i]]<-all_BPR_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD004_test)) {
  FD004_BPR_all[i+length(FD004_train),1:FD004_test[i]]<-all_BPR_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 BPR train and test separated
FD004_BPR_train<-FD004_BPR_all[1:length(FD004_train),]
FD004_BPR_test<-FD004_BPR_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_BPR_train)
dim(FD004_BPR_test)

###################################################

##FD001 W31 matrix generation
##FD001 W31 matrix generation
##FD001 W31 matrix generation
FD001_W31_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_W31_all)

##FD001 W31 all
for (i in 1:length(FD001_train)) {
  FD001_W31_all[i,1:FD001_train[i]]<-all_W31_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_W31_all[i+length(FD001_train),1:FD001_test[i]]<-all_W31_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 W31 train and test separated
FD001_W31_train<-FD001_W31_all[1:length(FD001_train),]
FD001_W31_test<-FD001_W31_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_W31_train)
dim(FD001_W31_test)

###################################################

##FD002 W31 matrix generation
##FD002 W31 matrix generation
##FD002 W31 matrix generation
FD002_W31_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_W31_all)

##FD002 W31 all
for (i in 1:length(FD002_train)) {
  FD002_W31_all[i,1:FD002_train[i]]<-all_W31_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD002_test)) {
  FD002_W31_all[i+length(FD002_train),1:FD002_test[i]]<-all_W31_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}
##FD002 W31 train and test separated
FD002_W31_train<-FD002_W31_all[1:length(FD002_train),]
FD002_W31_test<-FD002_W31_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_W31_train)
dim(FD002_W31_test)

###################################################

##FD003 W31 matrix generation
##FD003 W31 matrix generation
##FD003 W31 matrix generation
FD003_W31_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_W31_all)

##FD003 W31 all
for (i in 1:length(FD003_train)) {
  FD003_W31_all[i,1:FD003_train[i]]<-all_W31_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD003_test)) {
  FD003_W31_all[i+length(FD003_train),1:FD003_test[i]]<-all_W31_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}
##FD003 W31 train and test separated
FD003_W31_train<-FD003_W31_all[1:length(FD003_train),]
FD003_W31_test<-FD003_W31_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_W31_train)
dim(FD003_W31_test)

###################################################

##FD004 W31 matrix generation
##FD004 W31 matrix generation
##FD004 W31 matrix generation
FD004_W31_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_W31_all)

##FD004 W31 all
for (i in 1:length(FD004_train)) {
  FD004_W31_all[i,1:FD004_train[i]]<-all_W31_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD004_test)) {
  FD004_W31_all[i+length(FD004_train),1:FD004_test[i]]<-all_W31_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 W31 train and test separated
FD004_W31_train<-FD004_W31_all[1:length(FD004_train),]
FD004_W31_test<-FD004_W31_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_W31_train)
dim(FD004_W31_test)

###################################################

##FD001 W32 matrix generation
##FD001 W32 matrix generation
##FD001 W32 matrix generation
FD001_W32_all<-matrix(NA,length(FD001_train)+length(FD001_test),max(FD001_train,FD001_test))
dim(FD001_W32_all)

##FD001 W32 all
for (i in 1:length(FD001_train)) {
  FD001_W32_all[i,1:FD001_train[i]]<-all_W32_train[(FD001START_trn+(sum(FD001_train[0:(i-1)]))):(FD001START_trn-1+sum(FD001_train[0:i]))]
}
for (i in 1:length(FD001_test)) {
  FD001_W32_all[i+length(FD001_train),1:FD001_test[i]]<-all_W32_test[(FD001START_tst+(sum(FD001_test[0:(i-1)]))):(FD001START_tst-1+sum(FD001_test[0:i]))]
}
##FD001 W32 train and test separated
FD001_W32_train<-FD001_W32_all[1:length(FD001_train),]
FD001_W32_test<-FD001_W32_all[(length(FD001_train)+1):(length(FD001_train)+length(FD001_test)),]
dim(FD001_W32_train)
dim(FD001_W32_test)

###################################################

##FD002 W32 matrix generation
##FD002 W32 matrix generation
##FD002 W32 matrix generation
FD002_W32_all<-matrix(NA,length(FD002_train)+length(FD002_test),max(FD002_train,FD002_test))
dim(FD002_W32_all)

##FD002 W32 all
for (i in 1:length(FD002_train)) {
  FD002_W32_all[i,1:FD002_train[i]]<-all_W32_train[(FD002START_trn+(sum(FD002_train[0:(i-1)]))):(FD002START_trn-1+sum(FD002_train[0:i]))]
}
for (i in 1:length(FD002_test)) {
  FD002_W32_all[i+length(FD002_train),1:FD002_test[i]]<-all_W32_test[(FD002START_tst+(sum(FD002_test[0:(i-1)]))):(FD002START_tst-1+sum(FD002_test[0:i]))]
}
##FD002 W32 train and test separated
FD002_W32_train<-FD002_W32_all[1:length(FD002_train),]
FD002_W32_test<-FD002_W32_all[(length(FD002_train)+1):(length(FD002_train)+length(FD002_test)),]
dim(FD002_W32_train)
dim(FD002_W32_test)

###################################################

##FD003 W32 matrix generation
##FD003 W32 matrix generation
##FD003 W32 matrix generation
FD003_W32_all<-matrix(NA,length(FD003_train)+length(FD003_test),max(FD003_train,FD003_test))
dim(FD003_W32_all)

##FD003 W32 all
for (i in 1:length(FD003_train)) {
  FD003_W32_all[i,1:FD003_train[i]]<-all_W32_train[(FD003START_trn+(sum(FD003_train[0:(i-1)]))):(FD003START_trn-1+sum(FD003_train[0:i]))]
}
for (i in 1:length(FD003_test)) {
  FD003_W32_all[i+length(FD003_train),1:FD003_test[i]]<-all_W32_test[(FD003START_tst+(sum(FD003_test[0:(i-1)]))):(FD003START_tst-1+sum(FD003_test[0:i]))]
}
##FD003 W32 train and test separated
FD003_W32_train<-FD003_W32_all[1:length(FD003_train),]
FD003_W32_test<-FD003_W32_all[(length(FD003_train)+1):(length(FD003_train)+length(FD003_test)),]
dim(FD003_W32_train)
dim(FD003_W32_test)

###################################################

##FD004 W32 matrix generation
FD004_W32_all<-matrix(NA,length(FD004_train)+length(FD004_test),max(FD004_train,FD004_test))
dim(FD004_W32_all)

##FD004 W32 all
for (i in 1:length(FD004_train)) {
  FD004_W32_all[i,1:FD004_train[i]]<-all_W32_train[(FD004START_trn+(sum(FD004_train[0:(i-1)]))):(FD004START_trn-1+sum(FD004_train[0:i]))]
}
for (i in 1:length(FD004_test)) {
  FD004_W32_all[i+length(FD004_train),1:FD004_test[i]]<-all_W32_test[(FD004START_tst+(sum(FD004_test[0:(i-1)]))):(FD004START_tst-1+sum(FD004_test[0:i]))]
}
##FD004 W32 train and test separated
FD004_W32_train<-FD004_W32_all[1:length(FD004_train),]
FD004_W32_test<-FD004_W32_all[(length(FD004_train)+1):(length(FD004_train)+length(FD004_test)),]
dim(FD004_W32_train)
dim(FD004_W32_test)

###################################################

# -----------------------------------------------------------------------------
# 4. CONTROL THE DATA
# -----------------------------------------------------------------------------

FD001argvalssmooth
FD002argvalssmooth
FD003argvalssmooth
FD004argvalssmooth

FD001_T24_train
FD001_T30_train
FD001_T50_train
FD001_P30_train
FD001_ps30_train
FD001_phi_train
FD001_BPR_train
FD001_W31_train
FD001_W32_train

FD001_T24_test
FD001_T30_test
FD001_T50_test
FD001_P30_test
FD001_ps30_test
FD001_phi_test
FD001_BPR_test
FD001_W31_test
FD001_W32_test

FD002_T24_train
FD002_T30_train
FD002_T50_train
FD002_P30_train
FD002_ps30_train
FD002_phi_train
FD002_BPR_train
FD002_W31_train
FD002_W32_train

FD002_T24_test
FD002_T30_test
FD002_T50_test
FD002_P30_test
FD002_ps30_test
FD002_phi_test
FD002_BPR_test
FD002_W31_test
FD002_W32_test

FD003_T24_train
FD003_T30_train
FD003_T50_train
FD003_P30_train
FD003_ps30_train
FD003_phi_train
FD003_BPR_train
FD003_W31_train
FD003_W32_train

FD003_T24_test
FD003_T30_test
FD003_T50_test
FD003_P30_test
FD003_ps30_test
FD003_phi_test
FD003_BPR_test
FD003_W31_test
FD003_W32_test

FD004_T24_train
FD004_T30_train
FD004_T50_train
FD004_P30_train
FD004_ps30_train
FD004_phi_train
FD004_BPR_train
FD004_W31_train
FD004_W32_train

FD004_T24_test
FD004_T30_test
FD004_T50_test
FD004_P30_test
FD004_ps30_test
FD004_phi_test
FD004_BPR_test
FD004_W31_test
FD004_W32_test

###ALL CMAPPS DATA ready