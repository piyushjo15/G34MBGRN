library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(batchelor)
})

##load data -----
##cosine normed tumor data
load("cosnormG34.RData")

load("plotdataG34MB.RData")
plot.data.MB <- plot.data
rm(plot.data)

##load reference cerebllum data
load("referencedata.RData") ##this is a list of Cb RNA libraries, CBsce.list
rm(dec.list)

sams <- unique(plot.data.MB$Batch) ##tumor samples
#v2
common.genes <- readLines("Files/topHVGHUMCBnoccSX.txt")
com.hvg <- readLines("Files/topHVGG34new10xnoRPMTSX.txt")
com.hvg <- intersect(common.genes[1:5000],com.hvg)
# 
sams2 <- names(CBsce.list)
CB_cnnc <- c()
for(y in sams2){
  del <- logcounts(CBsce.list[[y]])
  del <- cosineNorm(del[com.hvg,])
  CB_cnnc <- cbind(CB_cnnc,del)
  rm(del)
}
rm(y)

Cluster <- plot.data$dev_state ##reference cluster annotation
rm(CBsce.list, plot.data, sams, all.sce)
## MB 
mdt <- t(CB_cnnc)
mdt_test <- t(com.sce_cnnc[com.hvg,])
rm(com.sce_cnnc)
##run SVM-----

repl_python()
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
import random
X=r.mdt
y=r.Cluster
X_test=r.mdt_test
CLFA = LinearSVC(max_iter=100000, class_weight='balanced', random_state=0)
SVM_c = CalibratedClassifierCV(CLFA, method='isotonic', n_jobs=-1)

#Class weight SVM
random.seed(145)
SVM_c.fit(X,y)

test_l1 = SVM_c.predict(X_test)
test_prob = SVM_c.predict_proba(X_test)

exit
a <-  py$test_l1
b <- py$test_prob
plot.data <- plot.data.MB
plot.data$SVM_prilab <- plot.data$SVM_finlab <- a
prob <- b
prob <- apply(prob, 1, max)
plot.data$SVM_prob <- prob
keep <- prob >=0.5
plot.data[!keep,]$SVM_finlab <- "ND"
save(plot.data,a,b, file = "SVMpredG34MB.RData")

q()
###Second round of SVM prediction ######
###running second round of prediction by using tumor cells with high confidence
## calls as reference to label tumor cells with lower confidence call

## using entropy and probability to identify cells labelled with high confidence
load("SVMpredG34MB.RData")
Entropy <- function(vec){
  vec <- vec[vec>0]
  return(-sum(vec*log2(vec)))
}

## entropy
plot.data$SVM_EN <- apply(b,1,Entropy)

keep1 <- plot.data$SVM_prob>0.7
table(keep1)
#plot.data$Col <- keep1

con <- quantile(plot.data$SVM_EN,0.4)
keep2 <- plot.data$SVM_EN<con
table(keep2)
#plot.data$Col2 <- keep2

keep3 <- keep1 & keep2

table(keep3)

pd_train <- plot.data[keep3,]

pd_test <- plot.data[!keep3,]

tb <- data.frame(table(pd_train$SVM_prilab))
tbx <- tb$Freq<10
tbx <- as.character(tb[tbx,1])

keep5 <- pd_train$SVM_prilab %in% tbx
table(keep5)
pd_trainx <- pd_train[keep5,]
pd_train <- pd_train[!keep5,]

com.hvg <- lreadLines("Files/topHVGG34new10xnoRPMTSX.txt")
load("cosnormG34.RData")


mdt_train <- t(com.sce_cnnc[com.hvg,row.names(pd_train)])
Cluster <- pd_train$SVM_prilab

mdt_test <- t(com.sce_cnnc[com.hvg,row.names(pd_test)])

rm(com.sce_cnnc, keep1,keep2, keep3, keep5, tbx)

##run SVM-----

repl_python()
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
import random

X=r.mdt_train
X_test = r.mdt_test
y=r.Cluster


CLFA = LinearSVC(max_iter=100000, class_weight='balanced', random_state=0)
SVM_c = CalibratedClassifierCV(CLFA, method='isotonic', n_jobs=-1)

#Class weight SVM
random.seed(145)
SVM_c.fit(X,y)

test_l1 = SVM_c.predict(X_test)
test_prob = SVM_c.predict_proba(X_test)

exit
a <-  py$test_l1
b <- py$test_prob
pd_test$SVM_EN <- apply(b,1,Entropy)

pd_test$SVM_prilab <-  a
prob <- b
prob <- apply(prob, 1, max)
pd_test$SVM_prob <- prob

pd <- rbind(pd_test,pd_train, pd_trainx)
pd <- pd[row.names(plot.data),]
save(plot.data, pd_test, pd,  file = "SVMpredG34MB_KD.RData")
table(row.names(plot.data)==row.names(pd))
plot.data$SVM_prilabKD <- pd$SVM_prilab
plot.data$SVM_probKD <- pd$SVM_prob
#### identify cells that are non-neuronal nromal
plot.data$Nor <- "No"
keep <- plot.data$SVM_prilabKD %in%  c("astrocyte",'mural/endoth',"immune","meningeal",
                                       "oligo_progenitor","erythroid","oligodendrocyte",
                                       "glioblast")
table(keep)
plot.data[keep,"Nor"] <- "Yes"
save(plot.data, file = "plotdataG34_SVM.RData")

q()
