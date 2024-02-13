### This script make obtains deconvolution data for subtype VII to plot
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
load("modscorecomzAUC.RData")
load("pd_for_mut.RData") ## this is metadata with mutation information
load("fpmall.RData")
fpm.all[1:4,1:4]
########### proption plot############
##selecting MYC affected tumors
metax <- meta[row.names(plot.data),]
keep1 <- plot.data$Subtype_M=="G34_VII"

subs <- read.csv("subtypeG34.csv", row.names = 1)
d <- sort(unique(pd$Subtype_M))
sc <- subs[d,1]

###Composition of Axes states
load("CR_thetaAnn.RData")
CR.theta <- theta_t
rm(theta_t, CR)
load("ICGC_thetaAnn.RData")
ICGC.theta <- theta_t
rm(theta_t, ICGC)
load("MDT_thetaAnn.RData")
MDT.theta <- theta_t
rm(theta_t, MDT)

theta <- rbind(CR.theta,ICGC.theta, MDT.theta)
theta <- theta[row.names(pd),]

####
lv2 <- c("PR_def","PR_early","PR_early_MYC",
         "MYC","MYC_CC","Prec_CC","Precursors","UBC_early","UBC_late")
##removing normal cells
theta <- theta[,lv2]
##converting to 1
sc <- rowSums(theta)
theta <- theta/sc

thetay <- data.frame(theta[row.names(pdy),])
thetay$Subtype <- pdy$Subtype_M
thetay$ID <- row.names(thetay)
thetay2 <- reshape2::melt(thetay)
head(thetay2)

##order sabple by high PRt to High UBCt score
ms2sel <- ms2[row.names(pdy),]
score <- ms2sel[,"CRX"]- ms2sel[,"LHX1"]
ord <- order(score, decreasing = TRUE)

lv <- row.names(pdy)[ord]


ann <- read.delim("Files/anncolor.txt",row.names = 1)

anncol <- ann[lv2,1]
thetay2$ID <- factor(thetay2$ID, levels = lv)
thetay2$variable <- factor(thetay2$variable, levels = lv2)
thetay2$Subtype <- factor(thetay2$Subtype, levels = lv3)

p <- ggplot(thetay2, aes(x=ID, y=value,fill=variable))+
  scale_fill_manual(values=anncol)+
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',linewidth = .8),
        axis.ticks = element_line(colour = 'black',linewidth=.8),
        axis.text.y= element_text(face="bold",colour = "black",size=10),
        axis.text.x = element_text(face="bold",colour = "black",size=8, angle = 90,vjust=0.5, hjust=1),
        axis.title=element_blank(),
        strip.background = element_blank(),
        legend.position = "None")+
  scale_y_continuous(breaks =c(0, 0.5, 1))

tiff("SVII_Ann_prop_bulk.tiff", width = 18, height = 3, units = "in", res = 300)
print(p)
dev.off()

