## Analyses of STS treatment success from Walter et al. (202X) in Pest Management Science
## Contact: Jonathan Walter, jaw3es@virginia.edu

rm(list=ls())

library(rgdal)
library(dplyr)
library(lubridate)
library(car)
library(ncf)
library(betareg)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

eval07<-read.csv("./data/TreatmentEvaluation_2007.csv", stringsAsFactors = F, na.strings="")
eval08<-read.csv("./data/TreatmentEvaluation_2008.csv", stringsAsFactors = F, na.strings="")
eval09<-read.csv("./data/TreatmentEvaluation_2009.csv", stringsAsFactors = F, na.strings="")
eval10<-read.csv("./data/TreatmentEvaluation_2010.csv", stringsAsFactors = F, na.strings="")
eval11<-read.csv("./data/TreatmentEvaluation_2011.csv", stringsAsFactors = F, na.strings="")
eval12<-read.csv("./data/TreatmentEvaluation_2012.csv", stringsAsFactors = F, na.strings="")
eval13<-read.csv("./data/TreatmentEvaluation_2013.csv", stringsAsFactors = F, na.strings="")
eval14<-read.csv("./data/TreatmentEvaluation_2014.csv", stringsAsFactors = F, na.strings="")
eval15<-read.csv("./data/TreatmentEvaluation_2015.csv", stringsAsFactors = F, na.strings="")
eval16<-read.csv("./data/TreatmentEvaluation_2016.csv", stringsAsFactors = F, na.strings="")
eval17<-read.csv("./data/TreatmentEvaluation_2017.csv", stringsAsFactors = F, na.strings="")
eval18<-read.csv("./data/TreatmentEvaluation_2018.csv", stringsAsFactors = F, na.strings="")

eval.all<-rbind(eval07,eval08,eval09,eval10,eval11,eval12,eval13,eval14,eval15,eval16,
                eval17,eval18)
eval.all$stblk<-paste(eval.all$State,eval.all$Blockname,sep="-")


#eval.va<-eval.all[eval.all$State=="VA",]


# eval.dups<-eval.all[duplicated(eval.all$stblk),]
# hist(eval.dups$T)

rm(eval07,eval08,eval09,eval10,eval11,eval12,eval13,eval14,eval15,eval16,eval17,eval18)

trtshp<-readOGR("./data/treatments_plus.shp", stringsAsFactors = F)
trtinf<-trtshp@data
trtinf<-trtinf[,colnames(trtinf) %in% c("OBJECTI","ACRES","ID","BLOCKNA","DOSAGE","PPA","STATE","YEAR",
                                        "TREATME","DATE_TR","DISCANC","rh100","pai","trpctch_v","trpctch_m",
                                        "hostba","hostsdi")]

# trthostba<-trtshp@data
# trthostba<-trthostba[,colnames(trthostba) %in% c("ID","BLOCKNA","STATE","hostba")]
trtinf$stblk<-paste(trtinf$STATE,trtinf$BLOCKNA,sep="-")
# trthostba<-aggregate(trthostba$hostba~trthostba$stblk, FUN="mean", na.rm=T)
# colnames(trthostba)<-c("stblk","hostba")

trt_ppa<-read.csv("./data/treatment_ppas_join.csv", stringsAsFactors = F)
trt_ppa<-trt_ppa[,colnames(trt_ppa) %in% c("BLOCKNAME","STATE","PRIORITY_I","RECOMMEND")]
trt_ppa$stblk<-paste(trt_ppa$STATE,trt_ppa$BLOCKNAME,sep="-")
#trt_ppa$PRIORITY_I[is.na(trt_ppa$PRIORITY_I)]<-0

eval.all<-left_join(eval.all,trtinf, by="stblk")
eval.all<-left_join(eval.all,trt_ppa, by="stblk")

eval.all$TreatYearFact<-as.factor(eval.all$TreatYear)

eval.all$Ta<-eval.all$T
eval.all$Ta["MD" %in% eval.all$Treatment]<-eval.all$T1["MD" %in% eval.all$Treatment]

eval.all$Ca<-eval.all$C
eval.all$Ca["MD" %in% eval.all$Treatment]<-eval.all$C1["MD" %in% eval.all$Treatment]

eval.all<-eval.all[!grepl("Study",eval.all$Blockname),]


## All years, initial treatments ------------------------------------------------------------------------

#eval.all$Ta[eval.all$Ta<0]<-0
eval.first<-eval.all[!duplicated(eval.all$stblk),] ## remove duplicates--look at those separately
#eval.first<-eval.first[complete.cases(eval.first),]
failures<-eval.first[eval.first$Ta < 0.33,]



table(failures$State, failures$TreatYear)
chisq.test(table(failures$State, failures$TreatYear))
sum(table(failures$State,failures$TreatYear))
1-107/1475 #92.7% of initial treatments are successful or partially successful.

cols=colorRampPalette(colors=c("white","darkred"))

layout(matrix(1:2,nrow=1,ncol=2),widths=c(0.825,0.175))

par(mar=c(4.1,5.1,1.1,0.1))
image(table(failures$State, failures$TreatYear), xaxt="n", yaxt="n", col=cols(11))
axis(1, seq(0,1,length.out=7), c("IL","IN","MN","NC","OH","VA","WI"))
axis(2, seq(0,1,length.out=12), as.character(2007:2018), las=2)
mtext("Year", 2, line=3.7)
mtext("State", 1, line=2.7)

par(mar=c(4.1,3.1,5.1,1.1))
image(t(matrix(1:11)) ,col=cols(11), xaxt="n", yaxt="n")
axis(2, seq(0,1,length.out=11), 0:10, las=2)

length(unique(eval.all$stblk[duplicated(eval.all$stblk)]))

length(unique(eval.all$stblk))

hist(eval.first$Ta)


table(eval.first$TreatYear)
table(eval.first$Treatment)
table(eval.first$State)



eval.first<-eval.first[eval.first$Treatment != "MIM" & eval.first$Treatment != "MIM-MD",]# & eval.first$Treatment != "DIM",]
eval.first<-eval.first[eval.first$State != "TN",]


#write.csv(eval.first$stblk,"studied_treatments_list.csv")

eval.first$Ta[eval.first$Ta==1]<-1-1e-5
eval.first$Ta[eval.first$Ta<=0]<-1e-5



testvars<-c("Ta","State","Treatment","DistanceKm","Acres","TreatYearFact","hostba","trpctch_m",
            "trpctch_v","hostsdi")
eval.first1<-eval.first[,colnames(eval.first) %in% testvars]
eval.first1<-eval.first1[complete.cases(eval.first1),]

sum(eval.first1$Ta==1-1e-5, na.rm=T)
sum(eval.first1$Ta==1e-5, na.rm=T)

tiff("fig1_histogram.tif",units="in",res=300,height=3.25,width=3.25)
par(mar=c(3.5,3.5,1.1,1.1), mgp=c(2,0.6,0))
hist(eval.first1$Ta, main="", xlab=expression(italic("T")),col="grey")
dev.off()

lm.Ta<-betareg(Ta ~ State + Treatment + DistanceKm + Acres + TreatYearFact +
          hostsdi + trpctch_m, data=eval.first1)
summary(lm.Ta)

lrt<-function(mod, submod){
  ll.mod<-mod$loglik
  df.mod<-length(unlist(mod$coefficients))
  ll.submod<-submod$loglik
  df.submod<-length(unlist(submod$coefficients))
  
  return(1-pchisq(2*(ll.mod-ll.submod), df=df.mod-df.submod))
}

lrt(lm.Ta, betareg(Ta ~ Treatment + DistanceKm + Acres + TreatYearFact + hostsdi + trpctch_m, data=eval.first1)) #test of State
lrt(lm.Ta, betareg(Ta ~ State + DistanceKm + Acres + TreatYearFact + hostsdi + trpctch_m, data=eval.first1)) #test of treatment
lrt(lm.Ta, betareg(Ta ~ State + Treatment + Acres + TreatYearFact + hostsdi + trpctch_m, data=eval.first1)) #test of distance
lrt(lm.Ta, betareg(Ta ~ State + Treatment + DistanceKm + TreatYearFact + hostsdi + trpctch_m, data=eval.first1)) #test of acres
lrt(lm.Ta, betareg(Ta ~ State + Treatment + DistanceKm + Acres + hostsdi + trpctch_m, data=eval.first1)) #test of treatment year
lrt(lm.Ta, betareg(Ta ~ State + Treatment + DistanceKm + Acres + TreatYearFact + trpctch_m, data=eval.first1)) #test of host basal area 
lrt(lm.Ta, betareg(Ta ~ State + Treatment + DistanceKm + Acres + TreatYearFact + hostsdi, data=eval.first1)) #test of trap catch


lm.Ta2<-betareg(Ta ~ Treatment + DistanceKm + Acres + TreatYearFact +
                 hostsdi + trpctch_m, data=eval.first1)
summary(lm.Ta2)

lm.Ta3<-betareg(Ta ~ State + DistanceKm + Acres + TreatYearFact +
                  hostsdi + trpctch_m, data=eval.first1)
summary(lm.Ta3)

tiff("fig2_model1boxplots.tif", units="in", res=300, width=3.25, height=7.5)
par(mfrow=c(3,1), mar=c(3.5,3.5,0.5,1.1), oma=c(0,0,0.6,0), mgp=c(2,0.6,0))
boxplot(Ta ~ State, data=eval.first1, ylab=expression(italic("T")), col="grey", pch=16, cex=0.7)
boxplot(Ta ~ TreatYearFact, data=eval.first1, ylab=expression(italic("T")), xlab="Treatment year", col="grey", pch=16, cex=0.7)
boxplot(Ta ~ Treatment, data=eval.first1, ylab=expression(italic("T")), 
        xlab="Treatment type",col="grey", pch=16, cex=0.7, cex.axis=0.9)
dev.off()


tiff("model1boxplots_horizontal.tif", units="in", res=300, width=8.5, height=3.25)
par(mfrow=c(1,3), mar=c(3.5,3.5,0.5,1.1), oma=c(0,0,0.6,0), mgp=c(2,0.6,0))
boxplot(Ta ~ State, data=eval.first1, ylab=expression(italic("T")), col="grey", pch=16, cex=0.7, cex.axis=0.85)
boxplot(Ta ~ TreatYearFact, data=eval.first1, ylab=expression(italic("T")), xlab="Treatment year", col="grey", pch=16, cex=0.7)
boxplot(Ta ~ Treatment, data=eval.first1, ylab=expression(italic("T")), 
        xlab="Treatment type",col="grey", pch=16, cex=0.7, cex.axis=0.9)
dev.off()

cor.test(eval.first1$DistanceKm, eval.first1$trpctch_m, method="spearman")
cor.test(eval.first1$DistanceKm, eval.first1$Acres, method="spearman")


plot(eval.first1$DistanceKm, eval.first1$Ta, ylab="Treatment success", xlab="Distance from invasion front")
lines(seq(min(eval.first1$DistanceKm), 
          max(eval.first1$DistanceKm), length.out=50), predict(betareg(Ta ~ DistanceKm, data=eval.first1),
                                      newdata=data.frame(DistanceKm=seq(min(eval.first1$DistanceKm), 
                                                                        max(eval.first1$DistanceKm), length.out=50))), col="red")



## Look at forest structure data ------------------------------------------------------
testvars2<-c("Ta","State","Treatment","DistanceKm","Acres","TreatYearFact","trpctch_m",
            "trpctch_v","pai","rh100")
eval.first2<-eval.first[,colnames(eval.first) %in% testvars2]
eval.first2<-eval.first2[complete.cases(eval.first2),]

cor(eval.first2$rh100, eval.first2$pai)

lm.Ta2a<-betareg(Ta ~ State + Treatment + DistanceKm + Acres + TreatYearFact +
                 rh100, data=eval.first2)

lm.Ta2b<-betareg(Ta ~ State + Treatment + DistanceKm + Acres + TreatYearFact +
                   pai, data=eval.first2)
AIC(lm.Ta2a)
AIC(lm.Ta2b)#slightly better

summary(lm.Ta2a)
summary(lm.Ta2b)

lrt<-function(mod, submod){
  ll.mod<-mod$loglik
  df.mod<-length(unlist(mod$coefficients))
  ll.submod<-submod$loglik
  df.submod<-length(unlist(submod$coefficients))
  
  return(1-pchisq(2*(ll.mod-ll.submod), df=df.mod-df.submod))
}

lrt(lm.Ta2b, betareg(Ta ~ Treatment + DistanceKm + Acres + TreatYearFact + pai, data=eval.first2)) #test on state
lrt(lm.Ta2b, betareg(Ta ~ State + DistanceKm + Acres + TreatYearFact + pai, data=eval.first2)) #test on treatment
lrt(lm.Ta2b, betareg(Ta ~ State + Treatment + Acres + TreatYearFact + pai, data=eval.first2)) #test on distance
lrt(lm.Ta2b, betareg(Ta ~ State + Treatment + DistanceKm + TreatYearFact + pai, data=eval.first2)) #test on acres
lrt(lm.Ta2b, betareg(Ta ~ State + Treatment + DistanceKm + TreatYearFact + Acres, data=eval.first2)) #test on pai


plot(eval.first2$pai, eval.first2$Ta, xlab="Plant area index", ylab="Treatment success")
lines(seq(min(eval.first2$pai), 
          max(eval.first2$pai), length.out=50), predict(betareg(Ta ~ pai, data=eval.first2),
                                                               newdata=data.frame(pai=seq(min(eval.first2$pai), 
                                                                                                 max(eval.first2$pai), length.out=50))), col="red")


## ------------------------------------------------------------------------------------
## Recent years, mating disruption ----------------------------------------------------
## here, test for dose effect, too




MDtrts.poly<-readOGR("./data/trt_MD_15_18.shp")
MDtrts<-MDtrts.poly@data

#recentMD<-right_join(eval.all,MDtrts, by=c("Blockname" = "BLOCKNA"))

MDtrts<-MDtrts[MDtrts$TREATME == "MD",]
MDtrts<-MDtrts[,colnames(MDtrts) %in% c("ID","STATE","BLOCKNA","hostba","AF_50","DATE_TR","prcp_p_"
                                        ,"rh100","pai","trpctch_v","trpctch_m")]
MDtrts$DATE_TREAT<-as.POSIXct(MDtrts$DATE_TR)
MDtrts$DOY_TREAT<-yday(MDtrts$DATE_TREAT)
MDtrts$AFdiff<-MDtrts$DOY_TREAT - MDtrts$AF_50
MDtrts$stblk<-paste(MDtrts$STATE,MDtrts$BLOCKNA,sep="-")


eval.MD<-right_join(eval.first, MDtrts, by="stblk")
eval.MD<-eval.MD[!duplicated(eval.MD$stblk),]
eval.MD<-eval.MD[eval.MD$TreatYear >= 2015,]
eval.MD<-eval.MD[eval.MD$Treatment == "MD",]
eval.MD<-eval.MD[eval.MD$Dosage != "4.5g",]
eval.MD<-eval.MD[!grepl("Study",eval.MD$Blockname),]

sum(eval.first$TreatYear >= 2015)
sum(eval.first$TreatYear >= 2015 & eval.first$Treatment=="MD")
sum(eval.first$TreatYear >= 2015 & eval.first$Treatment=="Btk")

lm.MD<-betareg(Ta ~ State + Dosage + DistanceKm + Acres + TreatYearFact + AFdiff + prcp_p_,
               data=eval.MD)
#vif(lm.MD)
summary(lm.MD)

lrt(lm.MD, betareg(Ta ~ Dosage + DistanceKm + Acres + TreatYearFact + AFdiff + prcp_p_, data=eval.MD)) #test on state
lrt(lm.MD, betareg(Ta ~ State + DistanceKm + Acres + TreatYearFact + AFdiff + prcp_p_, data=eval.MD)) #test on dosage
lrt(lm.MD, betareg(Ta ~ State + Dosage + Acres + TreatYearFact + AFdiff + prcp_p_, data=eval.MD)) #test on distance
lrt(lm.MD, betareg(Ta ~ State + Dosage + DistanceKm + TreatYearFact + AFdiff + prcp_p_, data=eval.MD)) #test on acres
lrt(lm.MD, betareg(Ta ~ State + Dosage + DistanceKm + AFdiff + prcp_p_, data=eval.MD)) #test on treatment year
lrt(lm.MD, betareg(Ta ~ State + Dosage + DistanceKm + Acres + TreatYearFact + prcp_p_, data=eval.MD)) #test on phenology
lrt(lm.MD, betareg(Ta ~ State + Dosage + DistanceKm + Acres + TreatYearFact + AFdiff, data=eval.MD)) #test on precip post trt

## Recent years, Btk ------------------------------------------------------------------

lethalTrts.poly<-readOGR("./data/trt_lethal_15_18.shp")
lethalTrts<-lethalTrts.poly@data
lethalTrts<-lethalTrts[,colnames(lethalTrts) %in% c("ID","STATE","BLOCKNA","hostba","L2_50","DATE_TR"
                                                    ,"prcp_p_","rh100","pai","trpctch_v","trpctch_m")]
lethalTrts$DATE_TREAT<-as.POSIXct(lethalTrts$DATE_TR)
lethalTrts$DOY_TREAT<-yday(lethalTrts$DATE_TREAT)
lethalTrts$L2diff<-lethalTrts$DOY_TREAT - lethalTrts$L2_50
lethalTrts$stblk<-paste(lethalTrts$STATE,lethalTrts$BLOCKNA,sep="-")


recentLethal<-right_join(eval.all,lethalTrts, by=c("Blockname" = "BLOCKNA"))
## here, test for dose effect, too

eval.Btk<-right_join(eval.first, lethalTrts, by="stblk")
eval.Btk<-eval.Btk[!duplicated(eval.Btk$stblk),]
eval.Btk<-eval.Btk[eval.Btk$TreatYear >= 2015,]
eval.Btk<-eval.Btk[eval.Btk$Treatment == "Btk",]

lm.Btk<-betareg(Ta ~ State + DistanceKm + Acres + TreatYearFact + L2diff + prcp_p_,
               data=eval.Btk)
vif(lm.MD)
summary(lm.Btk)

lrt(lm.Btk, betareg(Ta ~ DistanceKm + Acres + TreatYearFact + L2diff + prcp_p_, data=eval.Btk)) #test on state
lrt(lm.Btk, betareg(Ta ~ State + Acres + TreatYearFact + L2diff + prcp_p_, data=eval.Btk)) #test on distance
lrt(lm.Btk, betareg(Ta ~ State + DistanceKm + TreatYearFact + L2diff + prcp_p_, data=eval.Btk)) #test on acres
lrt(lm.Btk, betareg(Ta ~ State + DistanceKm + Acres + L2diff + prcp_p_, data=eval.Btk)) #test on treatment year
lrt(lm.Btk, betareg(Ta ~ State + DistanceKm + Acres + TreatYearFact + prcp_p_, data=eval.Btk)) #test on L2 diff
lrt(lm.Btk, betareg(Ta ~ State + DistanceKm + Acres + TreatYearFact + L2diff, data=eval.Btk)) #test on prcp



