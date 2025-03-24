library(lme4)
Moeteke=Moeteke[Moeteke$Post_Intervention==1,]
Moeteke[Moeteke$ParticipantID=="p060",]$StudyGroup=2 # correct typo (see email correspondence April 12, 2024)

dummPrimary=ifelse(Moeteke$HospitalLevelofHealthcare == 1, 1, 0)
dummSecondary=ifelse(Moeteke$HospitalLevelofHealthcare == 2, 1, 0)
dummTertiary=ifelse(Moeteke$HospitalLevelofHealthcare == 3, 1, 0)

Moeteke=cbind(Moeteke,dummPrimary,dummSecondary,dummTertiary)


######################################################################################################################
### models per outcome separately
######################################################################################################################

############ empty model

model1=lmer(Change_in_Knowledge~1+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Attitude~1+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Confidence~1+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Practice~1+(1|LGAofPractice),data=Moeteke)
summary(model1)


############ with study condition but no confounders

model1=lmer(Change_in_Knowledge~1+StudyGroup+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Attitude~1+StudyGroup+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Confidence~1+StudyGroup+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Practice~1+StudyGroup+(1|LGAofPractice),data=Moeteke)
summary(model1)



############ with study condition and confounders

model1=lmer(Change_in_Knowledge~1+StudyGroup+Age+YearsofPractice+factor(HospitalLevelofHealthcare)+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Attitude~1+StudyGroup+Age+YearsofPractice+factor(HospitalLevelofHealthcare)+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Confidence~1+StudyGroup+Age+YearsofPractice+factor(HospitalLevelofHealthcare)+(1|LGAofPractice),data=Moeteke)
summary(model1)


model1=lmer(Change_in_Practice~1+StudyGroup+Age+YearsofPractice+factor(HospitalLevelofHealthcare)+(1|LGAofPractice),data=Moeteke)
summary(model1)


######################################################################################################################
############ create long format of data set for multivariate response model NOT FINISHED YET
######################################################################################################################

Moeteke.covars=cbind(Moeteke$LGAofPractice, Moeteke$StudyGroup, Moeteke$Age, Moeteke$YearsofPractice, Moeteke$HospitalLevelofHealthcare)
ID=rep(seq(1,261),time=4)

dummy.K=c(rep(1,261),rep(0,261),rep(0,261),rep(0,261))
dummy.A=c(rep(0,261),rep(1,261),rep(0,261),rep(0,261))
dummy.C=c(rep(0,261),rep(0,261),rep(1,261),rep(0,261))
dummy.P=c(rep(0,261),rep(0,261),rep(0,261),rep(1,261))

outcome=c(Moeteke$Change_in_Knowledge, Moeteke$Change_in_Attitude, Moeteke$Change_in_Confidence, Moeteke$Change_in_Practice)
outcome.var=rep(c("Knowledge", "Attitude", "Confidence", "Practice"),each=261)

Moeteke.long=cbind(ID, rbind(Moeteke.covars,Moeteke.covars,Moeteke.covars,Moeteke.covars), outcome.var, dummy.K, dummy.A, dummy.C, dummy.P,outcome)
colnames(Moeteke.long)=c("id","cluster","condition","age","years","level", "ouctome.var", "dummyK","dummyA","dummyC","dummyP","outcome")
Moeteke.long[order(Moeteke.long[,1],decreasing=FALSE),]

Moeteke.long=as.data.frame((Moeteke.long))
Moeteke.long$outcome=as.numeric(Moeteke.long$outcom)
lmer(outcome~dummy.K+dummy.A+dummy.C+dummy.P-1+(1|cluster/id),data=Moeteke.long)

library(glmmTMB)
model1=glmmTMB(outcome~dummy.K+dummy.A+dummy.C+dummy.P-1+(as.factor(outcome.var)-1|cluster/id),data=Moeteke.long)
summary(model1)

model2=glmmTMB(outcome~dummy.K+dummy.A+dummy.C+dummy.P+dummy.K*condition+dummy.A*condition+dummy.C*condition+dummy.P*condition-1-condition+(as.factor(outcome.var)-1|cluster/id),data=Moeteke.long)
summary(model2)

model2=glmmTMB(outcome~dummy.K+dummy.A+dummy.C+dummy.P+dummy.K*condition+dummy.A*condition+dummy.C*condition+dummy.P*condition-1-condition+(as.factor(outcome.var)-1|cluster),data=Moeteke.long)
summary(model2)


######################################################################################################################
############ lavaan: only the two outcomes with significant effects from lmer
######################################################################################################################

library(lavaan)
#https://lavaan.ugent.be/tutorial/multilevel.html

############ lavaan: only the two outcomes with significant effects ##################################################################################
mod1 <- "
level: 1
Change_in_Knowledge ~~ Change_in_Confidence

level: 2
Change_in_Knowledge ~~ Change_in_Confidence
"

fit1 <- sem(model = mod1, data = Moeteke, cluster = "LGAofPractice")
summary(fit1, fit.measures = TRUE)
standardizedsolution(fit1)


############ lavaan: two outcomes with study group as predictor ##################################################################################
mod2 <- "
level: 1
Change_in_Knowledge ~~ Change_in_Confidence


level: 2
Change_in_Knowledge ~~ Change_in_Confidence
Change_in_Knowledge + Change_in_Confidence ~ StudyGroup
"

fit2 <- sem(model = mod2, data = Moeteke, cluster = "LGAofPractice")
summary(fit2, fit.measures = TRUE)
standardizedsolution(fit2)


############ lavaan: two outcomes with study group and confounders ##################################################################################
mod3 <- "
level: 1
#Change_in_Knowledge ~~ Change_in_Confidence
Change_in_Knowledge + Change_in_Confidence ~ Age+YearsofPractice+dummPrimary+dummSecondary

level: 2
#Change_in_Knowledge ~~ Change_in_Confidence
Change_in_Knowledge + Change_in_Confidence ~ StudyGroup
"

fit3 <- sem(model = mod3, data = Moeteke, cluster = "LGAofPractice")
summary(fit3, fit.measures = TRUE)
standardizedsolution(fit3)

#fit6 <- sem(model = mod1, data = Moeteke, cluster = "LGAofPractice",verbose = TRUE, optim.method = "em", em.iter.max = 2000,em.fx.tol = 1e-08, em.dx.tol = 1e-04)
summary(fit3, fit.measures = TRUE)

#########################################################################################################################
############ lavaan: all four outcomes ##################################################################################
#########################################################################################################################
mod4 <- "
level: 1
Change_in_Knowledge ~~ Change_in_Attitude + Change_in_Confidence + Change_in_Practice
Change_in_Attitude ~~ Change_in_Confidence + Change_in_Practice
Change_in_Confidence ~~ Change_in_Practice

level: 2
Change_in_Knowledge ~~ Change_in_Attitude + Change_in_Confidence + Change_in_Practice
Change_in_Attitude ~~ Change_in_Confidence + Change_in_Practice
Change_in_Confidence ~~ Change_in_Practice
"
fit4 <- sem(model = mod4, data = Moeteke, cluster = "LGAofPractice")
#fit4 <- sem(model = mod4, data = Moeteke, cluster = "LGAofPractice",verbose = TRUE, optim.method = "em", em.iter.max = 2000,em.fx.tol = 1e-08, em.dx.tol = 1e-04)
summary(fit4, fit.measures = TRUE)
standardizedsolution(fit4)


############ lavaan: all four outcomes with zero between variances for confidence and practice ##################################################################################
mod5 <- "
level: 1
Change_in_Knowledge ~~ Change_in_Attitude + Change_in_Confidence + Change_in_Practice
Change_in_Attitude ~~ Change_in_Confidence + Change_in_Practice
Change_in_Confidence ~~ Change_in_Practice

level: 2
Change_in_Knowledge ~~  Change_in_Confidence 
"

fit5 <- sem(model = mod5, data = Moeteke, cluster = "LGAofPractice")
summary(fit5, fit.measures = TRUE)
#fit5 <- sem(model = mod5, data = Moeteke, cluster = "LGAofPractice",verbose = TRUE, optim.method = "em", em.iter.max = 2000,em.fx.tol = 1e-08, em.dx.tol = 1e-04)
standardizedsolution(fit5)


############ lavaan: all four outcomes with zero between variances for confidence and practice and study group as predictor ##################################################################################
mod6 <- "
level: 1
Change_in_Knowledge ~~ Change_in_Attitude + Change_in_Confidence + Change_in_Practice
Change_in_Attitude ~~ Change_in_Confidence + Change_in_Practice
Change_in_Confidence ~~ Change_in_Practice

level: 2
Change_in_Knowledge ~~ 0*Change_in_Attitude + Change_in_Confidence + 0*Change_in_Practice
Change_in_Attitude ~~ 0*Change_in_Confidence + 0*Change_in_Practice
Change_in_Confidence ~~ 0*Change_in_Practice 


Change_in_Attitude ~~ 0*Change_in_Attitude
Change_in_Practice ~~ 0*Change_in_Practice

Change_in_Knowledge + Change_in_Confidence +Change_in_Attitude +Change_in_Practice  ~ StudyGroup
"

fit6 <- sem(model = mod6, data = Moeteke, cluster = "LGAofPractice")
summary(fit6, fit.measures = TRUE)
#fit6 <- sem(model = mod6, data = Moeteke, cluster = "LGAofPractice",verbose = TRUE, optim.method = "em", em.iter.max = 2000,em.fx.tol = 1e-08, em.dx.tol = 1e-04)
standardizedsolution(fit6)




