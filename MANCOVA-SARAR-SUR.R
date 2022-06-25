#Load libraries-------------------------------------------------------------------------
library(rgdal) #shp files management
library(sp) #spatial data management
library(spdep) #spatial dependence: weighting schemes
library(sphet) #spatial autoregressive models with heteroscedascticity
library(spatialreg) #spatial regression models
library(jtools) #for robust errors
library(walrus) #for robust errors


#Set session preferences-------------------------------------------------------------------------
# Working directory
setwd("C:/Users/Norbert/OneDrive/A02. Florin - Specializare sectoriala/02. Model")

#Sets the number of digits to 3, and scientific notation to 7
options(digits=3, scipen=3)

# Report Robust errors, VIFs and three digits
set_summ_defaults (digits = 3, robust=TRUE, vifs=TRUE)


#Import data-------------------------------------------------------------------------
modelanova <- read.csv("C:/Users/Norbert/OneDrive/A02. Florin - Specializare sectoriala/02. Model/model anova.csv")
summary (modelanova)
names(modelanova)

#Equation-------------------------------------------------------------------------
eq1 = D_Ln_PopAdj ~ Cluster
eq2 = D_Ln_PopAdj ~ Cluster + Region
eq3 = D_Ln_PopAdj ~ Cluster + Region  + LocalityType
eq4 = D_Ln_PopAdj ~ Cluster + Region  + LocalityType +  Cluster*LocalityType
eq5 = D_Ln_PopAdj ~ Cluster + Region  + LocalityType +  Cluster*LocalityType + D_Ln_Emp
eq6 = D_Ln_PopAdj ~ Cluster + Region  + LocalityType +  Cluster*LocalityType + D_Ln_Emp + D_Ln_Turnover
eq7 = D_Ln_PopAdj ~ Cluster + Region  + LocalityType +  Cluster*LocalityType + D_Ln_Emp + D_Ln_Turnover + D_Ln_PIT 
eq8 = D_Ln_PopAdj ~ Cluster + Region  + LocalityType +  Cluster*LocalityType + D_Ln_Emp + D_Ln_Turnover + D_Ln_PIT + CoreEconomy

eq4a = D_Ln_PopAdj ~ Cluster + Region  + LocalityType +  Cluster*LocalityType + CoreEconomy


#Spatial data preparation-------------------------------------------------------------------------
uat = readOGR(dsn = ".", layer = "UAT-uri") #citeste shp UAT-uri
names(uat) #afiseaza variabilele
uat$SIRUTA=as.character(uat$SIRUTA) # recodeaza SIRUTA in dimensiune
spplot(uat,"SIRUTA") #afiseaza harta dupa siruta

# Lists of spatial weights from the shape file DIFF------------------------------------------------------------
colnames(modelanova)[1] <- "SIRUTA" #make sure the key has the same name in both tables
spat.anova <- merge(uat, modelanova, by = "SIRUTA", all.x=F, all.y=F) #merge files
queen.diff = poly2nb(spat.anova) #neighbors list with single shared boundary point 
queen.listw = nb2listw(queen.diff) #converts neighbors list to listwise object type

# Non-spatial models OLS and ANOVA-------------------------------------------------------------------------
m4a <- lm(eq4a, modelanova)
print(anova(m4a), digits=10)

m7 <- lm(eq7, modelanova)
print(anova(m8), digits=10)

m8 <- lm(eq8, modelanova)
print(anova(m8), digits=10)

# Spatial models OLS and ANOVA-----------------------------------------------------------------------------------------------
m4aols <- lm(eq4a, spat.anova)
print(anova(m4aols), digits=10)

m7ols <- lm(eq7, spat.anova)
print(anova(m7ols), digits=10)

m8ols <- lm(eq8, spat.anova)
print(anova(m8ols), digits=10)

#Moran test & Lagrange multiplier test----------------------------------------------------------------------------------------
lm.morantest(m4aols,queen.listw)
lm.LMtests(m4aols,queen.listw,test="all")

lm.morantest(m7ols,queen.listw)
lm.LMtests(m7ols,queen.listw,test="all")

lm.morantest(m7ols,queen.listw)
lm.LMtests(m7ols,queen.listw,test="all")

lm.morantest(m8ols,queen.listw)
lm.LMtests(m8ols,queen.listw,test="all")


#OLS F-estimation------------------------------------------------------------

m1 <-lm(eq1, modelanova)
m2 <-lm(eq2, modelanova)
m3 <-lm(eq3, modelanova)
m4 <-lm(eq4, modelanova)
m5 <-lm(eq5, modelanova)
m6 <-lm(eq6, modelanova)
m7 <-lm(eq7, modelanova)
m8 <-lm(eq8, modelanova)

sstotal <-(sum((modelanova$D_Ln_PopAdj - mean(modelanova$D_Ln_PopAdj))^2))

#Model1
ssresidual1 <-sum(m1$residuals^2)
ssregression.Cluster <- sstotal-ssresidual1

#Model2
ssresidual2 <-sum(m2$residuals^2)
ssregression2 <- sstotal-ssresidual2
ssregression.Region <- ssregression2-ssregression1

#Model3
ssresidual3 <-sum(m3$residuals^2)
ssregression3 <- sstotal-ssresidual3
ssregression.LocalityType <- ssregression3-ssregression2
print(ssregression3, digits=8)

#Model4
ssresidual4 <-sum(m4$residuals^2)
ssregression4 <- sstotal-ssresidual4
ssregression.Cluster_LocalityType <- ssregression4-ssregression3

#Model5
ssresidual5 <-sum(m5$residuals^2)
ssregression5 <- sstotal-ssresidual5
ssregression.D_Ln_Emp <- ssregression5-ssregression4

#Model6
ssresidual6 <-sum(m6$residuals^2)
ssregression6 <- sstotal-ssresidual6
ssregression.D_Ln_Turnover <- ssregression6-ssregression5

#Model7
ssresidual7 <-sum(m7$residuals^2)
ssregression7 <- sstotal-ssresidual7
ssregression.D_Ln_PIT <- ssregression7-ssregression6

#Model8
ssresidual8 <-sum(m8$residuals^2)
ssregression8 <- sstotal-ssresidual8
ssregression.CoreEconomy <- ssregression8-ssregression7

print(anova(m8),digits=4)
ssregression.Region
ssregression.LocalityType
ssregression.Cluster
ssregression.CoreEconomy
ssregression.Cluster_LocalityType
ssregression.D_Ln_Emp
ssregression.D_Ln_PIT
ssregression.D_Ln_Turnover

#ANOVA SUR-SARAR------------------------------------------------------------


#Mode1
ms1 <-spatialreg::sacsarlm(formula = eq1, data = spat.anova, listw = queen.listw)
ms2 <-spatialreg::sacsarlm(formula = eq2, data = spat.anova, listw = queen.listw)
ms3 <-spatialreg::sacsarlm(formula = eq3, data = spat.anova, listw = queen.listw)
ms4 <-spatialreg::sacsarlm(formula = eq4, data = spat.anova, listw = queen.listw)
ms5 <-spatialreg::sacsarlm(formula = eq5, data = spat.anova, listw = queen.listw)
ms6 <-spatialreg::sacsarlm(formula = eq6, data = spat.anova, listw = queen.listw)
ms7 <-spatialreg::sacsarlm(formula = eq7, data = spat.anova, listw = queen.listw)
ms8 <-spatialreg::sacsarlm(formula = eq8, data = spat.anova, listw = queen.listw)

#Total residuals
s.sstotal <-(sum((spat.anova$D_Ln_PopAdj - mean(spat.anova$D_Ln_PopAdj))^2))

#Model1
s.ssresidual1 <-sum(ms1$residuals^2)
s.ssregression1 <- s.sstotal-s.ssresidual1
s.ssregression.Cluster <- s.sstotal-s.ssresidual1

#Model2
s.ssresidual2 <-sum(ms2$residuals^2)
s.ssregression2 <- s.sstotal-s.ssresidual2
s.ssregression.Region <- s.ssregression2-s.ssregression1

#Model3
s.ssresidual3 <-sum(ms3$residuals^2)
s.ssregression3 <- s.sstotal-s.ssresidual3
s.ssregression.LocalityType <- s.ssregression3-s.ssregression2
print(s.ssregression3, digits=8)

#Model4
s.ssresidual4 <-sum(ms4$residuals^2)
s.ssregression4 <- s.sstotal-s.ssresidual4
s.ssregression.CoreEconomy <- s.ssregression4-s.ssregression3

#Model5
s.ssresidual5 <-sum(ms5$residuals^2)
s.ssregression5 <- s.sstotal-s.ssresidual5
s.ssregression.Cluster_LocalityType <- s.ssregression5-s.ssregression4

#Model6
s.ssresidual6 <-sum(ms6$residuals^2)
s.ssregression6 <- s.sstotal-s.ssresidual6
s.ssregression.D_Ln_Emp <- s.ssregression6-s.ssregression5

#Model7
s.ssresidual7 <-sum(ms7$residuals^2)
s.ssregression7 <- s.sstotal-s.ssresidual7
s.ssregression.D_Ln_Turnover <- s.ssregression7-s.ssregression6


#Model8
s.ssresidual8 <-sum(ms8$residuals^2)
s.ssregression8 <- sstotal-s.ssresidual8
s.ssregression.D_Ln_PIT <- s.ssregression8-s.ssregression7

print(s.ssresidual1, digits=8)
print(s.ssresidual2, digits=8)
print(s.ssresidual3, digits=8)
print(s.ssresidual4, digits=8)
print(s.ssresidual5, digits=8)
print(s.ssresidual6, digits=8)
print(s.ssresidual7, digits=8)
print(s.ssresidual8, digits=8)

print(s.ssregression.Region, digits=8)
print(s.ssregression.LocalityType, digits=8)
print(s.ssregression.Cluster, digits=8)
print(s.ssregression.CoreEconomy, digits=8)
print(s.ssregression.Cluster_LocalityType, digits=8)
print(s.ssregression.D_Ln_Emp, digits=8)
print(s.ssregression.D_Ln_PIT, digits=8)
print(s.ssregression.D_Ln_Turnover, digits=8)

summary(m7)
summ(m7)

summary(ms7, Nagelkerke=TRUE, Hausman=TRUE)
summary(ms8, Nagelkerke=TRUE, Hausman=TRUE)



#Compare models-----------------------------------------------------------------------------------------------------------------
m4a <-lm(eq4a,modelanova)
summary(m7, Nagelkerke=TRUE)
print(anova(m7), digits=8)

library(performance)
compare_performance(m4a, m7, m8, ms7, ms8, rank = TRUE)

#Save Yhat-----------------------------------------------------------------------------------------------------------------

#Model 7 ANOVA-LM
m7.pred <- data.frame (modelanova$SIRUTA, m7$fitted.values)
colnames(m7.pred)[1] <- "SIRUTA" #first column renamed to SIRUTA
colnames(m7.pred)[2] <- "Yhat.M7" #first column renamed to SIRUTA
write.csv (m7.pred, file="C:/Users/Norbert/OneDrive/A02. Florin - Specializare sectoriala/02. Model/m7.yhat.csv")


#Model 7 SUR-SARAR
ms7.pred <- data.frame (spat.anova$SIRUTA, ms7$fitted.values)
colnames(ms7.pred)[1] <- "SIRUTA" #first column renamed to SIRUTA
colnames(ms7.pred)[2] <- "Yhat.MS7" #first column renamed to SIRUTA
write.csv (ms7.pred, file="C:/Users/Norbert/OneDrive/A02. Florin - Specializare sectoriala/02. Model/ms7.yhat.csv")



#Post-hoc Tests-----------------------------------------------------------------------------------------------------------------
# Nonspatial
library(multcompView)
m8aov <- aov(eq8, modelanova)
summary(m8aov)
summ (m8)
m8T <- TukeyHSD(m8aov, 'Cluster', conf.level=0.95)
m8T <- TukeyHSD(m8aov, 'Region', conf.level=0.95)
m8T <- TukeyHSD(m8aov, 'LocalityType', conf.level=0.95)
m8T <- TukeyHSD(m8aov, 'CoreEconomy', conf.level=0.95)
m8T <- TukeyHSD(m8aov, 'Cluster:LocalityType', conf.level=0.95)

print(m8T)


# Spatial
eqC <- ms7$fitted.values ~ Cluster + Region + LocalityType + LocalityType*Cluster
mC <- aov (eqC, spat.anova)
TukeyCluster7 <- TukeyHSD(mC, 'Cluster', conf.level=0.95)
TukeyCluster7 <- TukeyHSD(mC, 'Region', conf.level=0.95)
TukeyCluster7 <- TukeyHSD(mC, 'LocalityType', conf.level=0.95)
TukeyCluster7 <- TukeyHSD(mC, 'Cluster:LocalityType', conf.level=0.95)

print(TukeyCluster7)
TukeyCluster <- TukeyHSD(m8, 'Cluster', conf.level=0.95)
TukeyRegion <- TukeyHSD(m8, 'Region', conf.level=0.95)
TukeyLocalityType <- TukeyHSD(m8, 'LocalityType', conf.level=0.95)
TukeyCoreEconomy <- TukeyHSD(m8, 'CoreEconomy', conf.level=0.95)
TukeyCoreEconomy <- TukeyHSD(m8, 'CoreEconomy', conf.level=0.95)

