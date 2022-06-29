################################################################################
##  Naylor& Kawano - Mudskipper locomotion on complex media: contact variables ##
##                                                             June 13, 2022   ##
################################################################################



rm(list = ls()) #clear environment 

# 1) Load packages and data -----------------------------------------------
 
library(stats)
library(lme4)
library(performance)
library(emmeans)
library(MuMIn)


## read in 'Naylor&Kawano_mudskippercontactvariables.csv' 
contact.data <- read.csv(choose.files(caption = "select contact variables csv"))

#substrates
#FS = sand
#GM = mud
#SL = solid


#check classification of each column
sapply(contact.data,class)
#modify classifications as needed
contact.data$Incline <- as.character(contact.data$Incline)
contact.data$Fish.Weight.g <- as.numeric(contact.data$Fish.Weight.g)



# 2) Find best model  ------------------------------------------------------

#model averaging approach:
#compare linear mixed-effects models (LMMs) with different covariate combinations
#base structure: outcome variable ~ substrate*incline + (1|Indiv)
#REML = FALSE (maximum likelihood)
#covariates: Ave.speed.units.s (average stride speed; continuous), tail.body.use (binary), 
  # Premovement (binary)

##count occurrences of tail use, premovement per environmental condition
# contact.data %>% count(tailbodyprop, Substrate, Incline)
# contact.data %>% count(Premovement, Substrate, Incline)


#### pectoral fin duty factor ####

#a) Run LMMs

#no additional fixed terms/covariates
DFpect0 <- lmer(pectoralDF ~  Substrate * Incline + (1|Indiv),
                contact.data, REML=FALSE)

#one additional fixed term/covariate
DFpect1a <- lmer(pectoralDF ~  Substrate * Incline + Ave.speed.units.s +
                (1|Indiv), contact.data, REML=FALSE)

DFpect1b <- lmer(pectoralDF ~  Substrate * Incline + tail.body.use + 
                (1|Indiv), contact.data, REML=FALSE)

DFpect1c <- lmer(pectoralDF ~  Substrate * Incline + Premovement + 
                (1|Indiv), contact.data, REML=FALSE)

#two additional fixed terms/covariates
DFpect2a <- lmer(pectoralDF ~  Substrate * Incline + Ave.speed.units.s + 
                   tail.body.use + (1|Indiv), contact.data, REML=FALSE)

DFpect2b <- lmer(pectoralDF ~  Substrate * Incline + Ave.speed.units.s + 
                   Premovement + (1|Indiv), contact.data, REML=FALSE)

DFpect2c <- lmer(pectoralDF ~  Substrate * Incline + tail.body.use + 
                   Premovement + (1|Indiv), contact.data, REML=FALSE)

#three additional fixed terms/covariates
DFpect3 <- lmer(pectoralDF ~  Substrate * Incline + Ave.speed.units.s + 
                  tail.body.use + Premovement + (1|Indiv), contact.data, REML=FALSE)



#b) rank models by AICc score
model.sel(DFpect0, DFpect1a, DFpect1b, DFpect1c, DFpect2a, DFpect2b, DFpect2c, DFpect3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

anova(DFpect0, DFpect1c) #DFpect0 and DFpect1c have delta AICc< 2; select DFpect0



#### pelvic fin duty factor ####

#a) Run LMMs

#no additional fixed terms/covariates
DFpelvic0 <- lmer(pelvicDF.total ~  Substrate * Incline + (1|Indiv),
                contact.data, REML=FALSE)

#one additional fixed term/covariate
DFpelvic1a <- lmer(pelvicDF.total ~  Substrate * Incline + Ave.speed.units.s +
                   (1|Indiv), contact.data, REML=FALSE)

DFpelvic1b <- lmer(pelvicDF.total ~  Substrate * Incline + tail.body.use + 
                   (1|Indiv), contact.data, REML=FALSE)

DFpelvic1c <- lmer(pelvicDF.total ~  Substrate * Incline + Premovement + 
                   (1|Indiv), contact.data, REML=FALSE)

#two additional fixed terms/covariates
DFpelvic2a <- lmer(pelvicDF.total ~  Substrate * Incline + Ave.speed.units.s + 
                   tail.body.use + (1|Indiv), contact.data, REML=FALSE)

DFpelvic2b <- lmer(pelvicDF.total ~  Substrate * Incline + Ave.speed.units.s + 
                   Premovement + (1|Indiv), contact.data, REML=FALSE)

DFpelvic2c <- lmer(pelvicDF.total ~  Substrate * Incline + tail.body.use + 
                   Premovement + (1|Indiv), contact.data, REML=FALSE)

#three additional fixed terms/covariates
DFpelvic3 <- lmer(pelvicDF.total ~  Substrate * Incline + Ave.speed.units.s + 
                  tail.body.use + Premovement + (1|Indiv), contact.data, REML=FALSE)



#b) rank models by AICc score
model.sel(DFpelvic0, DFpelvic1a, DFpelvic1b, DFpelvic1c, DFpelvic2a, DFpelvic2b, DFpelvic2c, DFpelvic3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

anova(DFpelvic0, DFpelvic2c) #DFpelvic0 and DFpelvic2c have delta AICc< 2; select DFpelvic2c



#### proportion of pectoral fin propulsion ####

#a) Run LMMs

#no additional fixed terms/covariates
PectProp0 <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + (1|Indiv),
                contact.data, REML=FALSE)

#one additional fixed term/covariate
PectProp1a <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + Ave.speed.units.s +
                   (1|Indiv), contact.data, REML=FALSE)

PectProp1b <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + tail.body.use + 
                   (1|Indiv), contact.data, REML=FALSE)

PectProp1c <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + Premovement + 
                   (1|Indiv), contact.data, REML=FALSE)

#two additional fixed terms/covariates
PectProp2a <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + Ave.speed.units.s + 
                   tail.body.use + (1|Indiv), contact.data, REML=FALSE)

PectProp2b <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + Ave.speed.units.s + 
                   Premovement + (1|Indiv), contact.data, REML=FALSE)

PectProp2c <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + tail.body.use + 
                   Premovement + (1|Indiv), contact.data, REML=FALSE)

#three additional fixed terms/covariates
PectProp3 <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + Ave.speed.units.s + 
                  tail.body.use + Premovement + (1|Indiv), contact.data, REML=FALSE)



#b) rank models by AICc score
model.sel(PectProp0, PectProp1a, PectProp1b, PectProp1c, PectProp2a, PectProp2b, PectProp2c, PectProp3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

anova(PectProp0, PectProp1c) #PectProp0 and PectProp1c have delta AICc< 2; select PectProp0


#### maximum pectoral fin spread ####

#a) Run LMMs

#no additional fixed terms/covariates
FinSpread0 <- lmer(pectoral.spread.mm ~  Substrate * Incline + (1|Indiv),
                  contact.data, REML=FALSE)

#one additional fixed term/covariate
FinSpread1a <- lmer(pectoral.spread.mm ~  Substrate * Incline + Ave.speed.units.s +
                     (1|Indiv), contact.data, REML=FALSE)

FinSpread1b <- lmer(pectoral.spread.mm ~  Substrate * Incline + tail.body.use + 
                     (1|Indiv), contact.data, REML=FALSE)

FinSpread1c <- lmer(pectoral.spread.mm ~  Substrate * Incline + Premovement + 
                     (1|Indiv), contact.data, REML=FALSE)

#two additional fixed terms/covariates
FinSpread2a <- lmer(pectoral.spread.mm ~  Substrate * Incline + Ave.speed.units.s + 
                     tail.body.use + (1|Indiv), contact.data, REML=FALSE)

FinSpread2b <- lmer(pectoral.spread.mm ~  Substrate * Incline + Ave.speed.units.s + 
                     Premovement + (1|Indiv), contact.data, REML=FALSE)

FinSpread2c <- lmer(pectoral.spread.mm ~  Substrate * Incline + tail.body.use + 
                     Premovement + (1|Indiv), contact.data, REML=FALSE)

#three additional fixed terms/covariates
FinSpread3 <- lmer(pectoral.spread.mm ~  Substrate * Incline + Ave.speed.units.s + 
                    tail.body.use + Premovement + (1|Indiv), contact.data, REML=FALSE)



#b) rank models by AICc score
model.sel(FinSpread0, FinSpread1a, FinSpread1b, FinSpread1c, FinSpread2a, FinSpread2b, FinSpread2c, FinSpread3)


#c) use likelihood ratio tests (LRTs) to compare all top ranking models 
#(delta AICc < 2); determine if added complexity improves fit
#more reduced model is retained unless Pr(>Chi-sq) is less than 0.05

anova(FinSpread0, FinSpread1b) #FinSpread0 and FinSpread1b have delta AICc< 2; select FinSpread0





# 3) Calculate point and interval estimates  ------------------------------------------------------


#### pectoral fin duty factor ####

#a) run best LMM with REML=TRUE
  #less biased estimation of random effects and better for smaller sample sizes 
    #(Bates et al. 2015)

DFpect0 <- lmer(pectoralDF ~  Substrate * Incline + (1|Indiv),
                contact.data, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
  #*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(DFpect0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(DFpect0))
qqline(residuals(DFpect0), col="blue")
hist(residuals(DFpect0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(DFpect0, ~ Substrate*Incline)


#d) calculate R-squared as indicator of model fit 
    #(i.e., how much variance is explained by the model?)
    #R2conditional: fixed and random effects 
    #R2marginal: fixed effects only

r2_nakagawa(DFpect0)



#### pelvic fin duty factor ####

#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

DFpelvic2c <- lmer(pelvicDF.total ~  Substrate * Incline + tail.body.use + 
                     Premovement + (1|Indiv),
                contact.data, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(DFpelvic2c, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(DFpelvic2c))
qqline(residuals(DFpelvic2c), col="blue")
hist(residuals(DFpelvic2c), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(DFpelvic2c, ~ Substrate*Incline + tail.body.use + Premovement)
         

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(DFpelvic2c)


#### proportion of pectoral fin propulsion ####


#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

PectProp0 <- lmer(prop.pectoral.propulsion ~  Substrate * Incline + (1|Indiv),
                contact.data, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(PectProp0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(PectProp0))
qqline(residuals(PectProp0), col="blue")
hist(residuals(PectProp0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(PectProp0, ~ Substrate*Incline)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(PectProp0)


#### maximum pectoral fin spread ####


#a) run best LMM with REML=TRUE
#less biased estimation of random effects and better for smaller sample sizes 
#(Bates et al. 2015)

FinSpread0 <- lmer(pectoral.spread.mm ~  Substrate * Incline + (1|Indiv),
                contact.data, REML=TRUE)

#b) check for outliers and inspect for normal distribution* of model residuals
#*but see LeBeau et al. 2018 on normality for LMMs

check_outliers(FinSpread0, method= c("cook","pareto"), threshold = NULL)
qqnorm(residuals(FinSpread0))
qqline(residuals(FinSpread0), col="blue")
hist(residuals(FinSpread0), freq=FALSE, breaks=7, col=rainbow(22))

#c) calculate estimated marginal means (EMMs) and interval estimates

emmeans(FinSpread0, ~ Substrate*Incline)

#d) calculate R-squared as indicator of model fit 
#(i.e., how much variance is explained by the model?)
#R2conditional: fixed and random effects 
#R2marginal: fixed effects only

r2_nakagawa(FinSpread0)
