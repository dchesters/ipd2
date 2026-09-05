 

 ###########################################################
 # read and process table

t1 <- read.table("integrated_result.tsv", header=T, row.names=1,na.strings = "NA") 

# model stepping will crash due to na's, remove last column which is full of these, and rm rows with them
t1 <- t1[,1:21];t1 <- na.omit(t1)

# make multigene-multigene the baseline
t1[,"class_pair"] <- with(t1, factor(class_pair, levels = c("GG","CG","GM","GO","GP","MM","MO","OO")))

head(t1)


 ###########################################################
 # remove regression outliers

library(lme4)
model_lmer <- lmer(log(t1[,"TreeDistance"]+1)~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first")+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), REML = FALSE)
stud_resid <- rstudent(model_lmer)
cooks_dist <- cooks.distance(model_lmer)
flag <- cooks_dist > 4/length(stud_resid) & abs(stud_resid) > 2
t1 <- t1[!flag, ]
# remake model
model_lmer <- lmer(log(t1[,"TreeDistance"]+1)~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first")+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), REML = FALSE)


 ###########################################################
 # compare mixed effects to basic lm, needs fitting to same table
 # tree dist should not need sum overlap, nor need normalizing

model_lm <- lm(log(t1[,"TreeDistance"]+1)~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first"))

anova(model_lmer, model_lm)
#           Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
#model_lm   15 -3151.4 -3067.6 1590.7  -3181.4                             
#model_lmer 17 -3804.9 -3710.0 1919.5  -3838.9 657.55      2  < 2.2e-16 ***
# significance mean random effects matter, and lm could be inaccurate.


 ###########################################################
 # check lmer assumptions

 # normality
qqnorm(resid(model_lmer))

# Heteroskedasticity, equal variances. 
# ncvTest runs on lm(), workaround for lmer
library(car); 	# qqPlot, ncvTest
aux_model <- lm(I(residuals(model_lmer)^2) ~ fitted(model_lmer))
ncvTest(aux_model)

# Multicollinearity Check,
# difficulty getting vif running on lmer, so used vif for lm()


 ###########################################################
 # table needs more specific formatting or lmer stepping doesnt work
 # 'object is not a matrix' error, then you need to remake model

t1$year_diff <- sqrt(abs(t1[,"year"]-t1[,"year2"]))
t1$log_TreeDistance <- log(t1[,"TreeDistance"]+1)
t1$log_sumloci <- log(t1[,"sum_loci"])
t1$log_morph <- rank(t1[,"morph_chars"], ties.method = "first")
t1$normalized_matchingsplitdist <- log((t1[,"MatchingSplitDistance"]/t1[,"sum_overlap"])+1)
t1$normalized_MASTSize <- log((t1[,"MASTSize"]/t1[,"sum_overlap"])+1)


 ###########################################################
 # lmerTest gives p values

library(lmerTest)
model_lmer <- lmer(log_TreeDistance~t1[,"class_pair"]+t1[,"year"]+year_diff+t1[,"rank_numeric"]+log_sumloci+t1[,"prop_overlap"]+log_morph+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), data = t1, REML = FALSE)

aic_model_lmer <- lmerTest::step(model_lmer)
summary(get_model(aic_model_lmer))

 ###########################################################
 # repeat for MatchingSplitDistance

t1 <- read.table("integrated_result", header=T, row.names=1,na.strings = "NA") 
t1 <- t1[,1:21];t1 <- na.omit(t1)
t1[,"class_pair"] <- with(t1, factor(class_pair, levels = c("GG","CG","GM","GO","GP","MM","MO","OO")))
library(lme4)
model_lmer <- lmer(t1[,"MatchingSplitDistance"]/t1[,"sum_overlap"]~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first")+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), REML = FALSE)
stud_resid <- rstudent(model_lmer)
cooks_dist <- cooks.distance(model_lmer)
flag <- cooks_dist > 4/length(stud_resid) & abs(stud_resid) > 2
t1 <- t1[!flag, ]
model_lmer <- lmer(t1[,"MatchingSplitDistance"]/t1[,"sum_overlap"]~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first")+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), REML = FALSE)
t1$year_diff <- sqrt(abs(t1[,"year"]-t1[,"year2"]))
t1$log_sumloci <- log(t1[,"sum_loci"])
t1$log_morph <- rank(t1[,"morph_chars"], ties.method = "first")
t1$normalized_matchingsplitdist <- log((t1[,"MatchingSplitDistance"]/t1[,"sum_overlap"])+1)
t1$normalized_MASTSize <- log((t1[,"MASTSize"]/t1[,"sum_overlap"])+1)
library(lmerTest)
model_lmer <- lmer(normalized_matchingsplitdist~t1[,"class_pair"]+t1[,"year"]+year_diff+t1[,"rank_numeric"]+log_sumloci+t1[,"prop_overlap"]+log_morph+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), data = t1, REML = FALSE)
aic_model_lmer <- lmerTest::step(model_lmer)
summary(get_model(aic_model_lmer))

 ###########################################################
 # repeat for MASTSize

t1 <- read.table("integrated_result", header=T, row.names=1,na.strings = "NA") 
t1 <- t1[,1:21];t1 <- na.omit(t1)
t1[,"class_pair"] <- with(t1, factor(class_pair, levels = c("GG","CG","GM","GO","GP","MM","MO","OO")))
library(lme4)
model_lmer <- lmer(t1[,"MASTSize"]/t1[,"sum_overlap"]~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first")+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), REML = FALSE)
stud_resid <- rstudent(model_lmer)
cooks_dist <- cooks.distance(model_lmer)
flag <- cooks_dist > 4/length(stud_resid) & abs(stud_resid) > 2
t1 <- t1[!flag, ]
model_lmer <- lmer(t1[,"MASTSize"]/t1[,"sum_overlap"]~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first")+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), REML = FALSE)
t1$year_diff <- sqrt(abs(t1[,"year"]-t1[,"year2"]))
t1$log_sumloci <- log(t1[,"sum_loci"])
t1$log_morph <- rank(t1[,"morph_chars"], ties.method = "first")
t1$normalized_matchingsplitdist <- log((t1[,"MatchingSplitDistance"]/t1[,"sum_overlap"])+1)
t1$normalized_MASTSize <- log((t1[,"MASTSize"]/t1[,"sum_overlap"])+1)
library(lmerTest)
model_lmer <- lmer(normalized_MASTSize~t1[,"class_pair"]+t1[,"year"]+year_diff+t1[,"rank_numeric"]+log_sumloci+t1[,"prop_overlap"]+log_morph+ (1 | t1[, "nwk1"]) + (1 | t1[, "nwk2"]), data = t1, REML = FALSE)
aic_model_lmer <- lmerTest::step(model_lmer)
summary(get_model(aic_model_lmer))

 ###########################################################
 # simple linear model on full dataset, despite non independence. 

t1 <- read.table("integrated_result", header=T, row.names=1,na.strings = "NA") 
t1 <- t1[,1:21];t1 <- na.omit(t1)
t1[,"class_pair"] <- with(t1, factor(class_pair, levels = c("GG","CG","GM","GO","GP","MM","MO","OO")))
model_lm <- lm(log(t1[,"TreeDistance"]+1)~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first"))
# model_lm <- lm(t1[,"MatchingSplitDistance"]/t1[,"sum_overlap"]~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first"))
# model_lm <- lm(log((t1[,"MASTSize"]/t1[,"sum_overlap"])+1)~t1[,"class_pair"]+t1[,"year"]+sqrt(abs(t1[,"year"]-t1[,"year2"]))+log(t1[,"rank_numeric"])+log(t1[,"sum_loci"])+t1[,"prop_overlap"]+rank(t1[,"morph_chars"], ties.method = "first"))


stud_resid <- rstudent(model_lm)
cooks_dist <- cooks.distance(model_lm)
flag <- cooks_dist > 4/length(stud_resid) & abs(stud_resid) > 2
t1 <- t1[!flag, ]
# remake model (above)

# test for normality of the model
library(car); 	# qqPlot, ncvTest
qqPlot(model_lm)
# Heteroskedasticity Check , for equal variances. 
ncvTest(model_lm)
# p-value indicating that this model does not have a problem of unequal variances.
# Multicollinearity Check
vif(model_lm)
# AIC
library(MASS); 	# stepAIC
aic_model<-stepAIC(model_lm)
print(summary(aic_model))



