## Package and options
library("MedSurvey")
options(prompt = "R> ", continue = "+  ", width = 70,
  useFancyQuotes = FALSE)


## Data and related information
MedData
R <- 160
wgtnames <- paste("repwgt", seq(0,R,by=1), sep="")
mwgtname=wgtnames[1]
repwgtnames=wgtnames[2:(R+1)]


## Sepcify the model 1
model1 <- ' 
+ 	numcg ~ u0*1 + c*workban + b1*sp_adltban + b2*sp_kidsban
+ 	sp_adltban ~ u1*1 + a1*workban
+ 	sp_kidsban ~ u2*1 + a2*workban
+ 	sp_adltban ~~ sp_kidsban
+ 	a1b1 := a1*b1
+ 	a2b2 := a2*b2
+ 	total := c + (a1*b1) + (a2*b2)
+ '


## Fit the model 1
fit.BRR <- med.fit.BRR(model=model1, data=MedData, mwgtname=mwgtname,
+                   repwgtnames=repwgtnames, fayfactor=0.5)



## View the summary results of the mediation analysis
med.summary(fit=fit.BRR, med.eff=c('a1b1' , 'a2b2'))


## View the model fit statistics adjusted for the complex survey designs
model0 <- ' 
+      numcg ~ u0*1 + c*workban + b1*sp_adltban + b2*sp_kidsban
+      sp_adltban ~ u1*1 + a1*workban
+      sp_kidsban ~ u2*1 + a2*workban
+      a1b1 := a1*b1
+      a2b2 := a2*b2
+      total := c + (a1*b1) + (a2*b2)
+      '
require(lavaan)
fit <- lavaan::sem(model=model0, data=MedData, estimator='ML', test='standard')
chisq.BRR(model0,fit,MedData,mwgtname, repwgtnames)

## Read the data for model 2
R> PisaMed <- as.data.frame(read.csv(
+ url("https://raw.githubusercontent.com/YujiaoMai/MedSurvey/master/PisaMed.csv")))
R> R <- 80
R> repwgtnames <- paste("W_FSTURWT", seq(1,R,by=1), sep="")
R> mwgtname <- "W_FSTUWT"

## Specify the model 2
R> model2 <- '
+ StuMtv =~ ST119Q01NA+ST119Q02NA+ST119Q03NA+ST119Q04NA+ST119Q05NA
+ StuAnxt =~ ST118Q01NA+ST118Q02NA+ST118Q03NA+ST118Q04NA+ST118Q05NA
+ ParenSpt =~ ST123Q01NA+ST123Q02NA+ST123Q03NA+ST123Q04NA
+ Math =~ PV1MATH+PV2MATH+PV3MATH+PV4MATH+PV5MATH+PV6MATH+PV7MATH
+ +PV8MATH+PV9MATH+PV10MATH
+ Math ~ c0*ParenSpt + b1*StuMtv + AGE + SEX + StuAnxt
+ StuMtv ~ a1*ParenSpt + AGE + SEX
+ StuMtv ~~ StuAnxt
+ a1b1 := a1*b1
+ '

## Fit the model 2
R> fit.BRR2 <- med.fit.BRR(model=model2, data=PisaMed, mwgtname=mwgtname,
+ repwgtnames=repwgtnames, fayfactor=0.5)

## Check the summary results of the mediation analysis
R> med.summary(fit=fit.BRR2, med.eff = c('a1b1'))


