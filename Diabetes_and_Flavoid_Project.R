library(readr)
library(meta)
library(metafor)
library(ggplot2)
library(ggrepel)
library(forcats)
library(dplyr)
library(grid)
library(gridExtra)
library(revtools)

## Screening the database and fetch matching literature using machine learning

# Fetch literature in database
data_all <- read_bibliography(c("2001-2010pubmed.nbib", "2001-2010WOS.bib"))
class(data_all)
# screen_duplicates
matches <- find_duplicates(data_all)
data_unique <- extract_unique_references(data_all, matches)
# Export matching literature
data_export <- subset(data_unique, select = c("abstract", "title"))
write.csv(data_export, file = "data_cleaned.csv", row.names = FALSE)


## Glucose Assessment

# Data import
glucose <- read.csv("glucose_final.csv")
str(glucose)
colnames(glucose)[1]<- "Author"
m.hksj <- glucose(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = glucose,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")

# effect sizes
m.hksj
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()

# eggers test and funnel plot 
tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 
funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=200, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = glucose$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = glucose$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = glucose$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = glucose$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = glucose$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup_health.tiff", res=300, height=200, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = glucose$Health_status)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot_es.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 
tiff("sensitivity_plot_i2.tiff", res=300, height=150, width=255, units="mm")
i2 <- plot(inf.analysis, "I2")
dev.off() 
metareg(m.hksj,flavonoid)


## 2h_PPG Assessment

# Data import
PPG_2h <- read.csv("2h-PPG.csv")
str(PPG_2h)
colnames(PPG_2h)[1]<- "Author"
m.hksj <- PPG_2h(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = PPG_2h,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj

# effect sizes
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()

# eggers test and funnel plot 
tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 
funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=80, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = PPG_2h$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = PPG_2h$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = PPG_2h$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = PPG_2h$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = PPG_2h$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup.tiff", res=300, height=230, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = PPG_2h$flavonoid)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 
metareg(m.hksj,flavonoid)


## HbA1C Assessment

# Data import
HbA1C <- read.csv("HbA1c.csv")
str(HbA1C)
colnames(HbA1C)[1]<- "Author"

m.hksj <- HbA1C(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = HbA1C,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj

# effect sizes
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()




tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 


funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

trimfill(m.hksj)
m.hksj$TE.random
m.hksj.trimfill<-trimfill(m.hksj)
funnel(m.hksj.trimfill,xlab = "Hedges' g")
m2 <- eggers.test(x = m.hksj.trimfill)
sink("eggers.test_corrected.txt")
print(m2)
sink()

tiff("funnel_plot_corrected.tiff", res=300, height=180, width=195, units="mm")
fu2 <- funnel(m.hksj.trimfill,xlab = "Hedges' g",studlab = TRUE)
dev.off() 


# eggers test and funnel plot 


tiff("forest_plot.tiff", res=300, height=110, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 
# forest plot


sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HbA1C$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")
# dosage sub_analysis

sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HbA1C$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")
# duration sub_analysis

sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HbA1C$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")
# flavonoid sub_analysis

sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HbA1C$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")
# form sub_analysis

sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HbA1C$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")
# health sub_analysis


tiff("Forest_subgroup_health.tiff", res=300, height=150, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = HbA1C$Health_status)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   
# subgroup plot


inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot_es.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 

tiff("sensitivity_plot_i2.tiff", res=300, height=150, width=255, units="mm")
i2 <- plot(inf.analysis, "I2")
dev.off() 
# sensitivity analysis

metareg(m.hksj,flavonoid)


## HDL Assessment

HDL <- read.csv("HDL.csv")
str(HDL)
colnames(HDL)[1]<- "Author"

m.hksj <- HDL(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = HDL,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj


sink("results_radom_effect.txt")
print(m.hksj)
sink()
# effect sizes


m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()
# remove outliers

# eggers test and funnel plot
tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 


funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=150, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HDL$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HDL$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HDL$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HDL$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HDL$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup.tiff", res=300, height=300, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = HDL$flavonoid)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 
metareg(m.hksj,flavonoid)


## HOMA-IR Assessment

# Import Data
HOMA_IR <- read.csv("HOMA-IR.csv")
str(HOMA_IR)
colnames(HOMA_IR)[1]<- "Author"

m.hksj <- HOMA_IR(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = HOMA_IR,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj

# effect sizes
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()

# eggers test and funnel plot 
tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 


funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=140, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_IR$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_IR$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_IR$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_IR$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_IR$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup_flavonoid.tiff", res=300, height=230, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = HOMA_IR$flavonoid)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   
tiff("Forest_subgroup_form.tiff", res=300, height=180, width=195, units="mm")

sgame2 <- subgroup.analysis.mixed.effects(x = m.hksj,
                                          subgroups = HOMA_IR$Form)
forest(sgame2,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 
tiff("sensitivity_plot_i2.tiff", res=300, height=150, width=255, units="mm")
i2 <- plot(inf.analysis, "I2")
dev.off()
metareg(m.hksj,flavonoid)


## HOMA-Beta Assessment

# Import Data
HOMA_Beta <- read.csv("HOMA-beta.csv")
str(HOMA_Beta)
colnames(HOMA_Beta)[1]<- "Author"
m.hksj <- HOMA_Beta(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = HOMA_Beta,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj

# effect sizes
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()

# eggers test and funnel plot 
tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 


funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=80, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_Beta$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_Beta$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_Beta$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_Beta$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = HOMA_Beta$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup_flavonoid.tiff", res=300, height=230, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = HOMA_Beta$flavonoid)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   
tiff("Forest_subgroup_form.tiff", res=300, height=180, width=195, units="mm")

sgame2 <- subgroup.analysis.mixed.effects(x = m.hksj,
                                          subgroups = HOMA_Beta$Form)
forest(sgame2,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off()
metareg(m.hksj,flavonoid)


## Insulin Assessment

# Import Data
Insulin <- read.csv("Insulin.csv")
str(Insulin)
colnames(Insulin)[1]<- "Author"

m.hksj <- Insulin(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = Insulin,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj

# effect sizes
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()

# eggers test and funnel plot 
tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 


funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=150, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = Insulin$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = Insulin$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = Insulin$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = Insulin$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = Insulin$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup_health_mixed.tiff", res=300, height=200, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = Insulin$flavonoid)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()    

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot_es.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 
tiff("sensitivity_plot_i2.tiff", res=300, height=150, width=255, units="mm")
i2 <- plot(inf.analysis, "I2")
dev.off()  
metareg(m.hksj,flavonoid)


## LDL Assessment

# Import Data
LDL <- read.csv("LDL.csv")
str(LDL)
colnames(LDL)[1]<- "Author"

m.hksj <- LDL(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = LDL,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj

# effect sizes
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()

# eggers test and funnel plot 
tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 


funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=150, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = LDL$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = LDL$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = LDL$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = LDL$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = LDL$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup.tiff", res=300, height=300, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = LDL$flavonoid)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 
metareg(m.hksj,flavonoid)

## TC Assessment

# Import Data
TC <- read.csv("TC.csv")
str(TC)
colnames(TC)[1]<- "Author"
m.hksj <- TC(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = TC,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj

# effect sizes
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()

# eggers test and funnel plot 
tiff("funnel_plot.tiff", res=300, height=180, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 
funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=150, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TC$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TC$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TC$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TC$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TC$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup.tiff", res=300, height=300, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = TC$flavonoid)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 
metareg(m.hksj,flavonoid)


## TG Assessment

# Import Data
TG <- read.csv("TG.csv")
str(TG)
colnames(TG)[1]<- "Author"



m.hksj <- TG(Ne,
                   Me,
                   Se,
                   Nc,
                   Mc,
                   Sc,
                   data = TG,
                   studlab = paste(Author),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD")
m.hksj

# effect sizes
sink("results_radom_effect.txt")
print(m.hksj)
sink()

# remove outliers
m_remove <- find.outliers(m.hksj)
sink("results_remove_outliers.txt")
print(m_remove)
sink()

# eggers test and funnel plot 
tiff("funnel_plot.tiff", res=300, height=220, width=195, units="mm")
fu <- funnel(m.hksj,xlab = "Hedges' g",studlab = TRUE)
dev.off() 
funnel(m.hksj, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

m1 <- eggers.test(x = m.hksj)
sink("eggers.test.txt")
print(m1)
sink()

# forest plot
tiff("forest_plot.tiff", res=300, height=150, width=255, units="mm")
fo <- forest(m.hksj,
             col.diamond = "black",
             leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
             lab.e = "Intervention",
             text.random = "Overall effect",
             comb.random = TRUE,
             print.Q = TRUE,
             prediction = FALSE,
             col.predict = "black")
dev.off() 

# dosage sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TG$Dosage)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "dosage.csv")
write.csv(sub_analysis$subgroup.analysis.results, "dosage2.csv")

# duration sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TG$Duration)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "duration.csv")
write.csv(sub_analysis$subgroup.analysis.results, "duration2.csv")

# flavonoid sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TG$flavonoid)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "flavonoid.csv")
write.csv(sub_analysis$subgroup.analysis.results, "flavonoid2.csv")

# form sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TG$Form)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "form.csv")
write.csv(sub_analysis$subgroup.analysis.results, "form2.csv")

# health sub_analysis
sub_analysis <- subgroup.analysis.mixed.effects(x = m.hksj,
                                                subgroups = TG$Health_status)
sub_analysis
write.csv(sub_analysis$within.subgroup.results, "health.csv")
write.csv(sub_analysis$subgroup.analysis.results, "health2.csv")

# subgroup plot
tiff("Forest_subgroup.tiff", res=300, height=300, width=195, units="mm")
sgame <- subgroup.analysis.mixed.effects(x = m.hksj,
                                         subgroups = TG$flavonoid)
forest(sgame,
       col.diamond = "black",
       prediction = FALSE)               
dev.off()   

# sensitivity analysis
inf.analysis <- InfluenceAnalysis(x = m.hksj,
                                  random = TRUE)
summary(inf.analysis)
tiff("sensitivity_plot.tiff", res=300, height=150, width=255, units="mm")
se <- plot(inf.analysis, "es")
dev.off() 
metareg(m.hksj,flavonoid)