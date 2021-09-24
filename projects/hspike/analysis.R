# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
# setTimeLimit(100000); setSessionTimeLimit(10000)
library(ggplot2)
library("cowplot")
library("gridExtra")
library(ggpubr)
library(plyr)

# load data: Patients x Units x timewindow
data <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/alldata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")

# prepare data
data$hyplabel     <- factor(data$hyplabel_pow, ordered = TRUE, levels = c("NO_SCORE", "REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))
data$patient      <- factor(data$patient)
data$part         <- factor(data$part,)

# Select only one unit for analyses, effectively removing unit dimension for analyses are restricted to patient x timewindow
data_patient      <- data[data$unit == 1, ]

# to create p-values
detach(package:lmerTest)
library(lmerTest)

# effect of sleepstage and EEG power for IED rate
l1 <- lmer(IEDsum ~ hyplabel + delta + theta + alpha + beta + gamma + (1 | part) + (1 | patient), data_patient); 
summary(l1)

# effect of sleep stage and for firing rate units: strong negative correlation with depth of sleep
data_clean <- na.omit(data) # removing missing data
l1 <- lmer(trialfreq_corrected ~ hyplabel + (1 | part) + (1 | patient) + (1 | unit), data_clean); 
summary(l1)

data_PSG <- data[data$hyplabel != "NO_SCORE", ] # remove "NO_SCORE"
plot(trialfreq_corrected ~ hyplabel , data_PSG)

# effect of sleep stage on CV2/regularity: strong positive correlation with depth of sleep
l1 <- lmer(CV2_trial ~ hyplabel + (1 | part) + (1 | patient) + (1 | unit), data_clean); 
summary(l1)



l1 <- lmer(IEDsum ~ hyplabel + (1 | part) + (1 + IEDsum | patient) + (1 + IEDsum | unit), data); 
summary(l1)




l1 <- lmer(IEDsum ~ delta + theta + beta + gamma + (1 | part) + (1 | patient), data_patient); summary(l1)


data_clean <- na.omit(data)

l1 <- lmer(IEDsum ~ hyplabel + delta + theta + beta + gamma + trialfreq_corrected + (1 | part) + (1 | patient) + (1 | unit), data_clean); summary(l1)



# Remove IEDs first, then predict firing rate based on hypnogram
data_nanan <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/alldata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")


data_noIED <- data[data_clean$IEDsum == 0, ]

l1 <- lmer(trialfreq_corrected ~ hyplabel + delta + theta + beta + gamma + (1 | part) + (1 | patient) + (1 | unit), data_noIED); summary(l1)
















# data$SU <- data$percRPV < 0.25
data$SU <- factor(ifelse(data$percRPV > 0.25, "MUA", "SUA"))
data$SU <- c("MUA","SUA","MUA","MUA","MUA","MUA","SUA","MUA","MUA","MUA","MUA","SUA","MUA","MUA","SUA","MUA","MUA","MUA","SUA","MUA","SUA","SUA","MUA","SUA","MUA")
# data$PI <- factor(data$PI)
data$PI <- factor(ifelse(data$template_tp < 550, "Int", "Pyr"))
# data$PI <- revalue(data$PI, c("Int"="Interneuron", "Pyr"="Pyramidal cell"))




data <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/aurelie/statAH.csv", sep=';', header=TRUE, dec=',', na.strings = " ")
data$P <- data$ï..Patient

# get some p-values
detach(package:lmerTest)
library(lmerTest)
l1 <- lmer(EEG ~ NSE + S100 + (1 | P) + (1 | Time), data); summary(l1)
l1 <- lmer(EEG ~ NSE + S100 + (1 | P), data); summary(l1)

l1 <- lm(EEG ~ NSE + S100, data); summary(l1)





