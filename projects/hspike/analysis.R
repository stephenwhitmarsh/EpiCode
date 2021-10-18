#install.packages("circular")
#install.packages("ggplot2")
#install.packages("units")
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
# data_patient      <- data[data$unit == 1, ]

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







## POLAR ##

library("ggplot2")
library(RColorBrewer)

## plot power values



e <- evalq(aggregate(list(meanAlpha=alpha), list(patient=patient), mean), data)
f <- evalq(aggregate(list(stdAlpha=alpha), list(patient=patient), sd), data)

data <- merge(data, e)
data <- merge(data, f)
data$zAlpha <- (data$alpha-data$meanAlpha)/data$stdAlpha

data$bin <- cut(data$minute, seq(0, 24*60, by = 10))

rm(data_binned)
data_binned$bin <-  as.data.frame(seq(0, 24*60, by = 1))
data_binned$zAlpha <- tapply(data$zAlpha, cut(data$minute, seq(0, 24*60, by = 1)), mean)

data_binned <- as.data.frame(tapply(data$zAlpha, cut(data$minute, seq(0, 24*60, by = 1)), mean), col.names = c("bin", "zAlpha"))

#alpha$bin = cut(data$minute, seq(0, 24*60, by = 1))
bind_rows(alpha)

ggplot(data, aes(x = bin, col = patient)) +
  geom_point(aes( y = zAlpha)) +
  coord_polar(theta = "x")
  
  

ggplot(data, aes(x = bin, col = patient)) +
  geom_point(aes( y = zAlpha)) +
  coord_polar(theta = "x")

ggplot(data, aes(x = bin, col = patient)) +
  stat_summary(aes( y = zAlpha), fun= "median", geom = "point") +
  coord_polar(theta = "x")

ggplot(data=data, aes(fill=patient, col=patient)) +
  #ggtitle("Circadian interictal density & seizure occurance") +
  #scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  #geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
  #geom_segment(aes(x = 0, xend = 0, y = 1, yend = 0.1)) +
  #geom_hline(yintercept = seq(-0.25+0.025, -0.05, by = 0.025), colour = "grey90") +
  geom_histogram(stat = 'identity', aes(y = zAlpha, binwidth=1, center = 0.5)) # make sure binwidth is multiple of 24/60
  
  

ggplot(data=data, aes(x=minute, fill=patient, col=patient)) +
  #ggtitle("Circadian interictal density & seizure occurance") +
  #scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  #geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
  #geom_segment(aes(x = 0, xend = 0, y = 1, yend = 0.1)) +
  #geom_hline(yintercept = seq(-0.25+0.025, -0.05, by = 0.025), colour = "grey90") +
  geom_boxplot(aes(y = alpha, binwidth=0.1, center = 0.05)) +
  coord_polar(theta = "x")
  




                 
                 + # make sure binwidth is multiple of 24/60
 
  
   coord_polar(theta = "x") +
  scale_x_continuous(breaks=seq(0, 21, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
  coord_polar(theta = "x", start = 0, direction = 1) + 
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank(),
  ) + 
  geom_point(size=0.5, data = data_seizures, aes(x = hour, y = -0.05 - as.numeric(patient)*0.025+0.025, fill = patient, col = patient)) +
  ylim(-0.3, 0.85)




# load data: Patients x Units x time window
data_IED <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/hypdata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")

# prepare data
data_IED$patient      <- factor(data_IED$patient, levels = c(1:8))
#data_IED$radian       <- data_IED$theta / (pi*2) * 360
data_IED$hour         <- data_IED$minute / (60)
#data_IED              <- data_IED[data_IED$part > 3, ]
data_IED$part         <- factor(data_IED$part,)

# load data: Patients x Units x time window
data_seizures <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/seizuredata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")

# prepare data
data_seizures$patient  <- factor(data_seizures$patient, levels = c(1:8))
data_seizures$hour     <- data_seizures$minute / (60)
data_seizures$radian   <- units(data_seizures$minute / (60) / 24 * 360, "radians")

# graphics.off()

#display.brewer.all(colorblindFriendly = TRUE)

pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/polar_density_combined.pdf")

ggplot(data=data_IED, aes(x=hour, fill=patient, col=patient)) +
  ggtitle("Circadian interictal density & seizure occurance") +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
  #geom_segment(aes(x = 0, xend = 0, y = 1, yend = 0.1)) +
  geom_hline(yintercept = seq(-0.25+0.025, -0.05, by = 0.025), colour = "grey90") +
  geom_histogram(position = 'stack', aes(y = stat(density)), binwidth=0.1, center = 0.05) + # make sure binwidth is multiple of 24/60
  coord_polar(theta = "x") +
  scale_x_continuous(breaks=seq(0, 21, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
  coord_polar(theta = "x", start = 0, direction = 1) + 
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank(),
        ) + 
  geom_point(size=0.5, data = data_seizures, aes(x = hour, y = -0.05 - as.numeric(patient)*0.025+0.025, fill = patient, col = patient)) +
  ylim(-0.3, 0.85)

dev.off()










ggplot(data=data_IED, aes(x=hour, fill=patient, col=patient)) +
  ggtitle("Circadian interictal density & seizure occurance") +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
  #geom_segment(aes(x = 0, xend = 0, y = 1, yend = 0.1)) +
  geom_hline(yintercept = seq(-0.25+0.025, -0.05, by = 0.025), colour = "grey90") +
  geom_histogram(position = 'stack', aes(y = stat(density)), binwidth=0.1, center = 0.05) + # make sure binwidth is multiple of 24/60
  coord_polar(theta = "x") +
  scale_x_continuous(breaks=seq(0, 21, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
  coord_polar(theta = "x", start = 0, direction = 1) + 
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank(),
  ) + 
  geom_point(size=0.5, data = data_seizures, aes(x = hour, y = -0.05 - as.numeric(patient)*0.025+0.025, fill = patient, col = patient)) +
  ylim(-0.3, 0.85)













library(circular)
library(units)

plot.circular(data_seizures$radian, pch = 16, cex = 1, stack = TRUE,
     axes = TRUE, start.sep=0, sep = 0.025, shrink = 1,
     bins = 23, ticks = TRUE, tcl = 0.025, tcl.text = 0.125,
     col = NULL, tol = 0.04, uin = NULL,
     xlim = c(-1, 1), ylim = c(-1, 1), digits = 2, units = "rad",
     template = NULL, zero = 0.1, rotation = NULL,
     main = NULL, sub=NULL, xlab = "", ylab = "")


plot(sin, -pi, 2*pi) # see ?plot.function

pi_rad <- as_units(pi, "radians")



pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/polar_density_per patient.pdf")

ggplot(data=data, aes(x=hour, fill=patient, col=patient)) +
  ggtitle("Circadian interictal density per patient") +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_histogram(position = 'stack', aes(y = stat(density)), binwidth=0.2, center = 0.1) + # make sure binwidth is multiple of 24/60
  coord_polar(theta = "x") +
  geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
  scale_x_continuous(breaks=seq(0, 21, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
  coord_polar(theta = "x", start = 0, direction = 1) + 
  theme_bw() +
  
  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank(),
        strip.background = element_blank(),
  ) + facet_wrap(~patient)
dev.off()






# get ylim and xlim from plot
#layer_scales(p)$y$range$range
#layer_scales(p)$x$range$range





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



