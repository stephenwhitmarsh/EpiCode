#install.packages("circular")
#install.packages("ggplot2")
#install.packages("units")
#install.packages("reshape")
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
# setTimeLimit(100000); setSessionTimeLimit(10000)

# install.packages("rsleep")
# library("rsleep")
# e <- data.frame(begin = as.POSIXlt(c(1536967800,1536967830,1536967860),origin = "1970-01-01"))
# e$end <- as.POSIXlt(c(1536967830,1536967860,1536967890), origin = "1970-01-01")
# e$event = c("N3","N3","REM")
# plot_hypnogram(e)



library(ggplot2)
library("cowplot")
library("gridExtra")
library(ggpubr)
library(plyr)
library(reshape)


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
library("dplyr")

## plot power values PER HOUR

data$bin    <- as.integer(cut(data$minute, seq(0, 24*60, by = 60)))
data_binned <- setNames(aggregate(data$alpha, c(list(data$patient), list(data$bin)), mean), c("patient", "bin", "alpha"))
data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Nalpha = alpha/max(alpha))) 

ggplot(data=data_binned, aes(x = bin, y = alpha, fill=patient, col=patient)) +
  ggtitle("Circadian alpha power") +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_bar(stat="identity", width = 1) +
  geom_vline(xintercept = seq(0.5, 21.5, by = 3), colour = "grey90")  + 
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank(),
  ) + 
  scale_x_continuous(breaks=seq(0.5, 21.5, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
  coord_polar(theta = "x", start = 0)


## plot power values (use bindivision for hour divisions)
bindivision = 1
data$bin    <- as.integer(cut(data$minute, seq(0, 24*60, by = 60/bindivision)))
y1 <- setNames(aggregate(data$alpha, c(list(data$patient), list(data$bin)), median), c("patient", "bin", "alpha"))
y2 <- setNames(aggregate(data$theta, c(list(data$patient), list(data$bin)), median), c("patient", "bin", "theta"))
y3 <- setNames(aggregate(data$beta,  c(list(data$patient), list(data$bin)), median), c("patient", "bin", "beta"))
y4 <- setNames(aggregate(data$delta, c(list(data$patient), list(data$bin)), median), c("patient", "bin", "delta"))
y1sd <- setNames(aggregate(data$alpha, c(list(data$patient), list(data$bin)), sd), c("patient", "bin", "SDalpha"))
y2sd <- setNames(aggregate(data$theta, c(list(data$patient), list(data$bin)), sd), c("patient", "bin", "SDtheta"))
y3sd <- setNames(aggregate(data$beta,  c(list(data$patient), list(data$bin)), sd), c("patient", "bin", "SDbeta"))
y4sd <- setNames(aggregate(data$delta, c(list(data$patient), list(data$bin)), sd), c("patient", "bin", "SDdelta"))
data_binned <- merge(y1,y2)
data_binned <- merge(data_binned,y3)
data_binned <- merge(data_binned,y4)
data_binned <- merge(data_binned,y1sd)
data_binned <- merge(data_binned,y2sd)
data_binned <- merge(data_binned,y3sd)
data_binned <- merge(data_binned,y4sd)
data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Nalpha = (alpha-min(alpha)) / max(alpha-min(alpha))))
data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Ntheta = (theta-min(theta)) / max(theta-min(theta))))  
data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Ndelta = (delta-min(delta)) / max(delta-min(delta))))
data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Nbeta  = (beta-min(beta))   / max(beta-min(beta))))
data_binned_melt <-melt(data_binned[c("Ndelta","Ntheta","Nalpha","Nbeta","patient","bin")], id=c("patient","bin"))
data_binned_melt$variable <- factor(data_binned_melt$variable, labels = c("\u03B4 power (1-4Hz)","\u03B8 power (5-7Hz)","\u03B1 power (8-14Hz)","\u03B2 power (15-25Hz)"))
data_binned_melt$patient  <- factor(data_binned_melt$patient, levels = c(8:1))

pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/polar_power_combined.pdf")

  ggplot(data=data_binned_melt, aes(x = bin, y = value, fill=patient, col=patient)) +
    ggtitle("Normalized circadian LFP power") +
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    geom_vline(xintercept = seq(0.5, 21.5 * bindivision, by = 3 * bindivision), colour = "grey90")  + 
    geom_hline(yintercept = seq(1, 8, by = 1), colour = "grey90")  + 
    geom_bar(stat="identity", width = 1) +
    theme_bw() +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          panel.grid  = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
    ) + 
    scale_x_continuous(breaks=seq(0.5, 21.5*bindivision, by = 3 * bindivision), 
                       labels = c("0" = "00:00", "3" = "3:00", "6" = "6:00", "9" = "9:00", 
                                  "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
    coord_polar(theta = "x", start = 0) +
    ylim(-1, 7) +
    facet_wrap(~variable)

dev.off()

library(circular)

data$radian = data$minute / (24*60) * pi * 2 

stats_circ <- median.circular(data$radian)


# load data: Patients x Units x time window
data_IED <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/hypdata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")

# prepare data
data_IED$patient      <- factor(data_IED$patient, levels = c(8:1))
data_IED$hour         <- data_IED$minute / (60)
data_IED$part         <- factor(data_IED$part,)

# load data: Patients x Units x time window
data_seizures <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/seizuredata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")

# prepare data
data_seizures$patient  <- factor(data_seizures$patient, levels = c(8:1))
data_seizures$hour     <- data_seizures$minute / (60)

# graphics.off()

# display.brewer.all(colorblindFriendly = TRUE)

pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/polar_density_combined.pdf")

  ggplot(data=data_IED, aes(x=hour, fill=patient, col=patient)) +
    ggtitle("Circadian interictal density & seizure occurance") +
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
    #geom_segment(aes(x = 0, xend = 0, y = 1, yend = 0.1)) +
    geom_hline(yintercept = seq(-0.25+0.025, -0.05, by = 0.025), colour = "grey90") +
    geom_histogram(position = 'stack', aes(y = stat(density)), binwidth=0.1, center = 0.05) + # make sure binwidth is multiple of 24/60
    coord_polar(theta = "x") +
    scale_x_continuous(breaks=seq(0, 21, by = 3), labels = c("0" = "00:00", "3" = "3:00", "6" = "6:00", "9" = "9:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
    coord_polar(theta = "x", start = 0, direction = 1) + 
    theme_bw() +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          panel.grid  = element_blank(),
          ) + 
    geom_point(size=1, data = data_seizures, aes(x = hour, y = -0.05 - as.numeric(patient)*0.025+0.025, fill = patient, col = patient)) +
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



