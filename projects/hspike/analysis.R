#install.packages("circular")
#install.packages("ggplot2")
#install.packages("units")
#install.packages("reshape2")
#install.packages("circlize")
#install.packages("ggthemes")
#install.packages("lemon")
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#setTimeLimit(100000); setSessionTimeLimit(10000)
#devtools::install_github("strengejacke/sjPlot")

library(ggplot2)
#library("cowplot")
#library("gridExtra")
library(ggpubr)
library(plyr)
library(reshape)
library(RColorBrewer)
#library("ggthemes")
#library(lemon)
library("sjPlot")
library(emmeans)

#############################
# LFP power per sleep stage #
#############################

# prepare data
data          <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/power_table_long.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data          <- data[!data$part > 3, ] # hypnogram is only scored on first three nights
data$hyplabel <- factor(data$hyplabel, ordered = TRUE, levels = c("NO_SCORE", "REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))
data$band     <- factor(data$band, ordered = TRUE, levels = c("delta", "theta", "alpha", "beta"))
data$patient  <- factor(data$patient, levels = c(8:1))
data$part     <- factor(data$part)

# scale data over windows
m             <- setNames(aggregate(data$power, by = list(data$patient, data$band, data$hyplabel), median), c("patient", "band", "hyplabel", "power_median"))
data          <- merge(data, m)
m             <- setNames(aggregate(data$power, by = list(data$patient, data$band), mean), c("patient", "band", "power_avg"))
s             <- setNames(aggregate(data$power, by = list(data$patient, data$band), sd),   c("patient", "band", "power_sd"))
data          <- merge(data, m)
data          <- merge(data, s)
data$Zpower   <- (data$power-data$power_avg)/data$power_sd

# clean data
data          <- data[!data$hyplabel == "NO_SCORE", ]
data          <- data[!is.na(data$hyplabel), ]
data$hyplabel <- factor(data$hyplabel, ordered = TRUE, levels = c("REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))

# rename levels for plotting
levels(data$hyplabel)[levels(data$hyplabel) == "PHASE_1"] <- "S1"
levels(data$hyplabel)[levels(data$hyplabel) == "PHASE_2"] <- "S2"
levels(data$hyplabel)[levels(data$hyplabel) == "PHASE_3"] <- "S3"

# plot patients separately
# pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/boxplot_power_combined.pdf")
pdf(file="D:/Dropbox/Apps/Overleaf/Hspike/boxplot_power_combined.pdf")

  ggplot(data=data, aes(x = hyplabel, y = Zpower, fill = patient)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    theme_bw() +
    coord_cartesian(ylim = c(-1.2, 3)) +
    theme(
      panel.border = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.y = element_blank(),
          panel.grid.major.x  = element_blank(),
          #axis.title.x = element_blank(),
          #axis.title.y = element_blank(),
          strip.background = element_blank()) +
    xlab("Sleep stage") + ylab("Power (z-scored)") +
    facet_wrap(~band, scales = "fixed", labeller = labeller(band = c("delta" = "Delta (1-4Hz)",
                                                                     "theta" = "Theta (5-7Hz)",
                                                                     "alpha" = "Alpha (8-14Hz)",
                                                                     "beta"  = "Beta  (15-24Hz)")))
dev.off()
  

# plot average per patient
#pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/boxplot_power_patient.pdf")
pdf(file="D:/Dropbox/Apps/Overleaf/Hspike/boxplot_power_patient.pdf")

  ggplot(data=data, aes(x = hyplabel, y = Zpower)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    theme_bw() +
    coord_cartesian(ylim = c(-1.2, 3)) +
    theme(
      panel.border = element_blank(),
      legend.key = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.y = element_blank(),
      panel.grid.major.x  = element_blank(),
      #axis.title.x = element_blank(),
      #axis.title.y = element_blank(),
      strip.background = element_blank()) +
    xlab("Sleep stage") + ylab("Power (z-scored)") +
    facet_wrap(~band, scales = "fixed", labeller = labeller(band = c("delta" = "Delta (1-4Hz)",
                                                                     "theta" = "Theta (5-7Hz)",
                                                                     "alpha" = "Alpha (8-14Hz)",
                                                                     "beta"  = "Beta  (15-24Hz)")))
dev.off()

# plot average per patient
# pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/boxplot_power_band.pdf")
pdf(file="D:/Dropbox/Apps/Overleaf/Hspike/boxplot_power_band.pdf")

  ggplot(data=data, aes(x = hyplabel, y = Zpower, fill = band)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    theme_bw() +
    coord_cartesian(ylim = c(-1.2, 3)) +
    theme(
      panel.border = element_blank(),
      legend.key = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.y = element_blank(),
      panel.grid.major.x  = element_blank(),
      #axis.title.x = element_blank(),
      #axis.title.y = element_blank(),
      strip.background = element_blank()) +
    xlab("Sleep stage") + ylab("Power (z-scored)") 
dev.off()


#######################
# LFP power circadian #
#######################

# prepare data
data          <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/power_table_long.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data$hyplabel <- factor(data$hyplabel, ordered = TRUE, levels = c("NO_SCORE", "REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))
data$band     <- factor(data$band, ordered = TRUE, levels = c("delta", "theta", "alpha", "beta"))
data$patient  <- factor(data$patient, levels = c(8:1))
data$part     <- factor(data$part)

# bin for polar representation
data$bin    <- as.integer(cut(data$minute, seq(0, 24*60, by = 60)))
data_binned <- setNames(aggregate(data$power, c(list(data$patient), list(data$bin), list(data$band)), mean), c("patient", "bin", "band", "power"))
data_binned <- as.data.frame(data_binned %>% group_by(patient, band) %>% mutate(Npower = power/max(power))) 

# plot
# pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/polar_power_band.pdf")
pdf(file="D:/Dropbox/Apps/Overleaf/Hspike/polar_power_band.pdf")

  ggplot(data=data_binned, aes(x = bin, y = Npower, fill=patient, col=patient)) +
    ggtitle("Circadian power") +
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    geom_vline(xintercept = seq(0.5, 21.5, by = 3), colour = "grey90") + 
    geom_hline(yintercept = seq(1, 7, by = 1), colour = "grey90") + 
    geom_bar(stat="identity", width = 1) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      legend.key = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major.x  = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.background = element_blank()) +
    scale_x_continuous(breaks=seq(0.5, 21.5, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
    coord_polar(theta = "x", start = 0) +
    ylim(0, 7) +
    facet_wrap(~band, scales = "fixed", labeller = labeller(band = c("delta" = "Delta (1-4Hz)",
                                                                     "theta" = "Theta (5-7Hz)",
                                                                     "alpha" = "Alpha (8-14Hz)",
                                                                     "beta"  = "Beta  (15-24Hz)")))
dev.off()


#################################
# IED rate & seizures circadian #
#################################

data_IED                <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/IED_table_PSG.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_IED[data_IED$hyplabel == "PRE_SLEEP", ] = "AWAKE"
data_IED[data_IED$hyplabel == "POST_SLEEP", ] = "AWAKE"
data_IED$hyplabel       <- factor(data_IED$hyplabel, ordered = TRUE, levels = c("REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))
data_IED$patient        <- factor(data_IED$patient, levels = c(8:1))
data_IED$part           <- factor(data_IED$part)
data_IED$marker         <- factor(data_IED$marker)
data_IED$hour           <- data_IED$minute / (60)

# load data: Patients x Units x time window
data_seizures           <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/seizuredata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")

# prepare data
data_seizures$patient   <- factor(data_seizures$patient, levels = c(8:1))
data_seizures$hour      <- data_seizures$minute / (60)

# plot
pdf(file="D:/Dropbox/Apps/Overleaf/Hspike/polar_density_combined.pdf")

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
          legend.key   = element_blank(),
          axis.ticks   = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid   = element_blank(),
    ) + 
    geom_point(size=1, data = data_seizures, aes(x = hour, y = -0.05 - as.numeric(patient)*0.025+0.025, fill = patient, col = patient)) + ylim(-0.3, 0.85)

dev.off()


##########################
# IED rate & sleep stage #
##########################

# prepare data
data_IED                <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/IED_table_PSG.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_IED                <- data_IED[!data_IED$hyplabel == "NO_SCORE", ]
data_IED[data_IED$hyplabel == "PRE_SLEEP", ] = "AWAKE"
data_IED[data_IED$hyplabel == "POST_SLEEP", ] = "AWAKE"
data_IED$hyplabel       <- factor(data_IED$hyplabel, ordered = TRUE, levels = c("REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))
data_IED$patient        <- factor(data_IED$patient, levels = c(8:1))
data_IED$part           <- factor(data_IED$part)
data_IED$marker         <- factor(data_IED$marker)
data_IED$hour           <- data_IED$minute / (60)

# create table time spend in sleep stages

# normalize by time spend in 


ggplot(data=data_IED, aes(x=hyplabel, fill=patient, col=patient)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  #geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
  #geom_segment(aes(x = 0, xend = 0, y = 1, yend = 0.1)) +
  #geom_hline(yintercept = seq(-0.25+0.025, -0.05, by = 0.025), colour = "grey90") +
  geom_bar(position = 'stack', aes(stat = "count")) # make sure binwidth is multiple of 24/60
  


## STOPPED HERE




######################
# IED STATISTICS 
######################

# prepare data
data          <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/power_table_wide.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data          <- data[!data$hyplabel == "NO_SCORE", ]
data[data$hyplabel == "PRE_SLEEP", ] = "AWAKE"
data[data$hyplabel == "POST_SLEEP", ] = "AWAKE"
    
data$hyplabel <- factor(data$hyplabel, ordered = TRUE, levels = c("REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))
data$patient  <- factor(data$patient, levels = c(8:1))
data$part     <- factor(data$part)

# normalize 
y1   <- setNames(aggregate(data$alpha, by = c(list(data$patient)), mean), c("patient", "Malpha"))
y2   <- setNames(aggregate(data$theta, by = c(list(data$patient)), mean), c("patient", "Mtheta"))
y3   <- setNames(aggregate(data$beta,  by = c(list(data$patient)), mean), c("patient", "Mbeta"))
y4   <- setNames(aggregate(data$delta, by = c(list(data$patient)), mean), c("patient", "Mdelta"))
y1sd <- setNames(aggregate(data$alpha, by = c(list(data$patient)), sd),   c("patient", "SDalpha"))
y2sd <- setNames(aggregate(data$theta, by = c(list(data$patient)), sd),   c("patient", "SDtheta"))
y3sd <- setNames(aggregate(data$beta,  by = c(list(data$patient)), sd),   c("patient", "SDbeta"))
y4sd <- setNames(aggregate(data$delta, by = c(list(data$patient)), sd),   c("patient", "SDdelta"))

data <- merge(data,y1)
data <- merge(data,y2)
data <- merge(data,y3)
data <- merge(data,y4)
data <- merge(data,y1sd)
data <- merge(data,y2sd)
data <- merge(data,y3sd)
data <- merge(data,y4sd)

data$Zdelta <- (data$delta-data$Mdelta)/data$SDdelta
data$Ztheta <- (data$theta-data$Mtheta)/data$SDtheta
data$Zalpha <- (data$alpha-data$Malpha)/data$SDalpha
data$Zbeta  <- (data$beta -data$Mbeta)/data$SDbeta

# to create p-values
detach(package:lmerTest)
library(lmerTest)
library(lme4)
library(xtable)

# effect of sleep stage and EEG power for IED rate
l1 <- lmer(IEDsum ~ hyplabel + Zdelta + Ztheta + Zalpha + Zbeta + (1 | part) + (1 | patient), data)
plot_model(l1)

# write model coefficients to table for LaTeX
temp = summary(l1)
s1            <- temp$coefficients
rownames(s1)  <- c("Intercept", "Sleepstage (linear)", "Sleepstage (quadratic)", "Sleepstage (Cubic)", "Sleepstage (^4)", "Delta power", "Theta power", "Alpha power", "Beta power")
s1            <- data.frame(Predictor = row.names(s1), s1);
rownames(s1)  <- NULL
colnames(s1)  <- c("Predictor","Estimate","SD","df","t","p")
fname <- "D:/Dropbox/Apps/Overleaf/Hspike/IEDsum_coefficients.tex"
print(xtable(s1, type = "latex", 
             digits = c(0, 0,3,3,1,3,3), 
             caption = "Effect of sleep stage on IED rate",
             label = 'tab:IEDcoef'),
             caption.placement = "top", 
             include.rownames = FALSE,
             file = fname)

# posthoc tests
temp = emmeans(l1, list(pairwise ~ hyplabel), adjust = "tukey")
ph1 <- as.data.frame(temp$`pairwise differences of hyplabel`)
colnames(ph1) <- c("Comparison","Estimate","SE","df","Z ratio","p")
ph1 <- ph1[, -4] # remove df since they are at inf

# write posthoc comparisons to table for LaTeX
fname <- "D:/Dropbox/Apps/Overleaf/Hspike/IEDsum_coefficients_posthoc.tex"
print(xtable(ph1, type = "latex", 
             digits = c(0,0,3,3,3,3), 
             caption = "Posthoc Tukey comparison of sleep stages on IED rate",
             label = 'tab:IEDposthoc'),
      caption.placement = "top",
      include.rownames = FALSE,
      file = fname)




























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



