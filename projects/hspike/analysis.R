#install.packages("circular")
#install.packages("ggplot2")
#install.packages("units")
#install.packages("reshape2")
#install.packages("circlize")
#install.packages("ggthemes")
#install.packages("lemon")
#install.packages("egg")
#install.packages("readxl")
#install.packages('Rcpp')
#install.packages("equatiomatic")
#install.packages("kableExtra")
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#setTimeLimit(100000); setSessionTimeLimit(10000)
#devtools::install_github("strengejacke/sjPlot")
#install.packages("CircStats")

# ggpval

install.packages("spiralize")
library(spiralize)
library(ggplot2)
#library("cowplot")
#library("gridExtra")
library(plyr)
library(reshape)
library(RColorBrewer)
#library("ggthemes")
#library(lemon)
library("sjPlot")
library(emmeans)
library(gridExtra)
library(gtable)
library(grid)
library(egg)
library(ggpubr)
require(dplyr)
library(xtable)
library("readxl")
library(Rcpp)
library(equatiomatic)
# library(CircStats)
library(circular)
library(kableExtra)

#####################
## Support function #
#####################

# replace subsequent fields with NAN for visualization purposes
cleanf <- function(x){
  oldx <- c(FALSE, x[-1]==x[-length(x)])  # is the value equal to the previous?
  res <- x
  res[oldx] <- NA
  res}

###############################
# Latex table: Clinical table #
###############################

data_clinical <- read.csv("D:/Dropbox/Apps/Overleaf/Hspike/tables/clinical.csv")
data_clinical$ID <- NULL
data_clinical$Label <- NULL
colnames(data_clinical) <- c("Patient","Sex","Age","Onset","Type","SOZ","MRI","PET","SPECT","Medication","Implantation","Pre-implantation surgery" )
kbl(data_clinical, "latex", booktabs = T, label = "clinical",
    caption = "Clinical summary.
    SOZ: Seizure Onset Zone,
    MRI: Indications from Magnetic Resonance Imaging ,
    AED: Antiepileptic drugs,
    IEDs: Interictal Epileptiform discharges,
    FSWLA:  Focal seizures without loss of awareness,
    FSLA: Focal seizures with loss of awareness,
    FTBS: Focal to bilateral seizures, 
    PNH: Periventricular Nodular Heterotopia. 
    PMG: Polymicrogyria, 
    SNH: Subcortical Nodular Heterotopia,
    CBZ: Carbamazepine,
    LCS: Lacosamide,
    LTG: Lamotrigine, 
    ZNG: Zonisamide,
    ESL: Eslicarabazepine, 
    VPA: Valproic Acid,
    OXC: Oxicarbazepine,
    PER: Perampanel, 
    TPM: Topiramate",
    )%>%
  kable_styling(latex_options = c("scale_down"))%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  column_spec(1, width = "3em")  %>%
  column_spec(2, width = "1em")  %>%
  column_spec(3, width = "1em")  %>%
  column_spec(4, width = "2em")  %>%
  column_spec(5, width = "5em")  %>%
  column_spec(6, width = "5em")  %>%
  column_spec(7, width = "5em")  %>%
  column_spec(8, width = "5em")  %>%
  column_spec(9, width = "5em")  %>%
  column_spec(10, width = "5em")  %>%
  column_spec(11, width = "5em")  %>%
  collapse_rows(columns = 1) %>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/clinical.tex")

###############################################
# Latex table: electrode anatomical locations #
###############################################
options(knitr.kable.NA = '')

macro <- read.csv("D:/Dropbox/Apps/Overleaf/Hspike/tables/macro_anatomical.csv")
macro$ID <- NULL
macro$Label <- NULL
clean.cols <- c("Patient")
macro[clean.cols] <- lapply(macro[clean.cols], cleanf)

micro <- read.csv("D:/Dropbox/Apps/Overleaf/Hspike/tables/micro_anatomical.csv")
micro$ID <- NULL
micro$Label <- NULL
clean.cols <- c("Patient")
micro[clean.cols] <- lapply(micro[clean.cols], cleanf)

locations <- bind_rows(macro,micro)

kbl(locations, "latex", booktabs = T, linesep = "", label = 'anatomical',
    caption = "Anatomical locations of macro and micro electrodes")%>%
    kable_styling(latex_options = c("HOLD_position"))%>%
    pack_rows("Macro contacts", 1, 14) %>% # latex_gap_space = "2em"
    pack_rows("Micro electrodes", 15, 27) %>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/anatomical.tex")


#########################
# Detection performance #
#########################
options(knitr.kable.NA = '')
data  <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/performance.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data$Patient[9] = "\\textit{Mean}"
data$Patient[10] = "\\textit{Std.}"
data[,1]  = round(data[,1],digits=1)
data[,2]  = round(data[,2],digits=0)
data[,3]  = round(data[,3],digits=0)
data[,4]  = round(data[,4],digits=1)
data[,5]  = round(data[,5],digits=1)
data[,6]  = round(data[,6],digits=0)
data[,7]  = round(data[,7],digits=0)
data[,8]  = round(data[,8],digits=1)
data[,9]  = round(data[,9],digits=1)
data[,10] = round(data[,10],digits=0)
data[,11] = round(data[,11],digits=1)

kbl(data, "latex", booktabs = T, linesep = "", label = 'performance', escape = FALSE,
    col.names = c("Patient","24hrs", "24hrs","Hit (\\%)","FA (\\%)", "Total","24hrs","Hit (\\%)","FA (\\%)", "Total", "Total hrs."),
    caption = "Template detection performance")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  kable_styling(latex_options = c("scale_down"))%>%
  row_spec(8, hline_after = TRUE)%>%
  add_header_above(c(" " = 1, "Visual" = 1, "All templates" = 4, "Selected templates" = 4)) %>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/performance.tex")

########################################
# Latex table: Electrode locations MNI #
########################################

data  <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/MNI_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data  <- data[!data$color == 0, ] # only make table of used contacts
data  <- data[, c('patient','electrode','contact','X','Y','Z')]
colnames(data) = c('Patient','Electrode','Contact', 'X', 'Y', 'Z')
rownames(data) <- NULL
clean.cols <- c("Patient","Electrode")
data[clean.cols] <- lapply(data[clean.cols], cleanf)

kbl(data, "latex", booktabs = T, linesep = "", label = 'MNI',
    caption = "Anatomical locations of micro electrodes") %>%
    kable_styling(font_size = 6) %>%
    kable_styling(latex_options = c("HOLD_position"))%>%
  row_spec(5,  extra_latex_after = "\\cline{1-6}") %>%
  row_spec(10, extra_latex_after = "\\cline{2-6}") %>%
  row_spec(14, extra_latex_after = "\\cline{1-6}") %>%
  row_spec(19, extra_latex_after = "\\cline{1-6}") %>%
  row_spec(24, extra_latex_after = "\\cline{2-6}") %>%
  row_spec(29, extra_latex_after = "\\cline{1-6}") %>%
  row_spec(34, extra_latex_after = "\\cline{2-6}") %>%
  row_spec(39, extra_latex_after = "\\cline{1-6}") %>%
  row_spec(44, extra_latex_after = "\\cline{2-6}") %>%
  row_spec(49, extra_latex_after = "\\cline{1-6}") %>%
  row_spec(54, extra_latex_after = "\\cline{1-6}") %>%
  row_spec(59, extra_latex_after = "\\cline{2-6}") %>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/MNI.tex")

###################################
# Latex table: Time in sleepstage #
###################################

# normalize by time spend in sleep stages
data <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/hypnogram_duration.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data <- data[, c("patient", "part", "PHASE_3", "PHASE_2", "PHASE_1", "AWAKE", "REM")] 
colnames(data) = c("Patient", "Night","S3", "S2", "S1", "Wake", "REM")
clean.cols <- c("Patient")
data[clean.cols] <- lapply(data[clean.cols], cleanf)

kbl(data, "latex", booktabs = T, linesep = "", label = 'stageduration',
    caption = "Time spend in sleep stages (hrs.)", digits=2) %>%
    kable_styling(latex_options = c("HOLD_position"))%>%
    add_header_above(c(" " = 2, "Sleep stage" = 5)) %>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stageduration.tex")

#######################
# LFP power circadian #
#######################

# prepare data
data_power          <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/power_table_long.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_power$hyplabel <- factor(data_power$hyplabel, ordered = TRUE, levels = c("NO_SCORE", "REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))
# data_power$band     <- factor(data_power$band, ordered = TRUE, levels = c("delta", "theta", "alpha", "beta", "delta_div_alpha"))
data_power$band     <- factor(data_power$band, ordered = TRUE, levels = c("Delta1", "Delta2"))
data_power$Patient  <- factor(data_power$patient, levels = c(8:1))
data_power$part     <- factor(data_power$part)

# bin for polar representation
data_power$bin <- as.integer(cut(data_power$minute, seq(0, 24*60, by = 60)))
  
# duplicate midnight to connect in figure
temp <- subset(data_power, bin == 24)
temp$bin <- 0
data_power <- bind_rows(data_power, temp)

data_binned <- setNames(aggregate(data_power$power, c(list(data_power$Patient), list(data_power$bin), list(data_power$band)), mean), c("Patient", "bin", "band", "power"))
data_binned <- as.data.frame(data_binned %>% group_by(Patient, band) %>% mutate(Npower = (power-min(power))/max(power-min(power)))) 

# data_binned$rad <- data_binned$bin / 24 * pi
# c <- circular(control, units = "degrees", template = "geographics") 
# data_binned$rad
# d <- density.circular(data_binned$bin[data_binned$patient == 1], bw = 50)

###
# data_binned <- as.data.frame(data_binned %>% group_by(Patient, band) %>% mutate(Npower = (mean(power)-power)/sd(power))) 

####
# getCentroid <- function(x, width = 1) {
#   A  <- x * width                  # area of each bar
#   xc <- seq(width/2, length(x), 1) # x coordinates of center of bars
#   yc <- x/2                        # y coordinatey
#   
#   cx <- sum(xc * A) / sum(A)
#   cy <- sum(yc * A) / sum(A)
#   return(list(x = cx, y = cy))
# }
# points(getCentroid(x), col = 'red', pch = 19)
####

# plot
data_binned$title1 = "Delta1 (0.1-2.5 Hz)"
data_binned$title2 = "Delta2 (2.5-4 Hz)"
# data_binned$title2 = "Theta (5-7Hz)"
# data_binned$title3 = "Alpha (8-14Hz)"
# data_binned$title4 = "Delta (1-4Hz) / Alpha (8-14Hz)"

polarpowerplots <- list()

polarpowerplots[[1]] <- 
  ggplot(data=data_binned[data_binned$band == "delta", ], aes(x = bin, y = Npower, fill=Patient, col=Patient)) +
  scale_fill_brewer(palette = "Set2", direction = -1) + scale_color_brewer(palette = "Set2", direction = -1) +
  geom_vline(xintercept = seq(0, 24, by = 3), colour = "grey90") + 
  geom_hline(yintercept = seq(1, 8, by = 1), colour = "grey90") + 
  geom_ribbon(aes(ymin = (9-as.numeric(Patient)), ymax=Npower*2+(9-as.numeric(Patient))), alpha=0.8, colour = NA) +
  theme_article() +
  theme(
    panel.border = element_blank(),
    legend.text  = element_blank(),
    axis.ticks   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  scale_x_continuous(breaks=seq(0, 21, by = 3), 
                     labels = c("0" = "00:00", "3" = "", "6" = "06:00", "9" = "", "12" = "12:00", "15" = "", "18" = "18:00", "21" = "")) +
  coord_polar(theta = "x", start = 0, clip="off") +
  ylim(-2, 10) +
  facet_wrap(~title1)

polardelta1power <- 
  ggplot(data=data_binned[data_binned$band == "Delta1", ], aes(x = bin, y = Npower, fill=Patient, col=Patient)) +
  scale_fill_brewer(palette = "Set2", direction = -1) + scale_color_brewer(palette = "Set2", direction = -1) +
  geom_vline(xintercept = seq(0, 24, by = 3), colour = "grey90") + 
  geom_hline(yintercept = seq(1, 8, by = 1), colour = "grey90") + 
  geom_ribbon(aes(ymin = (9-as.numeric(Patient)), ymax=Npower*1.5+(9-as.numeric(Patient))), alpha=1, colour = NA) +
  theme_article() +
  theme(
    panel.border = element_blank(),
    legend.text  = element_blank(),
    axis.ticks   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  scale_x_continuous(breaks=seq(0, 21, by = 3), 
                     labels = c("0" = "00:00", "3" = "", "6" = "06:00", "9" = "", "12" = "12:00", "15" = "", "18" = "18:00", "21" = "")) +
  coord_polar(theta = "x", start = 0, clip="off") +
  ylim(-2, 10)

polardelta2power <- 
  ggplot(data=data_binned[data_binned$band == "Delta2", ], aes(x = bin, y = Npower, fill=Patient, col=Patient)) +
  scale_fill_brewer(palette = "Set2", direction = -1) + scale_color_brewer(palette = "Set2", direction = -1) +
  geom_vline(xintercept = seq(0, 24, by = 3), colour = "grey90") + 
  geom_hline(yintercept = seq(0, 8, by = 1), colour = "grey90") + 
  geom_ribbon(aes(ymin = (9-as.numeric(Patient)), ymax=Npower*1.5+(9-as.numeric(Patient))), alpha=1, colour = NA) +
  theme_article() +
  theme(
    panel.border = element_blank(),
    legend.text  = element_blank(),
    axis.ticks   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  scale_x_continuous(breaks=seq(0, 21, by = 3), 
                     labels = c("0" = "00:00", "3" = "", "6" = "06:00", "9" = "", "12" = "12:00", "15" = "", "18" = "18:00", "21" = "")) +
  coord_polar(theta = "x", start = 0, clip="off") +
  ylim(-1, 10)


polarpowerplots[[2]] <- 
  ggplot(data=data_binned[data_binned$band == "theta", ], aes(x = bin, y = Npower, fill=Patient, col=Patient)) +
  scale_fill_brewer(palette = "Set2", direction = -1) + scale_color_brewer(palette = "Set2", direction = -1) +
  geom_vline(xintercept = seq(0, 24, by = 3), colour = "grey90") + 
  geom_hline(yintercept = seq(0, 8, by = 1), colour = "grey90") + 
  geom_ribbon(aes(ymin = (9-as.numeric(Patient)), ymax=Npower*2+(9-as.numeric(Patient))), alpha=0.8, colour = NA) +
  theme_article() +
  theme(
    panel.border = element_blank(),
    legend.text  = element_blank(),
    axis.ticks   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  scale_x_continuous(breaks=seq(0, 21, by = 3), 
                     labels = c("0" = "00:00", "3" = "", "6" = "06:00", "9" = "", "12" = "12:00", "15" = "", "18" = "18:00", "21" = "")) +
  coord_polar(theta = "x", start = 0, clip="off") +
  ylim(0, 10) +
  facet_wrap(~title2)

polarpowerplots[[3]] <- 
  ggplot(data=data_binned[data_binned$band == "alpha", ], aes(x = bin, y = Npower, fill=Patient, col=Patient)) +
  scale_fill_brewer(palette = "Set2", direction = -1) + scale_color_brewer(palette = "Set2", direction = -1) +
  geom_vline(xintercept = seq(0, 24, by = 3), colour = "grey90") + 
  geom_hline(yintercept = seq(0, 8, by = 1), colour = "grey90") + 
  geom_ribbon(aes(ymin = (9-as.numeric(Patient)), ymax=Npower*2+(9-as.numeric(Patient))), alpha=0.8, colour = NA) +
  theme_article() +
  theme(
    panel.border = element_blank(),
    legend.text  = element_blank(),
    axis.ticks   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  scale_x_continuous(breaks=seq(0, 21, by = 3), 
                     labels = c("0" = "00:00", "3" = "", "6" = "06:00", "9" = "", "12" = "12:00", "15" = "", "18" = "18:00", "21" = "")) +
  coord_polar(theta = "x", start = 0, clip="off") +
  ylim(0, 10) +
  facet_wrap(~title3)

polarpowerplots[[4]] <-
  ggplot(data=data_binned[data_binned$band == "delta_div_alpha", ], aes(x = bin, y = Npower, fill=Patient, col=Patient)) +
  scale_fill_brewer(palette = "Set2", direction = -1) + scale_color_brewer(palette = "Set2", direction = -1) +
  geom_vline(xintercept = seq(0, 24, by = 3), colour = "grey90") + 
  geom_hline(yintercept = seq(0, 8, by = 1), colour = "grey90") + 
  geom_ribbon(aes(ymin = (9-as.numeric(Patient)), ymax=Npower*2+(9-as.numeric(Patient))), alpha=0.8, colour = NA) +
  theme_article() +
  theme(
    panel.border = element_blank(),
    legend.text  = element_blank(),
    axis.ticks   = element_blank(),
    axis.text.y  = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) +
  scale_x_continuous(breaks=seq(0, 21, by = 3), 
                     labels = c("0" = "00:00", "3" = "", "6" = "06:00", "9" = "", "12" = "12:00", "15" = "", "18" = "18:00", "21" = "")) +
  coord_polar(theta = "x", start = 0, clip="off") +
  ylim(0, 10) +
  facet_wrap(~title4)

powerplots <- ggarrange(plotlist=polarpowerplots[c(1,2,3)], widths = c(1,1,1,1), heights = c(1,1,1,1), nrow = 1, 
          labels = c("D","E","F"), vjust = 18, hjust = -1, 
          legend = "right", common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold"))

ggarrange(plotlist=polarpowerplots[c(1,2,3)], widths = c(1,1,1,1), heights = c(1,1,1,1), nrow = 1, 
                        labels = c("A","B","C"), vjust = 18, hjust = -1, 
                        legend = "right", common.legend = TRUE, 
                        font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/polar_power_band.pdf")


######################
# IED rate circadian #
######################

data_IED                <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/IED_table_PSG.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_IED[data_IED$hyplabel == "PRE_SLEEP", ] = "AWAKE"
data_IED[data_IED$hyplabel == "POST_SLEEP", ] = "AWAKE"
data_IED$hyplabel       <- factor(data_IED$hyplabel, ordered = TRUE, levels = c("REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))

data_IED$Patient        <- factor(data_IED$patient, levels = c(8:1)) # same order in plot
data_IED$part           <- factor(data_IED$part)
data_IED$marker         <- factor(data_IED$marker)
data_IED$hour           <- data_IED$minute / (60)
data_IED$rad            <- data_IED$theta

# extract distribution statistics
dist_IED <- data.frame()  
d <- list()
stat_IED <- list()
for (ipatient in 1:8) {
  d <- density.circular(data_IED$rad[data_IED$patient == ipatient], bw = 100)
  temp = list()
  temp$x = as.numeric(d$x)
  temp$y = d$y
  le <- lengths(temp)
  temp$patient <- rep(ipatient,le[1])
  tempdf <- as.data.frame(temp)
  colnames(tempdf)           = c('rad','density','Patient')
  dist_IED <- bind_rows(dist_IED, tempdf)
  temp                       = list()
  
  # Rayleigh Test of Uniformity: General Unimodal Alternative
  temp                       <- rayleigh.test(data_IED$rad[data_IED$patient == ipatient])
  stat_IED$Patient[ipatient]    = ipatient
  stat_IED$p[ipatient]          = temp$p.value
  stat_IED$Rayleigh_stat[ipatient] = temp$statistic
  
  # Rayleigh Test of Uniformity: General Unimodal Alternative
  stat_IED$median[ipatient] <- as.numeric(median(circular(data_IED$rad[data_IED$patient == ipatient]))) / (pi * 2) * 24
  if (stat_IED$median[ipatient] < 0) {
    stat_IED$median[ipatient] = stat_IED$median[ipatient] + 24
  }
}

# format data for plotting
dist_IED$Khour = dist_IED$rad / (pi * 2) * 24
dist_IED$Patient = factor(dist_IED$Patient, levels = c(1:8)) # reversed order in plot
stat_IED <- as.data.frame(stat_IED)
stat_IED$Patient <- factor(stat_IED$Patient, levels = c(1:8))

# add some offset in degrees so it shows up from behind the rest
stat_IED$medianplot <- stat_IED$median
stat_IED$medianplot[1] = stat_IED$median[1] - 0.18 
stat_IED$medianplot[6] = stat_IED$median[6] + 0.2
stat_IED$medianplot[stat_IED$p >= 0.05] = NA # all are significant though

IEDpolarplot <-
  ggplot(data=data_IED, aes(x=hour, y = as.numeric(Patient), fill=Patient)) + 
  geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
  geom_hline(yintercept = seq(1, 8, by = 1),  colour = "grey90") +  
  geom_ribbon(data=dist_IED, alpha = 1, colour = NA, aes(x = Khour, 
                                                   ymin = (9-as.numeric(Patient)),
                                                   ymax = density*3 + (9-as.numeric(Patient)), 
                                                   col = Patient), show.legend = TRUE) +
  
  geom_segment(data=stat_IED, aes(x=medianplot, y=9.5, xend=medianplot, yend=10, col=Patient), 
               arrow = arrow(length = unit(0.25, "cm"), type="closed"), size = 0.5, show.legend = FALSE) + 
  
  coord_polar(theta = "x", start = 0, direction = 1, clip = 'off') + 
  scale_x_continuous(breaks = seq(0, 21, by = 3), 
                     labels = c("0" = "00:00", "3" = "", "6" = "06:00", "9" = "", "12" = "12:00", "15" = "", "18" = "18:00", "21" = "")) +
  scale_fill_brewer(palette = "Set2", direction=-1) + scale_color_brewer(palette = "Set2", direction=-1) +
  theme_article() +
  theme(panel.border = element_blank(),
        #legend.key   = element_blank(),
        axis.ticks   = element_blank(),
        axis.text.y  = element_blank(),
        #axis.text.x  = element_blank(),
        panel.grid   = element_blank(),
        axis.title.x = element_blank(),
        #legend.text  = element_blank(),
        axis.title.y = element_blank()) + 
  ylim(-2, 10)

# LaTeX table
library(stringr)
stat_IED$time = paste(str_pad( floor(stat_IED$median), 2, pad = "0"), ":", str_pad(floor((stat_IED$median- floor(stat_IED$median)) * 60), 2, pad = "0"), sep = "")
stat_IED <- stat_IED[, c("Patient", "time", "Rayleigh_stat", "p")]
stat_IED[,4] = ifelse(stat_IED[,4] > .05, paste(round(stat_IED[,4],digits=2),sep=""), ifelse(stat_IED[,4] < .0001, "<.0001\\textsuperscript{***}", ifelse(stat_IED[,4] < .001,"<.001\\textsuperscript{**}", ifelse(stat_IED[,4] < .01, "<.01\\textsuperscript{*}", "<.05"))))

kbl(stat_IED, "latex", booktabs = T, linesep = "", label = 'circstat_IED',
    col.names = c("Patient","Median angle (HH:mm)","Rayleigh", "\\textit{p}"),
    escape = FALSE, digits = 2,
    caption = "Circular statistics of circadian IED rate")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/circstat_IED.tex")

######################
# Seizures circadian #
###################### 

# load data: Patients x Units x time window
data_seizures           <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/seizuredata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")

# prepare data
data_seizures$Patient   <- factor(data_seizures$patient, levels = c(1:8))
data_seizures$hour      <- data_seizures$minute / 60
data_seizures$rad       <- data_seizures$minute / 60 / 24 * pi * 2

# extract distribution statistics
dist_seizures <- data.frame()  
d <- list()
stat_seizures <- list()
for (ipatient in 1:8) {
  d <- density.circular(data_seizures$rad[data_seizures$patient == ipatient], bw = 50)
  temp = list()
  temp$x = as.numeric(d$x)
  temp$y = d$y
  le <- lengths(temp)
  temp$patient <- rep(ipatient,le[1])
  tempdf <- as.data.frame(temp)
  colnames(tempdf) = c('rad','density','Patient')
  dist_seizures <- bind_rows(dist_seizures, tempdf)
  
  # Rayleigh Test of Uniformity: General Unimodal Alternative
  temp = list()
  temp <- rayleigh.test(data_seizures$rad[data_seizures$patient == ipatient])
  stat_seizures$Patient[ipatient]    = ipatient
  stat_seizures$p[ipatient]          = temp$p.value
  stat_seizures$Rayleigh_stat[ipatient] = temp$statistic
  
  # Rayleigh Test of Uniformity: General Unimodal Alternative
  stat_seizures$median[ipatient] <- as.numeric(median(circular(data_seizures$rad[data_seizures$patient == ipatient]))) / (pi * 2) * 24
  if (stat_seizures$median[ipatient] < 0) {
    stat_seizures$median[ipatient] = stat_seizures$median[ipatient] + 24
  }
}

# format data for plotting
dist_seizures$Khour = dist_seizures$rad / (pi * 2) * 24
dist_seizures$Patient = factor(dist_seizures$Patient, levels = c(8:1))
stat_seizures <- as.data.frame(stat_seizures)
stat_seizures$Patient = factor(stat_seizures$Patient, levels = c(8:1))
stat_seizures$significant = stat_seizures$p < 0.05
stat_seizures$median_sel = stat_seizures$median
stat_seizures$median_sel[stat_seizures$p >= 0.05] = NA # all are significant though

# tiny adjustment to make axes line out properly and arrows not overlap
data_seizures$hourplot <- data_seizures$hour
data_seizures$hourplot[which(data_seizures$hour==max(data_seizures$hour))] = 24

# add jitter function for points
jitter <- position_jitter(width = 0, height = 0.4)

# plot
Seizurepolarplot  <-
  ggplot(data=data_seizures, aes(x=hourplot, y = as.numeric(Patient), fill=Patient)) + 
  geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
  geom_hline(yintercept = seq(1, 8, by = 1), colour = "grey90") +  
  geom_ribbon(data=dist_seizures, alpha = 1, colour = NA, aes(x=Khour, 
                                                   ymin = (9-as.numeric(Patient)),
                                                   ymax = density*1.3 + (9-as.numeric(Patient)), 
                                                   col = Patient), show.legend = FALSE) +
  
  geom_segment(data=stat_seizures, aes(x=median_sel, y=9.5, xend=median_sel, yend=10, col=Patient), 
               arrow = arrow(length = unit(0.25, "cm"), type="closed"), size = 1, show.legend = FALSE) + 
  
  geom_point(colour="black", pch=21, size=1, position = jitter) +
  coord_polar(theta = "x", start = 0, direction = 1, clip = 'off') + 
  scale_x_continuous(breaks = seq(0, 21, by = 3),
                     labels = c("0" = "00:00", "3" = "", "6" = "06:00", "9" = "", "12" = "12:00", "15" = "", "18" = "18:00", "21" = "")) +
  scale_fill_brewer(palette = "Set2", direction = -1) + scale_color_brewer(palette = "Set2", direction = -1) +
  theme_article() +
  theme(panel.border = element_blank(),
        #legend.key   = element_blank(),
        axis.ticks   = element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x  = element_blank(),
        panel.grid   = element_blank(),
        axis.title.x = element_blank(),
        #legend.text  = element_blank(),
        axis.title.y = element_blank()) + 
  ylim(-2, 10)

# save combined to pdf
ggarrange(IEDpolarplot, Seizurepolarplot, 
          labels = c("A","B"), 
          vjust = 15, hjust = -1, 
          legend = "right", 
          common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/polar_density_seizures.pdf") 


# save combined to pdf
ggarrange(IEDpolarplot, Seizurepolarplot, polardelta1power,
          labels = c("A","B","C"), 
          vjust = 3, hjust = -1, 
          legend = "right", 
          common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/polar.pdf") 


# LaTeX table
library(stringr)
stat_seizures$time = paste(str_pad( floor(stat_seizures$median), 2, pad = "0"), ":", str_pad(floor((stat_seizures$median- floor(stat_seizures$median)) * 60), 2, pad = "0"), sep = "")
stat_seizures <- stat_seizures[, c("Patient", "time", "Rayleigh_stat", "p")]
stat_seizures[,4] = ifelse(stat_seizures[,4] > .05, paste(round(stat_seizures[,4],digits=2),sep=""), ifelse(stat_seizures[,4] < .0001, "<.0001\\textsuperscript{***}", ifelse(stat_seizures[,4] < .001,"<.001\\textsuperscript{**}", ifelse(stat_seizures[,4] < .01, "<.01\\textsuperscript{*}", "<.05"))))

kbl(stat_seizures, "latex", booktabs = T, linesep = "", label = 'circstat_seizures',
    col.names = c("Patient","Median angle (HH:mm)","Rayleigh", "\\textit{p}"),
    escape = FALSE, digits = 2,
    caption = "Circular statistics of circadian seizure occurance")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/circstat_seizures.tex")

# combined table 
temp <- stat_seizures
colnames(temp) <- c("Patient","time2","Rayleigh_stat2","p2")
combined = merge(stat_IED,temp)

kbl(combined, "latex", booktabs = T, linesep = "", label = 'circstats',
    col.names = c("Patient","Time","Rayleigh", "\\textit{p}","Time","Rayleigh", "\\textit{p}"),
    escape = FALSE, digits = 2,
    caption = "Circular statistics of circadian epileptic activity")%>%
    add_header_above(c(" ", "Interictal activity" = 3, "Seizures" = 3)) %>%
    kable_styling(latex_options = c("HOLD_position"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/circstats.tex")

#############################
# LFP power per sleep stage #
#############################

# prepare data
data_pow  <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/power_table_long.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_pow  <- data_pow[!data_pow$part > 3, ] # hypnogram is only scored on first three nights
data_pow  <- data_pow[!data_pow$hyplabel == "NO_SCORE", ]
data_pow$hyplabel[data_pow$hyplabel == "PHASE_1"]   = "S1"
data_pow$hyplabel[data_pow$hyplabel == "PHASE_2"]   = "S2"
data_pow$hyplabel[data_pow$hyplabel == "PHASE_3"]   = "S3"
data_pow$hyplabel[data_pow$hyplabel == "AWAKE"]     = "Wake"
data_pow$hyplabel[data_pow$hyplabel == "PRESLEEP"]  = "Pre"
data_pow$hyplabel[data_pow$hyplabel == "POSTSLEEP"] = "Post"
data_pow$hyplabel <- factor(data_pow$hyplabel, levels = c("REM", "Post", "Pre", "Wake", "S1", "S2", "S3"))
# data_pow$band     <- factor(data_pow$band, ordered = TRUE, levels = c("delta", "theta", "alpha"))
data_pow$band     <- factor(data_pow$band, ordered = TRUE, levels = c("Delta1", "Delta2"))
data_pow$patient  <- factor(data_pow$patient, levels = c(8:1))
data_pow$part     <- factor(data_pow$part)

# add clinical data (not yet used here)
data_clinical$patient <- data_clinical$Patient
data_pow <- merge(data_pow, data_clinical)

# scale data over windows
m               <- setNames(aggregate(data_pow$power, by = list(data_pow$patient, data_pow$band, data_pow$hyplabel), median), c("patient", "band", "hyplabel", "power_median"))
data_pow        <- merge(data_pow, m)
m               <- setNames(aggregate(data_pow$power, by = list(data_pow$patient, data_pow$band), mean), c("patient", "band", "power_avg"))
s               <- setNames(aggregate(data_pow$power, by = list(data_pow$patient, data_pow$band), sd),   c("patient", "band", "power_sd"))
data_pow        <- merge(data_pow, m)
data_pow        <- merge(data_pow, s)
data_pow$Zpower <- (data_pow$power-data_pow$power_avg)/data_pow$power_sd
m               <- setNames(aggregate(data_pow$Zpower, by = list(data_pow$patient, data_pow$part, data_pow$band, data_pow$hyplabel), mean), c("patient", "part", "band", "hyplabel", "power_avg"))
m$patient       <- factor(m$patient,  levels = c(8:1))
m$hyplabel      <- factor(m$hyplabel, levels = c("REM", "Post", "Pre", "Wake", "S1", "S2", "S3")) # data$hyplabel <- factor(data$hyplabel, ordered = TRUE, levels = c("REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"))
m$part          <- factor(m$part)

# plot
data_pow$title1 = "Normalized power Delta1 (0.1-2.5 Hz)"
data_pow$title2 = "Normalized power Delta2 (2.5-4 Hz)"
# data_pow$title2 = "Normalized theta power (5-7Hz)"
# data_pow$title3 = "Normalized alpha power (8-14Hz)"

plots_pow <- list()
plots_delta1 <- 
  ggplot(data=data_pow[data_pow$band == "Delta1",], aes(y = hyplabel, x = Zpower)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_point(data = m[m$band == "Delta1",], aes(group = interaction(patient, part), x = power_avg, y = hyplabel, col = patient), alpha = 0.5,
             position=position_dodge(width=0.5)) +
  guides(colour = "none") +
  coord_cartesian(xlim = c(-1.2, 2)) +
  theme_article() +
  ylab(NULL) + xlab(NULL) +
  facet_wrap(~title1)

plots_delta2 <- 
  ggplot(data=data_pow[data_pow$band == "Delta2",], aes(y = hyplabel, x = Zpower)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_point(data = m[m$band == "Delta2",], aes(group = interaction(patient, part), x = power_avg, y = hyplabel, col = patient), alpha = 0.5,
             position=position_dodge(width=0.5)) +
  guides(colour = "none") +
  coord_cartesian(xlim = c(-1.2, 2)) +
  theme_article() +
  ylab(NULL) + xlab(NULL) +
  facet_wrap(~title2)

plots_pow[[2]] <- 
  ggplot(data=data_pow[data_pow$band == "theta",], aes(y = hyplabel, x = Zpower)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_point(data = m[m$band == "theta",], aes(group = interaction(patient, part), x = power_avg, y = hyplabel, col = patient), alpha = 0.5, 
             position=position_dodge(width=0.5)) +
  guides(colour = "none") +
  coord_cartesian(xlim = c(-1.2, 2)) +
  theme_article() +
  ylab(NULL) + xlab(NULL) +
  facet_wrap(~title2)

plots_pow[[3]] <- 
  ggplot(data=data_pow[data_pow$band == "alpha",], aes(y = hyplabel, x = Zpower)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_point(data = m[m$band == "alpha",], aes(group = interaction(patient, part), x = power_avg, y = hyplabel, col = patient), alpha = 0.5, 
             position=position_dodge(width=0.5)) +
  guides(colour = "none") +
  coord_cartesian(xlim = c(-1.2, 2)) +
  theme_article() +
  ylab(NULL) + xlab(NULL) +
  facet_wrap(~title3)

##########################
# IED rate & sleep stage #
##########################

# prepare data
data_IED <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/IED_table_PSG.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_IED <- data_IED[!data_IED$hyplabel == "NO_SCORE", ]
data_IED$hyplabel[data_IED$hyplabel == "PHASE_1"]   = "S1"
data_IED$hyplabel[data_IED$hyplabel == "PHASE_2"]   = "S2"
data_IED$hyplabel[data_IED$hyplabel == "PHASE_3"]   = "S3"
data_IED$hyplabel[data_IED$hyplabel == "AWAKE"]     = "Wake"
data_IED$hyplabel[data_IED$hyplabel == "PRESLEEP"]  = "Pre"
data_IED$hyplabel[data_IED$hyplabel == "POSTSLEEP"] = "Post"
data_IED$hyplabel <- factor(data_IED$hyplabel, levels = c("REM", "Post", "Pre", "Wake", "S1", "S2", "S3"))
data_IED$patient  <- factor(data_IED$patient, levels = c(8:1))
data_IED$part     <- factor(data_IED$part)
data_IED$marker   <- factor(data_IED$marker)
data_IED$hour     <- data_IED$minute / (60)

# normalize by time spend in sleep stages
data_duration <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/hypnogram_duration.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_duration <- data_duration %>% rename(
  'S1' = 'PHASE_1',
  'S2' = 'PHASE_2',
  'S3' = 'PHASE_3',
  'Wake' = 'AWAKE',
  'Pre' = 'PRESLEEP',
  'Post' = 'POSTSLEEP')

data_duration <- melt(data_duration, id = c('patient','part'))
colnames(data_duration) = c('patient','part','hyplabel','duration')

# count IEDs per sleepstage
s <- na.omit(data_IED %>% dplyr::count(patient, part, hyplabel)) # library conflict of count 

# normalize by time spend in sleep stages
s <- merge(s, data_duration)
s$IEDrate = s$n / s$duration / 60 # count is per 10 seconds

# normalize rate by REM rate 
i <- s[s$hyplabel == 'REM',]
i <- i[, c('patient','part','IEDrate')]
colnames(i) = c('patient','part','REMrate')
s <- merge(s, i)
s$IEDrateNorm <- s$IEDrate / s$REMrate

# plot
s$title1 = "IED count"
s$title2 = "IED rate (count/min)"
s$title3 = "IED normalized to REM"
plot_IED_count <- ggplot(data=s, aes(y=hyplabel, x=n)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(patient, part), col = patient), position=position_dodge(width=0.3), alpha = 0.5) +
  guides(colour = "none") +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  coord_cartesian(xlim = c(1, 4000)) +
  facet_wrap(~title1)

plot_IED_rate <- ggplot(data=s, aes(y=hyplabel, x=IEDrate)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(patient, part), col = patient), position=position_dodge(width=0.3), alpha = 0.5) +
  guides(colour = "none") +
  ylab(NULL) + xlab(NULL) +
  coord_cartesian(xlim = c(0, 20)) +
  theme_article() + 
  facet_wrap(~title2)

plot_IED_norm <- ggplot(data=s, aes(y=hyplabel, x=IEDrateNorm)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(patient, part), col = patient), position=position_dodge(width=0.3), alpha = 0.5) + 
  ylab(NULL) + xlab(NULL) +
  theme(axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(0, 20)) +
  theme_article() +
  theme(legend.position ="bottom") +
  guides(colour = guide_legend(ncol = 1)) +
  guides(fill=guide_legend(title="Patient")) + 
  labs(color='Patient') + 
  facet_wrap(~title3) 

#################################
# IED rate vs. power STATISTICS #
#################################

# prepare data
data_pow_wide <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/power_table_wide.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_pow_wide <- data_pow_wide[!data_pow_wide$hyplabel == "NO_SCORE", ]
data_pow_wide$hyplabel[data_pow_wide$hyplabel == "PHASE_1"]   = "S1"
data_pow_wide$hyplabel[data_pow_wide$hyplabel == "PHASE_2"]   = "S2"
data_pow_wide$hyplabel[data_pow_wide$hyplabel == "PHASE_3"]   = "S3"
data_pow_wide$hyplabel[data_pow_wide$hyplabel == "AWAKE"]     = "Wake"
data_pow_wide$hyplabel[data_pow_wide$hyplabel == "PRESLEEP"]  = "Pre"
data_pow_wide$hyplabel[data_pow_wide$hyplabel == "POSTSLEEP"] = "Post"
data_pow_wide$hyplabel <- factor(data_pow_wide$hyplabel, levels = c("Pre","S3", "S2", "S1", "Wake", "REM", "Post"))
data_pow_wide$stage    <- data_pow_wide$hyplabel
data_pow_wide$patient  <- factor(data_pow_wide$patient, levels = c(8:1))
data_pow_wide$night    <- factor(data_pow_wide$part)

# add clinical data (not yet used here)
data_clinical$patient <- data_clinical$Patient
data_pow_wide <- merge(data_pow_wide, data_clinical)
data_pow_wide$duration <- data_pow_wide$Age - data_pow_wide$Onset


# normalize 
d1   <- setNames(aggregate(data_pow_wide$Delta1, by = c(list(data_pow_wide$patient)), mean), c("patient", "Mdelta1"))
d1sd <- setNames(aggregate(data_pow_wide$Delta1, by = c(list(data_pow_wide$patient)), sd),   c("patient", "SDdelta1"))
d2   <- setNames(aggregate(data_pow_wide$Delta2, by = c(list(data_pow_wide$patient)), mean), c("patient", "Mdelta2"))
d2sd <- setNames(aggregate(data_pow_wide$Delta2, by = c(list(data_pow_wide$patient)), sd),   c("patient", "SDdelta2"))

data_pow_wide <- merge(data_pow_wide,d1)
data_pow_wide <- merge(data_pow_wide,d2)
data_pow_wide <- merge(data_pow_wide,d1sd)
data_pow_wide <- merge(data_pow_wide,d2sd)

data_pow_wide$Zdelta1 <- (data_pow_wide$Delta1-data_pow_wide$Mdelta1)/data_pow_wide$SDdelta1
data_pow_wide$Zdelta2 <- (data_pow_wide$Delta2-data_pow_wide$Mdelta2)/data_pow_wide$SDdelta2

#   
# # normalize 
# y1   <- setNames(aggregate(data_pow_wide$alpha, by = c(list(data_pow_wide$patient)), mean), c("patient", "Malpha"))
# y2   <- setNames(aggregate(data_pow_wide$theta, by = c(list(data_pow_wide$patient)), mean), c("patient", "Mtheta"))
# y3   <- setNames(aggregate(data_pow_wide$beta,  by = c(list(data_pow_wide$patient)), mean), c("patient", "Mbeta"))
# y4   <- setNames(aggregate(data_pow_wide$delta, by = c(list(data_pow_wide$patient)), mean), c("patient", "Mdelta"))
# y1sd <- setNames(aggregate(data_pow_wide$alpha, by = c(list(data_pow_wide$patient)), sd),   c("patient", "SDalpha"))
# y2sd <- setNames(aggregate(data_pow_wide$theta, by = c(list(data_pow_wide$patient)), sd),   c("patient", "SDtheta"))
# y3sd <- setNames(aggregate(data_pow_wide$beta,  by = c(list(data_pow_wide$patient)), sd),   c("patient", "SDbeta"))
# y4sd <- setNames(aggregate(data_pow_wide$delta, by = c(list(data_pow_wide$patient)), sd),   c("patient", "SDdelta"))
# 
# data_pow_wide <- merge(data_pow_wide,y1)
# data_pow_wide <- merge(data_pow_wide,y2)
# data_pow_wide <- merge(data_pow_wide,y3)
# data_pow_wide <- merge(data_pow_wide,y4)
# data_pow_wide <- merge(data_pow_wide,y1sd)
# data_pow_wide <- merge(data_pow_wide,y2sd)
# data_pow_wide <- merge(data_pow_wide,y3sd)
# data_pow_wide <- merge(data_pow_wide,y4sd)
# 
# data_pow_wide$Zdelta <- (data_pow_wide$delta-data_pow_wide$Mdelta)/data_pow_wide$SDdelta
# data_pow_wide$Ztheta <- (data_pow_wide$theta-data_pow_wide$Mtheta)/data_pow_wide$SDtheta
# data_pow_wide$Zalpha <- (data_pow_wide$alpha-data_pow_wide$Malpha)/data_pow_wide$SDalpha
# data_pow_wide$Zbeta  <- (data_pow_wide$beta -data_pow_wide$Mbeta)/data_pow_wide$SDbeta

# to create p-values
detach(package:lmerTest)
library(lmerTest)
library(lme4)

####
# Delta1 explained by sleep stage
####

l1 <- lmer(Zdelta1 ~ stage + (1 | night) + (1 | patient), data_pow_wide, control = lmerControl(optimizer ='Nelder_Mead'))
# l1 <- lm(Zdelta ~ stage + duration, data_pow_wide, control = lmerControl(optimizer ='Nelder_Mead'))
summary(l1)
plot_model(l1)

# Coefficients
temp = summary(l1)
coefs <- as.data.frame(temp$coefficients)
coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post")

# Post-hoc tests
temp = emmeans(l1, list(pairwise ~ stage), adjust = "tukey")
ph1 <- as.data.frame(temp$`pairwise differences of stage`)
ph1 <- ph1[, -4] # remove df since they are at inf
ph1[,5] = ifelse(ph1[,5] > .05, paste(round(ph1[,5],digits=2),sep=""), ifelse(ph1[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(ph1[,5] < .001,"<.001\\textsuperscript{**}", ifelse(ph1[,5] < .01, "<.01\\textsuperscript{*}", "<.05"))))

# Concatenate in one LaTeX table
coefs           <- data.frame(Predictor = row.names(coefs), coefs);
rownames(coefs) <- NULL
colnames(ph1)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","z", "\\textit{p}")
colnames(coefs) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)", "df", "z","\\textit{p}")
stats_delta1   <- bind_rows(coefs,ph1)

options(knitr.kable.NA = '')

kbl(stats_delta1, "latex", booktabs = T, linesep = "", label = 'stats_delta1',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stages on Delta1 (0.1-2.5Hz) power")%>%
    pack_rows("Coefficients", 1, 7) %>%
    pack_rows("Post-hoc comparisons", 8, 28) %>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_delta1.tex")

####
# Delta2 explained by sleep stage
####
l1 <- lmer(Zdelta2 ~ stage + (1 | night) + (1 | patient), data_pow_wide, control = lmerControl(optimizer ='Nelder_Mead'))
# l1 <- lm(Zdelta ~ stage + duration, data_pow_wide, control = lmerControl(optimizer ='Nelder_Mead'))
summary(l1)
plot_model(l1)

# Coefficients
temp = summary(l1)
coefs <- as.data.frame(temp$coefficients)
coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post")

# Post-hoc tests
temp = emmeans(l1, list(pairwise ~ stage), adjust = "tukey")
ph1 <- as.data.frame(temp$`pairwise differences of stage`)
ph1 <- ph1[, -4] # remove df since they are at inf
ph1[,5] = ifelse(ph1[,5] > .05, paste(round(ph1[,5],digits=2),sep=""), ifelse(ph1[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(ph1[,5] < .001,"<.001\\textsuperscript{**}", ifelse(ph1[,5] < .01, "<.01\\textsuperscript{*}", "<.05"))))

# Concatenate in one LaTeX table
coefs           <- data.frame(Predictor = row.names(coefs), coefs);
rownames(coefs) <- NULL
colnames(ph1)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","z", "\\textit{p}")
colnames(coefs) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)", "df", "z","\\textit{p}")
stats_delta2   <- bind_rows(coefs,ph1)

options(knitr.kable.NA = '')

kbl(stats_delta2, "latex", booktabs = T, linesep = "", label = 'stats_delta2',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stages on Delta2 (2.5-4 Hz) power")%>%
  pack_rows("Coefficients", 1, 7) %>%
  pack_rows("Post-hoc comparisons", 8, 28) %>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_delta2.tex")

# # THETA explained by sleep stage
# l1 <- lmer(Ztheta ~ stage + (1 | night) + (1 | patient), data_pow_wide, control = lmerControl(optimizer ='Nelder_Mead'))
# summary(l1)
# plot_model(l1)
# 
# # THETA Coefficients to LaTeX table
# temp = summary(l1)
# coefs <- as.data.frame(temp$coefficients)
# coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
# rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post")
# 
# # THETA Post-hoc tests to LaTeX table
# temp = emmeans(l1, list(pairwise ~ stage), adjust = "tukey")
# ph1 <- as.data.frame(temp$`pairwise differences of stage`)
# ph1 <- ph1[, -4] # remove df since they are at inf
# ph1[,5] = ifelse(ph1[,5] > .05, paste(round(ph1[,5],digits=2),sep=""), ifelse(ph1[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(ph1[,5] < .001,"<.001\\textsuperscript{**}", ifelse(ph1[,5] < .01, "<.01\\textsuperscript{*}", "<.05"))))
# 
# # Concatenate in one LaTeX table
# coefs           <- data.frame(Predictor = row.names(coefs), coefs);
# rownames(coefs) <- NULL
# colnames(ph1)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","z", "\\textit{p}")
# colnames(coefs) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)", "df", "z","\\textit{p}")
# stats_theta     <- bind_rows(coefs,ph1)
# 
# kbl(stats_theta, "latex", booktabs = T, linesep = "", label = 'stats_theta',
#     escape = FALSE, digits = 2,
#     caption = "Effect of sleep stages on Theta power")%>%
#   pack_rows("Coefficients", 1, 7) %>%
#   pack_rows("Post-hoc comparisons", 8, 28) %>%
#   kable_styling(latex_options = c("HOLD_position"))%>%
#   save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_theta.tex")
# 
# # ALPHA explained by sleep stage
# l1 <- lmer(Zalpha ~ stage + (1 | night) + (1 | patient), data_pow_wide, control = lmerControl(optimizer ='Nelder_Mead'))
# summary(l1)
# plot_model(l1)
# 
# # ALPHA Coefficients to LaTeX table
# temp = summary(l1)
# coefs <- as.data.frame(temp$coefficients)
# coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
# rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post")
# 
# # ALPHA Post-hoc tests to LaTeX table
# temp = emmeans(l1, list(pairwise ~ stage), adjust = "tukey")
# ph1 <- as.data.frame(temp$`pairwise differences of stage`)
# ph1 <- ph1[, -4] # remove df since they are at inf
# ph1[,5] = ifelse(ph1[,5] > .05, paste(round(ph1[,5],digits=2),sep=""), ifelse(ph1[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(ph1[,5] < .001,"<.001\\textsuperscript{**}", ifelse(ph1[,5] < .01, "<.01\\textsuperscript{*}", "<.05"))))
# 
# # Concatenate in one LaTeX table
# coefs           <- data.frame(Predictor = row.names(coefs), coefs);
# rownames(coefs) <- NULL
# colnames(ph1)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","z", "\\textit{p}")
# colnames(coefs) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)", "df", "z","\\textit{p}")
# stats_alpha     <- bind_rows(coefs,ph1)
# 
# kbl(stats_alpha, "latex", booktabs = T, linesep = "", label = 'stats_alpha',
#     escape = FALSE, digits = 2,
#     caption = "Effect of sleep stages on Alpha power")%>%
#   pack_rows("Coefficients", 1, 7) %>%
#   pack_rows("Post-hoc comparisons", 8, 28) %>%
#   kable_styling(latex_options = c("HOLD_position"))%>%
#   save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_alpha.tex")


###################################################################
# Mixed model with both sleep stage and power explaining IED rate #
###################################################################

# to create p-values
detach(package:lmerTest)
library(lmerTest)
library(lme4)

# determine reference level
data_pow_wide$hyplabel = relevel(data_pow_wide$stage, ref="Pre")
# l1 <- lmer(IEDsum ~ stage + duration + Zdelta + Ztheta + Zalpha + (1 | night) + (1 | patient), data_pow_wide)
# l1 <- lmer(IEDsum ~ stage + Zdelta + Ztheta + Zalpha + (1 | night) + (1 | patient), data_pow_wide)

l1 <- lmer(IEDsum ~ stage + Zdelta1 + Zdelta2 + (1 | night) + (1 | patient), data_pow_wide)
# l1 <- lmer(IEDsum ~ stage + (1 | night) + (1 | patient), data_pow_wide) # COMPARE MODELS AS THIS SHOWS DELTA EXPLAINS S2 vs S3

summary(l1)
plot_model(l1)

# # get mathematical description of the model and write to latex
# eq <- equatiomatic::extract_eq(l1)
# fileConn<-file("D:/Dropbox/Apps/Overleaf/Hspike/formula/model1.tex")
# writeLines(c("$$",eq,"$$"), fileConn)
# close(fileConn)

# Coefficients to LaTeX table
temp = summary(l1)
coefs <- as.data.frame(temp$coefficients)
coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post", "Delta", "Theta", "Alpha")

# Post-hoc tests to LaTeX table
temp = emmeans(l1, list(pairwise ~ stage), adjust = "tukey")
ph1 <- as.data.frame(temp$`pairwise differences of stage`)
ph1 <- ph1[, -4] # remove df since they are at inf
ph1[,5] = ifelse(ph1[,5] > .05, paste(round(ph1[,5],digits=2),sep=""), ifelse(ph1[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(ph1[,5] < .001,"<.001\\textsuperscript{**}", ifelse(ph1[,5] < .01, "<.01\\textsuperscript{*}", "<.05"))))

# Concatenate in one LaTeX table
coefs           <- data.frame(Predictor = row.names(coefs), coefs);
rownames(coefs) <- NULL
colnames(ph1)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","z", "\\textit{p}")
colnames(coefs) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)", "df", "z","\\textit{p}")
stats_IEDsum   <- bind_rows(coefs,ph1)

kbl(stats_IEDsum, "latex", booktabs = T, linesep = "", label = 'stats_IEDsum',
    escape = FALSE, digits = 2,
    caption = "Effect of sleepstages on IEDrate")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  pack_rows("Sleep stages", 1, 7) %>% # latex_gap_space = "2em"
  pack_rows("Power", 8, 9) %>%
  pack_rows("Post-hoc comparisons", 10, 30) %>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_IEDsum.tex")


###############################
# IED amplitude & sleep stage #
###############################

# prepare data
data_amp <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/amplitude_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_amp <- data_amp[!data_amp$hyplabel == "NO_SCORE", ]
data_amp <- data_amp[!data_amp$hyplabel == "NO_SCORE", ]
data_amp$hyplabel[data_amp$hyplabel == "PHASE_1"]   = "S1"
data_amp$hyplabel[data_amp$hyplabel == "PHASE_2"]   = "S2"
data_amp$hyplabel[data_amp$hyplabel == "PHASE_3"]   = "S3"
data_amp$hyplabel[data_amp$hyplabel == "AWAKE"]     = "Wake"
data_amp$hyplabel[data_amp$hyplabel == "PRESLEEP"]  = "Pre"
data_amp$hyplabel[data_amp$hyplabel == "POSTSLEEP"] = "Post"
data_amp$hyplabel <- factor(data_amp$hyplabel, levels = c("REM", "Post", "Pre", "Wake", "S1", "S2", "S3"))
data_amp$Patient  <- factor(data_amp$patient, levels = c(8:1))
data_amp$part     <- factor(data_amp$part)

# average over label for plotting
posamp   <- setNames(aggregate(data_amp$posamp, by = c(list(data_amp$patient, data_amp$part, data_amp$hyplabel)), mean), c("Patient", "part", "hyplabel", "posamp"))
negamp   <- setNames(aggregate(data_amp$negamp, by = c(list(data_amp$patient, data_amp$part, data_amp$hyplabel)), mean), c("Patient", "part", "hyplabel", "negamp"))
data_amp_mean <- merge(posamp,negamp)

y1   <- setNames(aggregate(data_amp_mean$posamp, by = c(list(data_amp_mean$Patient)), mean), c("Patient", "Mposamp"))
y1sd <- setNames(aggregate(data_amp_mean$posamp, by = c(list(data_amp_mean$Patient)), sd),   c("Patient", "SDposamp"))
y2   <- setNames(aggregate(data_amp_mean$negamp, by = c(list(data_amp_mean$Patient)), mean), c("Patient", "Mnegamp"))
y2sd <- setNames(aggregate(data_amp_mean$negamp, by = c(list(data_amp_mean$Patient)), sd),   c("Patient", "SDnegamp"))

data_amp_mean <- merge(data_amp_mean,y1)
data_amp_mean <- merge(data_amp_mean,y2)
data_amp_mean <- merge(data_amp_mean,y1sd)
data_amp_mean <- merge(data_amp_mean,y2sd)

data_amp_mean$Zposamp <- (data_amp_mean$posamp-data_amp_mean$Mposamp)/data_amp_mean$SDposamp
data_amp_mean$Znegamp <- (data_amp_mean$negamp-data_amp_mean$Mnegamp)/data_amp_mean$SDnegamp
data_amp_mean$Patient <- as.factor(data_amp_mean$Patient)
data_amp_mean$hyplabel <- factor(data_amp_mean$hyplabel, levels = c("REM", "Post", "Pre", "Wake", "S1", "S2", "S3"))

# plot
data_amp_mean$title1 = "Normalized amplitude positive peaks"
data_amp_mean$title2 = "Normalized amplitude negative peaks"

plot_amp_pos <- ggplot(data=data_amp_mean, aes(y=hyplabel, x=Zposamp)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, part), col = Patient), position=position_dodge(width=0.3), alpha = 0.5) +
  # guides(colour = "none") +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  # coord_cartesian(xlim = c(0, 1500)) +
  facet_wrap(~title1)

plot_amp_neg <- ggplot(data=data_amp_mean, aes(y=hyplabel, x=Znegamp)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, part), col = Patient), position=position_dodge(width=0.3), alpha = 0.5) +
  # guides(colour = "none") +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  # coord_cartesian(xlim = c(0, 1500)) +
  facet_wrap(~title2)

###############################
# Boxplots together in Figure #
###############################

ggarrange(plots_delta1, plot_IED_count, plots_delta2, plot_IED_rate, plot_amp_neg, plot_amp_pos, 
          ncol = 2, nrow = 3,
          vjust = 1.5, hjust = -1, 
          labels = c("A","B","C","D","E","F","G"), 
          legend = "right", 
          common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/all_boxplots.pdf")

##############################
## STATISTICS Positive peaks #
##############################

# Normalize 
y1   <- setNames(aggregate(data_amp$posamp, by = c(list(data_amp$patient)), mean), c("patient", "Mposamp"))
y2   <- setNames(aggregate(data_amp$negamp, by = c(list(data_amp$patient)), mean), c("patient", "Mnegamp"))
y1sd <- setNames(aggregate(data_amp$posamp, by = c(list(data_amp$patient)), sd),   c("patient", "SDposamp"))
y2sd <- setNames(aggregate(data_amp$negamp, by = c(list(data_amp$patient)), sd),   c("patient", "SDnegamp"))
data_amp <- merge(data_amp,y1)
data_amp <- merge(data_amp,y2)
data_amp <- merge(data_amp,y1sd)
data_amp <- merge(data_amp,y2sd)
data_amp$Zposamp <- (data_amp$posamp-data_amp$Mposamp)/data_amp$SDposamp
data_amp$Znegamp <- (data_amp$negamp-data_amp$Mnegamp)/data_amp$SDnegamp

# Rename for plotting
data_amp$night <- factor(data_amp$part)

# Reorder for table
data_amp$stage <- factor(data_amp$hyplabel, levels =  c("Pre","S3", "S2", "S1", "Wake", "Post", "REM"))

# fit model
data_amp$stage = relevel(data_amp$stage, ref="Pre")

lpos <- lmer(Zposamp ~ stage + (1 | night) + (1 | patient), data_amp, control = lmerControl(optimizer ='Nelder_Mead'))

summary(lpos)
plot_model(lpos)

# Mathematical description of the model and write to latex
# eq <- equatiomatic::extract_eq(lpos)
# fileConn<-file("D:/Dropbox/Apps/Overleaf/Hspike/formula/model_posamp.tex")
# writeLines(c("$$",eq,"$$"), fileConn)
# close(fileConn)

# Model coefficients to table for LaTeX
temp = summary(lpos)
spos            <- temp$coefficients
spos            <- data.frame(Predictor = row.names(spos), spos);
rownames(spos)  <- NULL
colnames(spos)  <- c("Predictor","Estimate","SD","df","t","p")

# Post-hoc tests
temp = emmeans(lpos, list(pairwise ~ stage), adjust = "tukey")
phpos <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phpos) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phpos$df <- NA

##############################
## STATISTICS Negative peaks #
##############################

lneg <- lmer(Znegamp ~ stage + (1 | night) + (1 | patient), data_amp, control = lmerControl(optimizer ='Nelder_Mead'))

summary(lneg)
plot_model(lneg)

# Mathematical description of the model and write to latex
# eq <- equatiomatic::extract_eq(lneg)
# fileConn<-file("D:/Dropbox/Apps/Overleaf/Hspike/formula/model_negamp.tex")
# writeLines(c("$$",eq,"$$"), fileConn)
# close(fileConn)

# Coefficients for table
temp = summary(lneg)
sneg            <- temp$coefficients
sneg            <- data.frame(Predictor = row.names(sneg), sneg)
colnames(sneg)  <- c("Predictor","EstimateNeg","SDNeg","dfNeg","tNeg","pNeg")

# Coefficients to LaTeX table (and reorder)
sneg$id  <- 1:nrow(sneg)
coef <- merge(sneg, spos)
coef <- coef[order(coef$id), ]
coef <- coef[, c(-1, -7)] 
rownames(coef)  <- NULL

coef[,3] = round(coef[,3],digits=0)
coef[,8] = round(coef[,8],digits=0)

coef$Predictor  <- c("\\textit{Intercept}","S3", "S2", "S1", "Wake", "Post", "REM")
coef <- coef[, c(11,1,2,3,4,5,6,7,8,9,10)]

# Post-hoc comparisons to LaTeX table 
temp = emmeans(lneg, list(pairwise ~ stage), adjust = "tukey")
phneg <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phneg) <- c("Comparison","EstimateNeg","SENeg","dfNeg","Z ratioNeg","pNeg")

phneg$id  <- 1:nrow(phneg)
phneg$dfNeg <- NA

ph <- merge(phneg,phpos)
ph <- ph[order(ph$id), ]
ph <- ph[, -7] 
rownames(ph)  <- NULL

# Concatenate in one LaTeX table
colnames(ph)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")
colnames(coef) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")

stats_amp      <- bind_rows(coef,ph)
colnames(stats_amp) <- c("", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")

stats_amp[,6] = ifelse(stats_amp[,6] > .05, paste(round(stats_amp[,6],digits=2),sep=""), ifelse(stats_amp[,6] < .0001, "<.0001\\textsuperscript{***}", ifelse(stats_amp[,6] < .001,"<.001\\textsuperscript{**}", ifelse(stats_amp[,6] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
stats_amp[,11] = ifelse(stats_amp[,11] > .05, paste(round(stats_amp[,11],digits=2),sep=""), ifelse(stats_amp[,11] < .0001, "<.0001\\textsuperscript{***}", ifelse(stats_amp[,11] < .001,"<.001\\textsuperscript{**}", ifelse(stats_amp[,11] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))

options(knitr.kable.NA = '')

kbl(stats_amp, "latex", booktabs = T, linesep = "", label = 'stats_amp',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stage on ERP peak amplitude")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  pack_rows("Sleep stages", 1, 7) %>% # latex_gap_space = "2em"
  pack_rows("Post-hoc comparisons", 8, 28) %>%
  add_header_above(c(" ", "Negative peaks" = 5, "Positive peaks" = 5)) %>%
  kable_styling(latex_options = c("scale_down"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_amp.tex")

#########
# Units #
#########

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



