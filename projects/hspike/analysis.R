install.packages("circular")
install.packages("ggplot2")
install.packages("units")
install.packages("reshape2")
install.packages("circlize")
install.packages("ggthemes")
install.packages("lemon")
install.packages("egg")
install.packages("readxl")
install.packages('Rcpp')
install.packages("equatiomatic")
install.packages("kableExtra")
install.packages("CircStats")
install.packages("emmeans")
install.packages("lmerTest")
install.packages("devtools")
install.packages("sjPlot")
install.packages("ggpubr")


#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#setTimeLimit(100000); setSessionTimeLimit(10000)
#devtools::install_github("strengejacke/sjPlot")

# ggpval

# install.packages("spiralize")
# library(spiralize)
library(ggplot2)
#library("cowplot")
#library("gridExtra")
library(plyr)
library(reshape2)
library(RColorBrewer)
#library("ggthemes")
#library(lemon)
# library("sjPlot")
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
library("sjPlot") # plot_model


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
    PER: Percntanel, 
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
data[,2]  = round(data[,2],digits=0)
data[,3]  = round(data[,3],digits=0)
data[,4]  = round(data[,4],digits=1)
data[,5]  = round(data[,5],digits=1)
data[,6]  = round(data[,6],digits=0)
data[,7]  = round(data[,7],digits=1)
#data[,8]  = round(data[,8],digits=1)
#data[,9]  = round(data[,9],digits=1)
#data[,10] = round(data[,10],digits=0)
#data[,11] = round(data[,11],digits=1)


kbl(data, "latex", booktabs = T, linesep = "", label = 'performance', escape = FALSE,
    col.names = c("Patient","24hrs", "24hrs","Hit (\\%)","FA (\\%)", "Total", "Total hrs."),
    # col.names = c("Patient","24hrs", "24hrs","Hit (\\%)","FA (\\%)", "Total","24hrs","Hit (\\%)","FA (\\%)", "Total", "Total hrs."),
    caption = "Template detection performance")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
 # kable_styling(latex_options = c("scale_down"))%>%
  row_spec(8, hline_after = TRUE)%>%
  add_header_above(c(" " = 1, "Visual" = 1, "Automatic detection" = 4)) %>%
#  add_header_above(c(" " = 1, "Visual" = 1, "All templates" = 4, "Selected templates" = 4)) %>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/performance.tex")


########################################
# Latex table: Electrode locations MNI #
########################################

options(knitr.kable.NA = '')
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

options(knitr.kable.NA = '')

# normalize by time spend in sleep stages
data <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/hypnogram_duration.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data <- data[, c("patient", "part", "PHASE_3", "PHASE_2", "PHASE_1", "AWAKE", "REM")] 
colnames(data) = c("Patient", "Night","S3", "S2", "S1", "Wake", "REM")

df <- c("\\textit{Mean}", NA, mean(data$S3), mean(data$S2), mean(data$S1), mean(data$Wake), mean(data$REM))
data <- rbind(data, df)
df <- c("\\textit{Std.}", NA, sd(data$S3), sd(data$S2), sd(data$S1), sd(data$Wake), sd(data$REM))
data <- rbind(data, df)

clean.cols <- c("Patient")
data[clean.cols] <- lapply(data[clean.cols], cleanf)

kbl(data, "latex", booktabs = T, linesep = "", label = 'stageduration', escape = FALSE,
    caption = "Time spend in sleep stages (hrs.)", digits=2) %>%
    kable_styling(latex_options = c("HOLD_position"))%>%
    add_header_above(c(" " = 2, "Sleep stage" = 5)) %>%
    row_spec(24, hline_after = TRUE)%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stageduration.tex")


####################################
# Latex table: number of SUA / MUA #
####################################

options(knitr.kable.NA = '')

data <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/DataMUASUA.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data <- data[, c("PatientNr", "Part", "nrSUA", "nrMUA")]
colnames(data) = c("Patient", "Night","SUA", "MUA")
df <- c("\\textit{Sum}", NA, sum(data$SUA), sum(data$MUA))
data <- rbind(data, df)

clean.cols <- c("Patient")
data[clean.cols] <- lapply(data[clean.cols], cleanf)

kbl(data, "latex", booktabs = T, linesep = "", label = 'SUAMUA', escape = FALSE,
    caption = "Number of putatively isolated single units (SUA) and number of multiunits (MUA)")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  row_spec(24, hline_after = TRUE)%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/SUAMUA.tex")


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


# for finding the median for plotting
# data_power$bin <- as.integer(cut(data_power$minute, seq(0, 60*24, by = 1)))
# data_binned_fine <- setNames(aggregate(data_power$power, c(list(data_power$Patient), list(data_power$bin), list(data_power$band)), mean), c("Patient", "bin", "band", "power"))
# data_binned_fine <- as.data.frame(data_binned_fine %>% group_by(Patient, band) %>% mutate(Npower = (power-min(power))/max(power-min(power)))) 
# ggplot(data=data_binned_fine[data_binned_fine$band == "Delta1" & data_binned_fine$Patient == 7, ], aes(x = bin, y = power, fill=Patient, col=Patient)) + 
  geom_smooth()

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
# ggarrange(IEDpolarplot, Seizurepolarplot, 
#          labels = c("A","B"), 
#          vjust = 15, hjust = -1, 
#          legend = "right", 
#          common.legend = TRUE, 
#          font.label = list(size = 14, color = "black", face = "bold")) %>%
#  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/polar_density_seizures.pdf") 


# save combined to pdf
ggarrange(IEDpolarplot, Seizurepolarplot, polardelta1power, polardelta2power,
          labels = c("A","B","C", "D"), 
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
data_pow$hyplabel <- factor(data_pow$hyplabel, levels = c("Pre", "Post", "REM", "Wake", "S1", "S2", "S3"))
data_pow$band     <- factor(data_pow$band, ordered = TRUE, levels = c("Delta1", "Delta2"))
data_pow$patient  <- factor(data_pow$patient, levels = c(8:1))
data_pow$part     <- factor(data_pow$part)

# relative to pre-sleep
temp            <- setNames(aggregate(data_pow$power, by = list(data_pow$patient, data_pow$part, data_pow$band, data_pow$hyplabel), mean), c("patient", "part", "band", "hyplabel", "Pre"))
temp            <- temp[temp$hyplabel=="Pre", ]
temp            <- subset(temp, select = -c(hyplabel))
data_pow        <- merge(data_pow, temp)
data_pow$Zpower <- (data_pow$power-data_pow$Pre) / (data_pow$power+data_pow$Pre)
data_pow$Zpower <- (data_pow$power/data_pow$Pre)
data_pow_sel    <- data_pow[!data_pow$hyplabel=="Pre", ]

# averages for plotting
data_pow_avg    <- setNames(aggregate(data_pow_sel$Zpower, by = list(data_pow_sel$patient, data_pow_sel$part, data_pow_sel$band, data_pow_sel$hyplabel), mean), c("patient", "part", "band", "hyplabel", "power_avg"))
data_pow_avg$patient       <- factor(data_pow_avg$patient, levels = c(8:1))
data_pow_avg$part          <- factor(data_pow_avg$part)

# plot
data_pow_sel$title1 = "Delta1 (0.1-2.5 Hz) power (ratio vs. Pre-Sleep)"
data_pow_sel$title2 = "Delta2 (2.5-4.0 Hz) power (ratio vs. Pre-Sleep)"

plot_delta1 <- 
  ggplot(data=data_pow_sel[data_pow_sel$band == "Delta1",], aes(y = hyplabel, x = Zpower)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_point(data = data_pow_avg[data_pow_avg$band == "Delta1",], aes(group = interaction(patient, part), x = power_avg, y = hyplabel, col = patient),
             position=position_dodge(width=0.5)) +
  guides(colour = "none") +
  coord_cartesian(xlim = c(0, 20)) +
  scale_x_continuous(breaks=c(0, 20)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  theme_article() +
  ylab(NULL) + xlab(NULL) +
  facet_wrap(~title1)

plot_delta2 <- 
  ggplot(data=data_pow_sel[data_pow_sel$band == "Delta2",], aes(y = hyplabel, x = Zpower)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_point(data = data_pow_avg[data_pow_avg$band == "Delta2",], aes(group = interaction(patient, part), x = power_avg, y = hyplabel, col = patient),
             position=position_dodge(width=0.5)) +
  guides(colour = "none") +
  coord_cartesian(xlim = c(0, 10)) +
  scale_x_continuous(breaks=c(0, 10)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  theme_article() +
  ylab(NULL) + xlab(NULL) +
  facet_wrap(~title2)


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
data_IED$hyplabel <- factor(data_IED$hyplabel, levels = c("Pre", "Post", "REM", "Wake", "S1", "S2", "S3"))
data_IED$Patient  <- factor(data_IED$patient, levels = c(8:1))
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
colnames(data_duration) = c('Patient','part','hyplabel','duration')

# count IEDs per sleepstage
data_IEDrate <- na.omit(data_IED %>% dplyr::count(Patient, part, hyplabel))  

# normalize by time spend in sleep stages
data_IEDrate <- merge(data_IEDrate, data_duration)
data_IEDrate$IEDrate = data_IEDrate$n / data_IEDrate$duration / 60 # count is per 10 seconds

# normalize rate by Pre rate 
i <- data_IEDrate[data_IEDrate$hyplabel == 'Pre',]
i <- i[, c('Patient','part','IEDrate')]
colnames(i) = c('Patient','part','Prerate')
data_IEDrate <- merge(data_IEDrate, i)
data_IEDrate$IEDrateNorm <- data_IEDrate$IEDrate - data_IEDrate$Prerate
data_IEDrate_sel <- data_IEDrate[!data_IEDrate$hyplabel=="Pre", ]

# plot
data_IEDrate_sel$title1 = "IED count"
data_IEDrate_sel$title2 = "IED rate (count/minute)"
data_IEDrate_sel$title3 = "IED rate - IED rate Pre-sleep (count/min)"

plot_IED_count <- ggplot(data=data_IEDrate_sel, aes(y=hyplabel, x=n)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, part), col = Patient), position=position_dodge(width=0.5)) +
  # guides(colour = "none") +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  coord_cartesian(xlim = c(1, 3000)) +
  scale_x_continuous(breaks=c(0, 3000)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  theme(legend.position="right") +
  facet_wrap(~title1)

plot_IED_rate <- ggplot(data=data_IEDrate_sel, aes(y=hyplabel, x=IEDrate)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, part), col = Patient), position=position_dodge(width=0.5)) +
  # guides(colour = "none") +
  ylab(NULL) + xlab(NULL) +
  theme_article() + 
  coord_cartesian(xlim = c(0, 20)) +
  scale_x_continuous(breaks=c(0, 20)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  theme(legend.position="right") +
  facet_wrap(~title2)

plot_IED_norm <- ggplot(data=data_IEDrate_sel, aes(y=hyplabel, x=IEDrateNorm)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, part), col = Patient), position=position_dodge(width=0.5)) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.text.y = element_blank()) +
  theme_article() +
  coord_cartesian(xlim = c(0, 20)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  scale_x_continuous(breaks=c(0, 20)) +
  # theme(legend.position ="bottom") +
  # guides(colour = guide_legend(ncol = 1)) +
  # guides(fill=guide_legend(title="Patient")) +
  theme(legend.position="right") +
  labs(color='Patient') +
  facet_wrap(~title3)


#####################################
# IED rate explained by sleep stage #
#####################################

# to create p-values
detach(package:lmerTest)
library(lmerTest)
library(lme4)

# determine reference level
data_IEDrate$hyplabel <- factor(data_IEDrate$hyplabel, levels = c("Pre","S3", "S2", "S1", "Wake", "REM", "Post"))
data_IEDrate$hyplabel = relevel(data_IEDrate$hyplabel, ref="Pre")

lIEDrate <- lmer(IEDrate ~ hyplabel + (1 | part) + (1 | Patient), data_IEDrate, control = lmerControl(optimizer ='Nelder_Mead'))
summary(lIEDrate)
plot_model(lIEDrate)

# Coefficients
temp = summary(lDelta1)
coefs <- as.data.frame(temp$coefficients)
coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post")

# Post-hoc tests
temp = emmeans(lDelta1, list(pairwise ~ stage), adjust = "tukey")
phIEDrate <- as.data.frame(temp$`pairwise differences of stage`)
phIEDrate <- phIEDrate[, -4] # remove df since they are at inf
phIEDrate[,5] = ifelse(phIEDrate[,5] > .05, paste(round(phIEDrate[,5],digits=2),sep=""), ifelse(phIEDrate[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(phIEDrate[,5] < .001,"<.001\\textsuperscript{**}", ifelse(phIEDrate[,5] < .01, "<.01\\textsuperscript{*}", "<.05"))))

# Concatenate in one LaTeX table
coefs           <- data.frame(Predictor = row.names(coefs), coefs);
rownames(coefs) <- NULL
colnames(phIEDrate)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","z", "\\textit{p}")
colnames(coefs) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)", "df", "z","\\textit{p}")
stats_IEDrate   <- bind_rows(coefs,phIEDrate)

options(knitr.kable.NA = '')

kbl(stats_IEDrate, "latex", booktabs = T, linesep = "", label = 'stats_IEDrate',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stages on Delta1 (0.1-2.5Hz) power")%>%
  pack_rows("Coefficients", 1, 7) %>%
  pack_rows("Post-hoc comparisons", 8, 28) %>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_IEDrate.tex")


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
data_pow_wide$Patient  <- factor(data_pow_wide$patient, levels = c(8:1))
data_pow_wide$night    <- factor(data_pow_wide$part)

# to create p-values
detach(package:lmerTest)
library(lmerTest)
library(lme4)

###################################
# Delta1 explained by sleep stage #
###################################

lDelta1 <- lmer(Delta1 ~ stage + (1 | night) + (1 | Patient), data_pow_wide, control = lmerControl(optimizer ='Nelder_Mead'))
summary(lDelta1)
plot_model(lDelta1)

# Coefficients
temp = summary(lDelta1)
coefs <- as.data.frame(temp$coefficients)
coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post")

# Post-hoc tests
temp = emmeans(lDelta1, list(pairwise ~ stage), adjust = "tukey")
phDelta1 <- as.data.frame(temp$`pairwise differences of stage`)
phDelta1 <- phDelta1[, -4] # remove df since they are at inf
phDelta1[,5] = ifelse(phDelta1[,5] > .05, paste(round(phDelta1[,5],digits=2),sep=""), ifelse(phDelta1[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(phDelta1[,5] < .001,"<.001\\textsuperscript{**}", ifelse(phDelta1[,5] < .01, "<.01\\textsuperscript{*}", "<.05"))))

# Concatenate in one LaTeX table
coefs           <- data.frame(Predictor = row.names(coefs), coefs);
rownames(coefs) <- NULL
colnames(phDelta1)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","z", "\\textit{p}")
colnames(coefs) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)", "df", "z","\\textit{p}")
stats_delta1   <- bind_rows(coefs,phDelta1)

options(knitr.kable.NA = '')

kbl(stats_delta1, "latex", booktabs = T, linesep = "", label = 'stats_delta1',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stages on Delta1 (0.1-2.5Hz) power")%>%
    pack_rows("Coefficients", 1, 7) %>%
    pack_rows("Post-hoc comparisons", 8, 28) %>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_delta1.tex")


###################################
# Delta2 explained by sleep stage #
###################################

lDelta2 <- lmer(Delta2 ~ stage + (1 | night) + (1 | Patient), data_pow_wide, control = lmerControl(optimizer ='Nelder_Mead'))
summary(lDelta2)
plot_model(lDelta2)

# Coefficients
temp = summary(lDelta2)
coefs <- as.data.frame(temp$coefficients)
coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post")

# Post-hoc tests
temp = emmeans(lDelta2, list(pairwise ~ stage), adjust = "tukey")
phDelta2 <- as.data.frame(temp$`pairwise differences of stage`)
phDelta2 <- phDelta2[, -4] # remove df since they are at inf
phDelta2[,5] = ifelse(phDelta2[,5] > .05, paste(round(phDelta2[,5],digits=2),sep=""), ifelse(phDelta2[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(phDelta2[,5] < .001,"<.001\\textsuperscript{**}", ifelse(phDelta2[,5] < .01, "<.01\\textsuperscript{*}", "<.05"))))

# Concatenate in one LaTeX table
coefs           <- data.frame(Predictor = row.names(coefs), coefs);
rownames(coefs) <- NULL
colnames(phDelta2)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","z", "\\textit{p}")
colnames(coefs) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)", "df", "z","\\textit{p}")
stats_delta2   <- bind_rows(coefs,phDelta2)

options(knitr.kable.NA = '')

kbl(stats_delta2, "latex", booktabs = T, linesep = "", label = 'stats_delta2',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stages on Delta2 (2.5-4 Hz) power")%>%
  pack_rows("Coefficients", 1, 7) %>%
  pack_rows("Post-hoc comparisons", 8, 28) %>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_delta2.tex")


# ###################################################################
# # Mixed model with both sleep stage and power explaining IED rate #
# ###################################################################
# 
# # to create p-values
# detach(package:lmerTest)
# library(lmerTest)
# library(lme4)
# 
# # determine reference level
# data_pow_wide$hyplabel = relevel(data_pow_wide$stage, ref="Pre")
# l1 <- lmer(IEDsum ~ stage + Delta1 + Delta2 + (1 | night) + (1 | patient), data_pow_wide)
# 
# summary(l1)
# plot_model(l1)
# 
# # # get mathematical description of the model and write to latex
# # eq <- equatiomatic::extract_eq(l1)
# # fileConn<-file("D:/Dropbox/Apps/Overleaf/Hspike/formula/model1.tex")
# # writeLines(c("$$",eq,"$$"), fileConn)
# # close(fileConn)
# 
# # Coefficients to LaTeX table
# temp = summary(l1)
# coefs <- as.data.frame(temp$coefficients)
# coefs[,5] = ifelse(coefs[,5] > .05, paste(round(coefs[,5],digits=2),sep=""), ifelse(coefs[,5] < .0001, "<.0001\\textsuperscript{***}", ifelse(coefs[,5] < .001,"<.001\\textsuperscript{**}", ifelse(coefs[,5] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
# rownames(coefs) <- c("\\textit{Intercept}", "S3", "S2", "S1", "Wake", "REM", "Post", "Delta", "Theta", "Alpha")
# 
# # Post-hoc tests to LaTeX table
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
# stats_IEDsum   <- bind_rows(coefs,ph1)
# 
# kbl(stats_IEDsum, "latex", booktabs = T, linesep = "", label = 'stats_IEDsum',
#     escape = FALSE, digits = 2,
#     caption = "Effect of sleepstages on IEDrate")%>%
#   kable_styling(latex_options = c("HOLD_position"))%>%
#   pack_rows("Sleep stages", 1, 7) %>% # latex_gap_space = "2em"
#   pack_rows("Power", 8, 9) %>%
#   pack_rows("Post-hoc comparisons", 10, 30) %>%
#   save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_IEDsum.tex")

#################
## Plot models ##
#################

set_theme(
  base = theme_article(),
  # panel.bordercol = NA
)

plot_IEDrate_model <- plot_model(
  lIEDrate,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  ylim(-3, 9) +
  font_size(labels.x = 9, labels.y = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3")) +
  scale_y_continuous(breaks=c(-2, 8))
plot_IEDrate_model$data$title = "IED rate (count/minute) model"
plot_IEDrate_model <- plot_IEDrate_model + facet_wrap(~title, scales="free_y")

plot_delta1_model <- plot_model(
  lDelta1,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  # ylim(0, 30) +
  font_size(labels.x = 9, labels.y = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3")) +
  scale_y_continuous(breaks=c(0, 250))
plot_delta1_model$data$title = "Delta1 (0.1-2.5 Hz) model"
plot_delta1_model <- plot_delta1_model + facet_wrap(~title, scales="free_y")

plot_delta2_model <- plot_model(
  lDelta2,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  ylim(-1, 41) +
  font_size(labels.x = 9, labels.y = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3")) +
  scale_y_continuous(breaks=c(0, 40))
plot_delta2_model$data$title = "Delta2 (2.5-4.0 Hz) model"
plot_delta2_model <- plot_delta2_model + facet_wrap(~title, scales="free_y")


###############################
# Boxplots together in Figure #
###############################

ggarrange(plot_IED_rate, plot_IEDrate_model, plot_delta1, plot_delta1_model, plot_delta2, plot_delta2_model,
          ncol = 2, nrow = 3,
          vjust = 1.5, hjust = -1, 
          labels = c("A","B","C","D","E","F","G","H"), 
          legend = "right", 
          common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/IEDrate_delta_boxplots.pdf")


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
data_amp$stage    <- factor(data_amp$hyplabel, levels = c("Pre","S3", "S2", "S1", "Wake", "REM", "Post"))
data_amp$Patient  <- factor(data_amp$patient, levels = c(8:1))
data_amp$night    <- factor(data_amp$part)
data_amp$Template <- factor(data_amp$template)

# average over patient/part/sleepstage/template for plotting
posamp   <- setNames(aggregate(data_amp$posamp, by = c(list(data_amp$patient, data_amp$part, data_amp$Template, data_amp$hyplabel)), mean), c("Patient", "part", "Template", "hyplabel", "posamp"))
negamp   <- setNames(aggregate(data_amp$negamp, by = c(list(data_amp$patient, data_amp$part, data_amp$Template, data_amp$hyplabel)), mean), c("Patient", "part", "Template", "hyplabel", "negamp"))
data_amp_mean <- merge(posamp, negamp)

# # relative to pre-sleep
# temp_pos        <- setNames(aggregate(data_amp$posamp, by = list(data_amp$patient, data_amp$part, data_amp$Template, data_amp$hyplabel), mean), c("Patient", "part", "Template", "hyplabel", "Pre_pos"))
# temp_neg        <- setNames(aggregate(data_amp$negamp, by = list(data_amp$patient, data_amp$part, data_amp$Template, data_amp$hyplabel), mean), c("Patient", "part", "Template", "hyplabel", "Pre_neg"))
# temp_pos        <- temp_pos[temp_pos$hyplabel=="Pre", ]
# temp_neg        <- temp_neg[temp_neg$hyplabel=="Pre", ]
# temp_pos        <- subset(temp_pos, select = -c(hyplabel))
# temp_neg        <- subset(temp_neg, select = -c(hyplabel))
# data_amp_mean   <- merge(data_amp_mean, temp_pos)
# data_amp_mean   <- merge(data_amp_mean, temp_neg)
# data_amp_mean$Zposamp <- (data_amp_mean$posamp-data_amp_mean$Pre_pos)
# data_amp_mean$Znegamp <- (data_amp_mean$negamp-data_amp_mean$Pre_neg)
# data_amp_mean$diffamp <- data_amp_mean$Zposamp + data_amp_mean$Znegamp

# Normalized over all trials
y1   <- setNames(aggregate(data_amp_mean$posamp, by = c(list(data_amp_mean$Patient, data_amp_mean$Template)), mean), c("Patient", "Template", "Mposamp"))
y1sd <- setNames(aggregate(data_amp_mean$posamp, by = c(list(data_amp_mean$Patient, data_amp_mean$Template)), sd),   c("Patient", "Template", "SDposamp"))
y2   <- setNames(aggregate(data_amp_mean$negamp, by = c(list(data_amp_mean$Patient, data_amp_mean$Template)), mean), c("Patient", "Template", "Mnegamp"))
y2sd <- setNames(aggregate(data_amp_mean$negamp, by = c(list(data_amp_mean$Patient, data_amp_mean$Template)), sd),   c("Patient", "Template", "SDnegamp"))

data_amp_mean <- merge(data_amp_mean,y1)
data_amp_mean <- merge(data_amp_mean,y2)
data_amp_mean <- merge(data_amp_mean,y1sd)
data_amp_mean <- merge(data_amp_mean,y2sd)

data_amp_mean$Zposamp <- (data_amp_mean$posamp-data_amp_mean$Mposamp)/data_amp_mean$SDposamp
data_amp_mean$Znegamp <- (data_amp_mean$negamp-data_amp_mean$Mnegamp)/data_amp_mean$SDnegamp
data_amp_mean$diffamp <- data_amp_mean$Zposamp - data_amp_mean$Znegamp

data_amp_mean$Patient <- as.factor(data_amp_mean$Patient)
data_amp_mean$hyplabel <- factor(data_amp_mean$hyplabel, levels = c("Pre","Post", "REM", "Wake", "S1", "S2", "S3"))

data_amp_mean_sel    <- data_amp_mean[!data_amp_mean$hyplabel=="Pre", ]

# plot
data_amp_mean_sel$title1 = "Standardized spike amplitude"
data_amp_mean_sel$title2 = "Standardized slow-wave amplitude"
data_amp_mean_sel$title3 = "Difference standardized spike vs. slow-wave amplitude"

plot_amp_pos <- ggplot(data=data_amp_mean_sel, aes(y=hyplabel, x=Zposamp)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, part), col = Patient), position=position_dodge(width=0.5), size = 1) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  facet_wrap(~title1)

plot_amp_neg <- ggplot(data=data_amp_mean_sel, aes(y=hyplabel, x=Znegamp)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, part), col = Patient), position=position_dodge(width=0.5),size = 1) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  coord_cartesian(xlim = c(-3, 3)) +
  # coord_cartesian(xlim = c(-300, 100)) +
  facet_wrap(~title2)

plot_amp_diff <- ggplot(data=data_amp_mean_sel, aes(y=hyplabel, x=diffamp)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, part), col = Patient), position=position_dodge(width=0.5), size = 1) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  coord_cartesian(xlim = c(-5, 5)) +
  scale_x_continuous(breaks=c(-5, 0, 5)) +
  facet_wrap(~title3)


#########################
## STATISTICS LFP peaks #
#########################

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
data_amp$stage    <- factor(data_amp$hyplabel, levels = c("Pre","S3", "S2", "S1", "Wake", "REM", "Post"))
data_amp$Patient  <- factor(data_amp$patient, levels = c(8:1))
data_amp$night    <- factor(data_amp$part)
data_amp$template <- factor(data_amp$template)

# Standardize data 
y1   <- setNames(aggregate(data_amp$posamp, by = c(list(data_amp$patient, data_amp$template)), mean), c("patient", "template", "Mposamp"))
y2   <- setNames(aggregate(data_amp$negamp, by = c(list(data_amp$patient, data_amp$template)), mean), c("patient", "template", "Mnegamp"))
y1sd <- setNames(aggregate(data_amp$posamp, by = c(list(data_amp$patient, data_amp$template)), sd),   c("patient", "template", "SDposamp"))
y2sd <- setNames(aggregate(data_amp$negamp, by = c(list(data_amp$patient, data_amp$template)), sd),   c("patient", "template", "SDnegamp"))
data_amp <- merge(data_amp,y1)
data_amp <- merge(data_amp,y2)
data_amp <- merge(data_amp,y1sd)
data_amp <- merge(data_amp,y2sd)
data_amp$Zposamp <- (data_amp$posamp-data_amp$Mposamp)/data_amp$SDposamp
data_amp$Znegamp <- (data_amp$negamp-data_amp$Mnegamp)/data_amp$SDnegamp
data_amp$Zdiffamp <- data_amp$Zposamp - data_amp$Znegamp

# fit model
data_amp$stage = relevel(data_amp$stage, ref="Pre")

lpos  <- lmer(posamp  ~ stage + (1 | night) + (1 | template) + (1 | patient), data_amp, control = lmerControl(optimizer ='Nelder_Mead'))
lneg  <- lmer(negamp  ~ stage + (1 | night) + (1 | template) + (1 | patient), data_amp, control = lmerControl(optimizer ='Nelder_Mead'))
ldiff <- lmer(Zdiffamp ~ stage + (1 | night) + (1 | template) + (1 | patient), data_amp, control = lmerControl(optimizer ='Nelder_Mead'))

summary(lpos)
plot_model(lpos)
summary(lneg)
plot_model(lneg)
summary(ldiff)
plot_model(ldiff)

# Post-hoc tests
temp = emmeans(lpos, list(pairwise ~ stage), adjust = "tukey")
phpos <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phpos) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phpos$df <- NA

temp = emmeans(lneg, list(pairwise ~ stage), adjust = "tukey")
phneg <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phneg) <- c("Comparison","EstimateNeg","SENeg","dfNeg","Z ratioNeg","pNeg")
phneg$dfNeg <- NA

temp = emmeans(ldiff, list(pairwise ~ stage), adjust = "tukey")
phdiff <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phdiff) <- c("Comparison","EstimateDiff","SEDiff","dfDiff","Z ratioDiff","pDiff")
phdiff$dfDiff <- NA

# Mathematical description of the model and write to latex
# eq <- equatiomatic::extract_eq(lpos)
# fileConn<-file("D:/Dropbox/Apps/Overleaf/Hspike/formula/model_posamp.tex")
# writeLines(c("$$",eq,"$$"), fileConn)
# close(fileConn)


#############################################
# Coefficients to LaTeX table (and reorder) #
#############################################

# Model coefficients for LaTeX table
temp = summary(lpos)
spos            <- temp$coefficients
spos            <- data.frame(Predictor = row.names(spos), spos);
rownames(spos)  <- NULL
colnames(spos)  <- c("Predictor","Estimate","SD","df","t","p")

temp = summary(lneg)
sneg            <- temp$coefficients
sneg            <- data.frame(Predictor = row.names(sneg), sneg)
colnames(sneg)  <- c("Predictor","EstimateNeg","SDNeg","dfNeg","tNeg","pNeg")

temp = summary(lneg)
sdiff           <- temp$coefficients
sdiff           <- data.frame(Predictor = row.names(sdiff), sdiff)
colnames(sdiff) <- c("Predictor","EstimatePos","SDDiff","dfDiff","tDiff","pDiff")

sneg$id  <- 1:nrow(sneg)
coef <- merge(spos, sneg)
coef <- merge(coef, sdiff)
coef <- coef[order(coef$id), ]
coef <- coef[, c(-1, -12)] 
rownames(coef)  <- NULL

coef[, 3] = round(coef[,3],digits=0)
coef[, 8] = round(coef[,8],digits=0)
coef[,13] = round(coef[,8],digits=0)

coef$Predictor  <- c("\\textit{Intercept}","S3", "S2", "S1", "Wake", "Post", "REM")
coef <- coef[, c(16,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)]

# Post-hoc comparisons
phneg$id  <- 1:nrow(phneg)
ph <- merge(phpos,phneg)
ph <- merge(ph,phdiff)
ph <- ph[order(ph$id), ]
ph <- ph[, -12] 
rownames(ph)  <- NULL

# Concatenate in one LaTeX table
colnames(ph)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", 
                    "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", 
                    "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")
colnames(coef) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", 
                    "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", 
                    "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")

stats_amp      <- bind_rows(coef,ph)
colnames(stats_amp) <- c("", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", 
                         "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", 
                         "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")

stats_amp[,6] = ifelse(stats_amp[,6] > .05, paste(round(stats_amp[,6],digits=2),sep=""), 
                       ifelse(stats_amp[,6] < .0001, "<.0001\\textsuperscript{***}", 
                              ifelse(stats_amp[,6] < .001,"<.001\\textsuperscript{**}", 
                                     ifelse(stats_amp[,6] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
stats_amp[,11] = ifelse(stats_amp[,11] > .05, paste(round(stats_amp[,11],digits=2),sep=""), 
                        ifelse(stats_amp[,11] < .0001, "<.0001\\textsuperscript{***}", 
                               ifelse(stats_amp[,11] < .001,"<.001\\textsuperscript{**}", 
                                      ifelse(stats_amp[,11] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
stats_amp[,16] = ifelse(stats_amp[,16] > .05, paste(round(stats_amp[,16],digits=2),sep=""), 
                        ifelse(stats_amp[,16] < .0001, "<.0001\\textsuperscript{***}", 
                               ifelse(stats_amp[,16] < .001,"<.001\\textsuperscript{**}", 
                                      ifelse(stats_amp[,16] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))

options(knitr.kable.NA = '')

kbl(stats_amp, "latex", booktabs = T, linesep = "", label = 'stats_amp',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stage on ERP peak amplitude")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  pack_rows("Sleep stages", 1, 7) %>% # latex_gap_space = "2em"
  pack_rows("Post-hoc comparisons", 8, 28) %>%
  add_header_above(c(" ", "Positive peaks" = 5, "Negative peaks" = 5, "Peak difference" = 5)) %>%
  kable_styling(latex_options = c("scale_down"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/stats_amp.tex")


#################
## Plot models ##
#################

set_theme(
  base = theme_article(),
  # panel.bordercol = NA
)

plot_amp_pos_model <- plot_model(
  lpos,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  font_size(labels.x = 9, labels.y = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3")) +
  scale_y_continuous(breaks=c(-20, 0, 50, 100, 150))
plot_amp_pos_model$data$title = "Model fixed effects"
plot_amp_pos_model <- plot_amp_pos_model + facet_wrap(~title, scales="free_y")

plot_amp_neg_model <- plot_model(
  lneg,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  font_size(labels.x = 9, labels.y = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3")) +
  scale_y_continuous(breaks=c(0, 50, 100, 150))
plot_amp_neg_model$data$title = "Model fixed effects"
plot_amp_neg_model <- plot_amp_neg_model + facet_wrap(~title, scales="free_y")

plot_amp_diff_model <- plot_model(
  ldiff,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  font_size(labels.x = 9, labels.y = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3")) +
  scale_y_continuous(breaks=c(-0.4, -0.2, 0, 0.2, 0.4))
plot_amp_diff_model$data$title = "Model fixed effects"
plot_amp_diff_model <- plot_amp_diff_model + facet_wrap(~title, scales="free_y")


###############################
# Boxplots together in Figure #
###############################

ggarrange(plot_amp_pos, plot_amp_pos_model, plot_amp_neg, plot_amp_neg_model, plot_amp_diff, plot_amp_diff_model,
          ncol = 2, nrow = 3,
          vjust = 1.5, hjust = -1, 
          labels = c("A","B","C","D","E","F","G","H"), 
          legend = "right", 
          common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/amp_boxplots.pdf")


######################################
# Plot normalized PSTH & sleep stage #
######################################

# prepare data
data_psth <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/psth_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_psth$hyplabel <- factor(data_psth$hyplabel, levels = c("Pre", "Post", "REM", "Wake", "S1", "S2", "S3"))
data_psth$Patient  <- factor(data_psth$Patient, levels = c(8:1))
data_psth$part     <- factor(data_psth$part)
data_psth$template <- factor(data_psth$template)
data_psth$unit     <- factor(data_psth$unit)
data_psth$Type     <- factor(data_psth$SUA)

# select responsive units
data_psth <- data_psth[c(data_psth$responsive == 1), ]

# remove unresponsive patient
data_psth <- data_psth[-c(data_psth$Patient == 7), ]

# select SUA
# data_psth <- data_psth[c(data_psth$SUA == 1), ]

posrate   <- setNames(aggregate(data_psth$posrate, by = c(list(data_psth$Patient, data_psth$unit, data_psth$template, data_psth$hyplabel, data_psth$Type)), mean), c("Patient", "unit", "template", "hyplabel", "Type", "posrate"))
negrate   <- setNames(aggregate(data_psth$negrate, by = c(list(data_psth$Patient, data_psth$unit, data_psth$template, data_psth$hyplabel, data_psth$Type)), mean), c("Patient", "unit", "template", "hyplabel", "Type", "negrate"))
data_psth_mean <- merge(posrate, negrate)

y1   <- setNames(aggregate(data_psth_mean$posrate, by = c(list(data_psth_mean$Patient, data_psth_mean$unit, data_psth_mean$template, data_psth_mean$Type)), mean), c("Patient", "unit", "template", "Type", "Mposrate"))
y1sd <- setNames(aggregate(data_psth_mean$posrate, by = c(list(data_psth_mean$Patient, data_psth_mean$unit, data_psth_mean$template, data_psth_mean$Type)), sd),   c("Patient", "unit", "template", "Type", "SDposrate"))
y2   <- setNames(aggregate(data_psth_mean$negrate, by = c(list(data_psth_mean$Patient, data_psth_mean$unit, data_psth_mean$template, data_psth_mean$Type)), mean), c("Patient", "unit", "template", "Type", "Mnegrate"))
y2sd <- setNames(aggregate(data_psth_mean$negrate, by = c(list(data_psth_mean$Patient, data_psth_mean$unit, data_psth_mean$template, data_psth_mean$Type)), sd),   c("Patient", "unit", "template", "Type", "SDnegrate"))

data_psth_mean <- merge(data_psth_mean,y1)
data_psth_mean <- merge(data_psth_mean,y2)
data_psth_mean <- merge(data_psth_mean,y1sd)
data_psth_mean <- merge(data_psth_mean,y2sd)

data_psth_mean$Zposrate <- (data_psth_mean$posrate-data_psth_mean$Mposrate)/data_psth_mean$SDposrate
data_psth_mean$Znegrate <- (data_psth_mean$negrate-data_psth_mean$Mnegrate)/data_psth_mean$SDnegrate
data_psth_mean$diffrate <- data_psth_mean$Zposrate - data_psth_mean$Znegrate

# plot
data_psth_mean$title1 = "Standardized spike count during LFP spike"
data_psth_mean$title2 = "Standardized spike count during slow wave"
data_psth_mean$title3 = "Relative difference"

data_psth_mean_sel    <- data_psth_mean[!data_psth_mean$hyplabel=="Pre", ]


plot_cnt_pos <- ggplot(data=data_psth_mean_sel, aes(y=hyplabel, x=Zposrate)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, unit), col = Patient, shape = Type), position=position_dodge(width=0.5), alpha = 0.2, size = 0.8) +
  geom_point(data=data_psth_mean_sel[data_psth_mean_sel$Type == 1, ], 
             aes(group = interaction(Patient, unit), col = Patient, shape = Type), position=position_dodge(width=0.5), size = 0.8) +
  scale_shape_discrete(label = c("MUA", "SUA"))+
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  # coord_cartesian(xlim = c(0, 1500)) +
  facet_wrap(~title1)
  
plot_cnt_neg <- ggplot(data=data_psth_mean_sel, aes(y=hyplabel, x=Znegrate)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, unit), col = Patient, shape = Type), position=position_dodge(width=0.5), alpha = 0.2, size = 0.8) +
  geom_point(data=data_psth_mean_sel[data_psth_mean_sel$Type == 1, ], 
             aes(group = interaction(Patient, unit), col = Patient, shape = Type), position=position_dodge(width=0.5), size = 0.8) +
  scale_shape_discrete(label = c("MUA", "SUA"))+
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  # coord_cartesian(xlim = c(0, 1500)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  facet_wrap(~title2)

plot_cnt_diff <- ggplot(data=data_psth_mean_sel, aes(y=hyplabel, x=diffrate)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, unit), col = Patient, shape = Type), position=position_dodge(width=0.5), alpha = 0.2, size = 0.8) +
  geom_point(data=data_psth_mean_sel[data_psth_mean_sel$Type == 1, ], 
             aes(group = interaction(Patient, unit), col = Patient, shape = Type), position=position_dodge(width=0.5), size = 0.8) +
  scale_shape_discrete(label = c("MUA", "SUA"))+
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  # coord_cartesian(xlim = c(-2.5, 2.5)) +
  facet_wrap(~title3)


###################################
## STATISTICS Positive peaks PSTH #
###################################

# to create p-values
detach(package:lmerTest)
library(lmerTest)
library(lme4)

# prepare data
data_psth            <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/psth_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_psth$hyplabel   <- factor(data_psth$hyplabel, levels = c("Pre", "Post", "REM", "Wake", "S1", "S2", "S3"))
data_psth$Patient    <- factor(data_psth$Patient,  levels = c(8:1))
data_psth$part       <- factor(data_psth$part)
data_psth$template   <- factor(data_psth$template)
data_psth$unit       <- factor(data_psth$unit)
data_psth$responsive <- factor(data_psth$responsive)
data_psth$diffrate   <- data_psth$posrate - data_psth$negrate

# select responsive units
data_psth <- data_psth[c(data_psth$responsive == 1), ]

# remove unresponsive patient
data_psth <- data_psth[-c(data_psth$Patient == 7), ]

# remove NaN
data_psth <- data_psth[!is.nan(data_psth$posrate), ]

# Rename for plotting
data_psth$night <- factor(data_psth$part)

# Reorder for table
data_psth$stage <- factor(data_psth$hyplabel, levels =  c("Pre","S3", "S2", "S1", "Wake", "REM", "Post"))

# fit model
data_psth$stage = relevel(data_psth$stage, ref="Pre")

lpos  <- lmer(posrate  ~ stage + (1 | Patient) + (1 | template) + (1 | unit), data_psth, control = lmerControl(optimizer ='Nelder_Mead'))
lneg  <- lmer(negrate  ~ stage + (1 | Patient) + (1 | template) + (1 | unit), data_psth, control = lmerControl(optimizer ='Nelder_Mead'))
ldiff <- lmer(diffrate ~ stage + (1 | Patient) + (1 | template) + (1 | unit), data_psth, control = lmerControl(optimizer ='Nelder_Mead'))

summary(lpos)
plot_model(lpos)
summary(lneg)
plot_model(lneg)
summary(ldiff)
plot_model(ldiff)

# Post-hoc tests
temp = emmeans(lpos, list(pairwise ~ stage), adjust = "tukey")
phpos <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phpos) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phpos$df <- NA

temp = emmeans(lneg, list(pairwise ~ stage), adjust = "tukey")
phneg <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phneg) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phneg$df <- NA

temp = emmeans(ldiff, list(pairwise ~ stage), adjust = "tukey")
phdiff <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phdiff) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phdiff$df <- NA

# Post-hoc comparisons to LaTeX table 
temp = emmeans(lneg, list(pairwise ~ stage), adjust = "tukey")
phneg <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phneg) <- c("Comparison","EstimateNeg","SENeg","dfNeg","Z ratioNeg","pNeg")


#########################################
# Model coefficients to table for LaTeX #
#########################################

# Model coefficients for table
temp = summary(lpos)
spos            <- temp$coefficients
spos            <- data.frame(Predictor = row.names(spos), spos);
rownames(spos)  <- NULL
colnames(spos)  <- c("Predictor","Estimate","SD","df","t","p")

temp = summary(lneg)
sneg            <- temp$coefficients
sneg            <- data.frame(Predictor = row.names(sneg), sneg);
rownames(sneg)  <- NULL
colnames(sneg)  <- c("Predictor","Estimate","SD","df","t","p")

temp = summary(ldiff)
sdiff            <- temp$coefficients
sdiff            <- data.frame(Predictor = row.names(sdiff), sdiff);
rownames(sdiff)  <- NULL
colnames(sdiff)  <- c("Predictor","Estimate","SD","df","t","p")

# to LaTeX table (and reorder)
sneg$id  <- 1:nrow(sneg)
coef <- merge(sneg, spos, by="Predictor")
coef <- merge(coef, sdiff, by="Predictor")

coef <- coef[order(coef$id), ]
coef <- coef[, c(-1, -7)] 
rownames(coef)  <- NULL

coef[,3] = round(coef[,3],digits=0)
coef[,8] = round(coef[,8],digits=0)
coef[,8] = round(coef[,13],digits=0)

coef$Predictor  <- c("\\textit{Intercept}","S3", "S2", "S1", "Wake", "Post", "REM")
coef <- coef[, c(16,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)]

# posthoc coefficients
phneg$id  <- 1:nrow(phneg)
phneg$dfNeg <- NA

ph <- merge(phneg,phpos)
ph <- merge(ph,phdiff, by = "Comparison")
ph <- ph[order(ph$id), ]
ph <- ph[, -7] 
rownames(ph)  <- NULL

# Concatenate in one LaTeX table
colnames(ph)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")
colnames(coef) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")

stats_cnt      <- bind_rows(coef,ph)
colnames(stats_cnt) <- c("", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")

stats_cnt[,6]  = ifelse(stats_cnt[,6] > .05, paste(round(stats_cnt[,6],digits=2),sep=""), ifelse(stats_cnt[,6] < .0001, "<.0001\\textsuperscript{***}", ifelse(stats_cnt[,6] < .001,"<.001\\textsuperscript{**}", ifelse(stats_cnt[,6] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
stats_cnt[,11] = ifelse(stats_cnt[,11] > .05, paste(round(stats_cnt[,11],digits=2),sep=""), ifelse(stats_cnt[,11] < .0001, "<.0001\\textsuperscript{***}", ifelse(stats_cnt[,11] < .001,"<.001\\textsuperscript{**}", ifelse(stats_cnt[,11] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
stats_cnt[,16] = ifelse(stats_cnt[,16] > .05, paste(round(stats_cnt[,16],digits=2),sep=""), ifelse(stats_cnt[,16] < .0001, "<.0001\\textsuperscript{***}", ifelse(stats_cnt[,16] < .001,"<.001\\textsuperscript{**}", ifelse(stats_cnt[,16] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))

options(knitr.kable.NA = '')

kbl(stats_cnt, "latex", booktabs = T, linesep = "", label = 'unit_stats',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stage on firing rates during IED")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  pack_rows("Sleep stages", 1, 7) %>% # latex_gap_space = "2em"
  pack_rows("Post-hoc comparisons", 8, 28) %>%
  add_header_above(c(" ", "Slow wave" = 5, "Peak" = 5, "Ratio" = 5)) %>%
  kable_styling(latex_options = c("scale_down"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/unit_stats.tex")


#################
## Plot models ##
#################

set_theme(
  base = theme_article(),
  # panel.bordercol = NA
)

plot_cnt_pos_model <- plot_model(
  lpos,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  ylim(-5, 6) +
  font_size(labels.x = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2))) +
  theme( axis.text.y=element_blank())
plot_cnt_pos_model$data$title = "Model of spike count during LFP spike (vs. Pre-sleep)"
plot_cnt_pos_model$data$title = "Model fixed effects"
plot_cnt_pos_model <- plot_cnt_pos_model + facet_wrap(~title, scales="free_y")

plot_cnt_neg_model <- plot_model(
  lneg,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  ylim(-3, 2) +
  font_size(labels.x = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2))) +
  theme( axis.text.y=element_blank())
plot_cnt_neg_model$data$title = "Model of spike count during slow-wave (vs. Pre-sleep)"
plot_cnt_neg_model$data$title = "Model fixed effects"
plot_cnt_neg_model <- plot_cnt_neg_model + facet_wrap(~title, scales="free_y")

plot_cnt_diff_model <- plot_model(
  ldiff,
  title = "",
  colors = "bw",
  axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  ylim(-5, 7) +
  font_size(labels.x = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2))) +
  theme( axis.text.y=element_blank())
plot_cnt_diff_model$data$title = "Model of spike count difference between LFP spike & slow wave (vs. Pre-sleep)"
plot_cnt_diff_model$data$title = "Model fixed effects"
plot_cnt_diff_model <- plot_cnt_diff_model + facet_wrap(~title, scales="free_y")

ggarrange(plot_cnt_pos, plot_cnt_pos_model, plot_cnt_neg, plot_cnt_neg_model, plot_cnt_diff, plot_cnt_diff_model,
          ncol = 2, nrow = 3,
          vjust = 1.5, hjust = -1, 
          labels = c("A","B","C","D","E","F","G","H"), 
          legend = "right", 
          common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/unit_boxplots.pdf")


#########################################
# Plot correlation between LFP and PSTH #
#########################################

# prepare data
data_rho <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/rho_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_rho$Patient  <- factor(data_rho$Patient, levels = c(8:1))
data_rho$part     <- factor(data_rho$part)
data_rho$template <- factor(data_rho$template)
data_rho$unit     <- factor(data_rho$unit)

# select responsive units
# data_rho <- data_rho[c(data_rho$responsive == 1), ]

plot_rho <- ggplot(data=data_rho, aes(y=Patient, x=corr_rho)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = interaction(Patient, unit), col = Patient), position=position_dodge(width=0.3), alpha = 0.5) +
  # guides(colour = "none") +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  theme(legend.position="bottom") +
  coord_cartesian(xlim = c(-0.6, 0.6)) 

final_plot = ggarrange(plot_rho, 
                       nrow = 1, ncol = 1,
                       vjust = 1.5, hjust = -1, 
                       legend = "bottom", 
                       font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/rho_boxplots.pdf")



##########
Combine:
  
data_psth 
data_amp

###########################
# Unit baseline behaviour #
###########################

data_window            <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/window_spike_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_window$hyplabel   <- factor(data_window$hyplabel, levels = c("Pre", "Post", "REM", "Wake", "S1", "S2", "S3"))
data_window$Patient    <- factor(data_window$patient,  levels = c(8:1))
data_window$Night      <- factor(data_window$part)
data_window$Unit       <- factor(data_window$unit)
data_window$responsive <- factor(data_window$responsive)
data_window$Type       <- factor(data_window$Type)
levels(data_window$Type)[levels(data_window$Type)=="good"] <- "SUA"
levels(data_window$Type)[levels(data_window$Type)=="mua"] <- "MUA"
data_window$Type       <- factor(ordered(data_window$Type, levels = c("MUA", "SUA")))

data_sel <- na.omit(data_window)              # removing rows with missing data
data_sel <- data_sel[data_sel$BAD_cnt == 0, ] # removing rows with artefacts
data_sel <- data_sel[data_sel$IEDsum == 0, ]  # removing windows with IEDs

data_CV2   <- setNames(aggregate(data_sel$CV2_trial,           by = c(list(data_sel$Patient, data_sel$Unit, data_sel$Type, data_sel$hyplabel)), mean), c("Patient", "unit", "Type", "hyplabel", "CV2"))
data_CV2_burst <- setNames(aggregate(data_sel$CV2_intraburst_trial, 
                                     by = c(list(data_sel$Patient, data_sel$Unit, data_sel$Type, data_sel$hyplabel)), mean), c("Patient", "unit", "Type", "hyplabel", "CV2_burst"))
data_FR    <- setNames(aggregate(data_sel$trialfreq,           by = c(list(data_sel$Patient, data_sel$Unit, data_sel$Type, data_sel$hyplabel)), mean), c("Patient", "unit", "Type", "hyplabel", "FR"))
data_FRcor <- setNames(aggregate(data_sel$trialfreq_corrected, by = c(list(data_sel$Patient, data_sel$Unit, data_sel$Type, data_sel$hyplabel)), mean), c("Patient", "unit", "Type", "hyplabel", "FRcor"))
data_burst <- setNames(aggregate(data_sel$burst_trialsum,      by = c(list(data_sel$Patient, data_sel$Unit, data_sel$Type, data_sel$hyplabel)), mean), c("Patient", "unit", "Type", "hyplabel", "Bursts"))
data_amp   <- setNames(aggregate(data_sel$amplitude,           by = c(list(data_sel$Patient, data_sel$Unit, data_sel$Type, data_sel$hyplabel)), mean), c("Patient", "unit", "Type", "hyplabel", "Amplitude"))
data_rpv   <- setNames(aggregate(data_sel$RPV,                 by = c(list(data_sel$Patient, data_sel$Unit, data_sel$Type, data_sel$hyplabel)), mean), c("Patient", "unit", "Type", "hyplabel", "RPV"))
data_mean  <- merge(data_CV2,  data_FR)
data_mean  <- merge(data_mean, data_FRcor)
data_mean  <- merge(data_mean, data_burst)
data_mean  <- merge(data_mean, data_amp)
data_mean  <- merge(data_mean, data_rpv)
data_mean  <- merge(data_mean, data_CV2_burst)

# relative to pre-sleep
temp              <- data_CV2
names(temp)[names(temp) == "CV2"] <- "Pre_CV2"
temp              <- temp[temp$hyplabel=="Pre", ]
temp              <- subset(temp, select = -c(hyplabel))
data_mean         <- merge(data_mean, temp)
data_mean$CV2_rel <- (data_mean$CV2-data_mean$Pre_CV2) / (data_mean$CV2+data_mean$Pre_CV2)

temp              <- data_CV2_burst
names(temp)[names(temp) == "CV2_burst"] <- "Pre_CV2_burst"
temp              <- temp[temp$hyplabel=="Pre", ]
temp              <- subset(temp, select = -c(hyplabel))
data_mean         <- merge(data_mean, temp)
data_mean$CV2_burst_rel <- (data_mean$CV2_burst-data_mean$Pre_CV2_burst) / (data_mean$CV2_burst+data_mean$Pre_CV2_burst)

temp              <- data_FR
names(temp)[names(temp) == "FR"] <- "Pre_FR"
temp              <- temp[temp$hyplabel=="Pre", ]
temp              <- subset(temp, select = -c(hyplabel))
data_mean         <- merge(data_mean, temp)
data_mean$FR_rel <- (data_mean$FR-data_mean$Pre_FR) / (data_mean$FR+data_mean$Pre_FR)

temp              <- data_FRcor
names(temp)[names(temp) == "FRcor"] <- "Pre_FRcor"
temp              <- temp[temp$hyplabel=="Pre", ]
temp              <- subset(temp, select = -c(hyplabel))
data_mean         <- merge(data_mean, temp)
data_mean$FRcor_rel <- (data_mean$FRcor-data_mean$Pre_FRcor) / (data_mean$FRcor+data_mean$Pre_FRcor)

temp              <- data_amp
names(temp)[names(temp) == "Amplitude"] <- "Pre_Amplitude"
temp              <- temp[temp$hyplabel=="Pre", ]
temp              <- subset(temp, select = -c(hyplabel))
data_mean         <- merge(data_mean, temp)
data_mean$Amplitude_rel <- (data_mean$Amplitude-data_mean$Pre_Amplitude) / (data_mean$Amplitude+data_mean$Pre_Amplitude)

temp              <- data_burst
names(temp)[names(temp) == "Bursts"] <- "Pre_Bursts"
temp              <- temp[temp$hyplabel=="Pre", ]
temp              <- subset(temp, select = -c(hyplabel))
data_mean         <- merge(data_mean, temp)
data_mean$Bursts_rel <- (data_mean$Bursts-data_mean$Pre_Bursts) / (data_mean$Bursts+data_mean$Pre_Bursts)

# Normalize per unit
# FRmean <- setNames(aggregate(data_mean_sel$FR_corrected, by = c(list(data_mean_sel$unit)), mean), c("unit", "FRmean"))
# FRsd   <- setNames(aggregate(data_mean_sel$FR_corrected, by = c(list(data_mean_sel$unit)), sd),   c("unit", "FRsd"))
# data_mean_sel <- merge(data_mean_sel,FRmean)
# data_mean_sel <- merge(data_mean_sel,FRsd)
# data_mean_sel$FR_standardized <- (data_mean_sel$FR_corrected-data_mean_sel$FRmean)/data_mean_sel$FRsd
# 
# Ampmean <- setNames(aggregate(data_mean_sel$Amplitude, by = c(list(data_mean_sel$unit)), mean), c("unit", "Ampmean"))
# Ampsd   <- setNames(aggregate(data_mean_sel$Amplitude, by = c(list(data_mean_sel$unit)), sd),   c("unit", "Ampsd"))
# data_mean_sel <- merge(data_mean_sel,Ampmean)
# data_mean_sel <- merge(data_mean_sel,Ampsd)
# data_mean_sel$Amp_standardized <- (data_mean_sel$Amplitude-data_mean_sel$Ampmean)/data_mean_sel$Ampsd


# plot
data_mean_sel <- data_mean[!data_mean$hyplabel == "Pre", ] # Remove for plotting
data_mean_sel$title_CV2   = "CV2 (spike intervals)"
data_mean_sel$title_CV2_burst   = "CV2 (burst intervals)"
data_mean_sel$title_FR    = "Firing rate"
data_mean_sel$title_FRcor = "Firing rate corrected"
data_mean_sel$title_SUA_FRcor = "Firing rate corrected (SUA)"
data_mean_sel$title_MUA_FRcor = "Firing rate corrected (MUA)"
data_mean_sel$title_Amp   = "Amplitude"
data_mean_sel$title_Burst = "Burst rate"

FR_plot <- 
  ggplot(data=data_mean_sel, aes(y=hyplabel, x=FR_rel )) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.8,
             aes(group = interaction(Patient, unit), col = Patient, shape = Type, fill = Patient), 
             position = position_dodge(width = 0.5), alpha = data_mean_sel$alphaval) +
  scale_shape_manual(name = "Unit type", labels = c("MUA", "SUA"), values = c(15, 17)) +
  scale_alpha(guide = 'none') +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  coord_cartesian(xlim = c(-1, 1)) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  facet_wrap(~title_FR)

FRcor_SUA_plot <- 
  ggplot(data=data_mean_sel[data_mean_sel$Type=="SUA", ], aes(y=hyplabel, x=FRcor_rel )) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.8,
             aes(group = interaction(Patient, unit), col = Patient, fill = Patient), 
             position = position_dodge(width = 0.5)) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  coord_cartesian(xlim = c(-1, 1)) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  facet_wrap(~title_SUA_FRcor)

FRcor_MUA_plot <- 
  ggplot(data=data_mean_sel[data_mean_sel$Type=="MUA", ], aes(y=hyplabel, x=FRcor_rel )) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.8,
             aes(group = interaction(Patient, unit), col = Patient, fill = Patient), 
             position = position_dodge(width = 0.5)) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  coord_cartesian(xlim = c(-1, 1)) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  facet_wrap(~title_MUA_FRcor)

Amp_plot <- 
  ggplot(data=data_mean_sel[data_mean_sel$Type=="SUA", ], aes(y=hyplabel, x=Amplitude_rel)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.8,
             aes(group = interaction(Patient, unit), col = Patient, fill = Patient), 
             position = position_dodge(width = 0.5)) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  coord_cartesian(xlim = c(-0.2, 0.2)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  facet_wrap(~title_Amp)

CV2_plot <- 
  ggplot(data=data_mean_sel[data_mean_sel$Type == "SUA", ], aes(y=hyplabel, x=CV2_rel)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.8, # shape = 17,
             aes(group = interaction(Patient, unit), col = Patient, fill = Patient), 
             position = position_dodge(width = 0.5)) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  coord_cartesian(xlim = c(-0.1, 0.1)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  facet_wrap(~title_CV2)

CV2_burst_plot <- 
  ggplot(data=data_mean_sel[data_mean_sel$Type == "SUA", ], aes(y=hyplabel, x=CV2_burst_rel)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.8, # shape = 17,
             aes(group = interaction(Patient, unit), col = Patient, fill = Patient), 
             position = position_dodge(width = 0.5)) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  coord_cartesian(xlim = c(-0.6, 0.6)) +
  scale_x_continuous(breaks = c(-0.6, 0, 0.6)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  facet_wrap(~title_CV2_burst)

Burst_plot <- ggplot(data=data_mean_sel[data_mean_sel$Type == "SUA", ], aes(y=hyplabel, x=Bursts_rel)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.8, # shape = 17,
             aes(group = interaction(Patient, unit), col = Patient, fill = Patient), 
             position = position_dodge(width = 0.5)) +
  # scale_shape_manual(name = "Unit type", labels = c("MUA", "SUA"), values = c(15, 17)) +
  # scale_shape_discrete(label = c("MUA", "SUA")) +
  ylab(NULL) + xlab(NULL) + 
  theme_article() + 
  # coord_cartesian(xlim = c(-1, 1)) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_discrete(expand=expansion(mult=c(0.1,0.2))) +
  facet_wrap(~title_Burst)


################
## STATISTICS ##
################

# to create p-values
detach(package:lmerTest)
library(lmerTest)
library(lme4)

# load and prepare data
data_window            <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/window_spike_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
data_window$hyplabel   <- factor(data_window$hyplabel, levels = c("Pre", "Post", "REM", "Wake", "S1", "S2", "S3"))
data_window$Patient    <- factor(data_window$patient,  levels = c(8:1))
data_window$Night      <- factor(data_window$part)
data_window$Unit       <- factor(data_window$unit)
data_window$responsive <- factor(data_window$responsive)
data_window$Type       <- factor(data_window$Type)
data_window$FR         <- data_window$trialfreq
data_window$FRcor      <- data_window$trialfreq_corrected
levels(data_window$Type)[levels(data_window$Type)=="good"] <- "SUA"
levels(data_window$Type)[levels(data_window$Type)=="mua"] <- "MUA"
data_window$Type       <- ordered(data_window$Type, levels = c("MUA", "SUA"))

data_sel <- na.omit(data_window)              # remove rows with missing data
data_sel <- data_sel[data_sel$BAD_cnt == 0, ] # remove rows with artefacts
data_sel <- data_sel[data_sel$IEDsum == 0, ]  # remove windows with IEDs

# Rename for plotting
data_sel$night <- factor(data_sel$part)

# Reorder for table
data_sel$stage <- factor(data_sel$hyplabel, levels =  c("Pre","S3", "S2", "S1", "Wake", "REM", "Post"))

# fit model
data_sel$stage = relevel(data_sel$stage, ref="Pre")

lFR    <- lmer(FR             ~ stage + (1 | Patient) + (1 | Unit), data_sel[data_sel$Type=="MUA", ])
lFRcor_SUA <- lmer(FRcor      ~ stage + (1 | Patient) + (1 | Unit), data_sel[data_sel$Type=="SUA", ])
lFRcor_MUA <- lmer(FRcor      ~ stage + (1 | Patient) + (1 | Unit), data_sel[data_sel$Type=="MUA", ])
lAmp   <- lmer(amplitude      ~ stage + (1 | Patient) + (1 | Unit), data_sel)
lCV2   <- lmer(CV2_trial      ~ stage + (1 | Patient) + (1 | Unit), data_sel)
lCV2_burst   <- lmer(CV2_intraburst_trial ~ stage + (1 | Patient) + (1 | Unit), data_sel)
lBurst <- lmer(burst_trialsum ~ stage + (1 | Patient) + (1 | Unit), data_sel)

summary(lFR)
summary(lFRcor_SUA)
summary(lFRcor_MUA)
summary(lAmp)
summary(lCV2)
summary(lCV2_burst)
summary(lBurst)
plot_model(lFR)
plot_model(lFRcor_SUA)
plot_model(lFRcor_MUA)
plot_model(lAmp)
plot_model(lCV2)
plot_model(lCV2_burst)
plot_model(lBurst)

# Post-hoc tests
temp = emmeans(lFR, list(pairwise ~ stage), adjust = "tukey")
phFR <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phFR) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phFR$df <- NA

temp = emmeans(lFRcor_SUA, list(pairwise ~ stage), adjust = "tukey")
phFR_SUA <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phFR_SUA) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phFR_SUA$df <- NA

temp = emmeans(lFRcor_MUA, list(pairwise ~ stage), adjust = "tukey")
phFR_MUA <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phFR_MUA) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phFR_MUA$df <- NA

temp = emmeans(lFRcor, list(pairwise ~ stage), adjust = "tukey")
phFRcor <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phFRcor) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phFRcor$df <- NA

temp = emmeans(lAmp, list(pairwise ~ stage), adjust = "tukey")
phAmp <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phAmp) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phAmp$df <- NA

temp = emmeans(lCV2, list(pairwise ~ stage), adjust = "tukey")
phCV2 <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phCV2) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phCV2$df <- NA

temp = emmeans(lBurst, list(pairwise ~ stage), adjust = "tukey")
phBurst <- as.data.frame(temp$`pairwise differences of stage`)
colnames(phBurst) <- c("Comparison","Estimate","SE","df","Z ratio","p")
phBurst$df <- NA

#########################################
# Model coefficients to table for LaTeX #
#########################################

# Model coefficients for table
temp = summary(lFR)
sFR            <- temp$coefficients
sFR            <- data.frame(Predictor = row.names(sFR), sFR);
rownames(sFR)  <- NULL
colnames(sFR)  <- c("Predictor","Estimate","SD","df","t","p")

temp = summary(lFRcor_SUA)
sFRcor_SUA           <- temp$coefficients
sFRcor_SUA           <- data.frame(Predictor = row.names(sFRcor_SUA), sFRcor_SUA);
rownames(sFRcor_SUA) <- NULL
colnames(sFRcor_SUA) <- c("Predictor","Estimate","SD","df","t","p")

temp = summary(lFRcor_MUA)
sFRcor_MUA           <- temp$coefficients
sFRcor_MUA           <- data.frame(Predictor = row.names(sFRcor_MUA), sFRcor_MUA);
rownames(sFRcor_MUA) <- NULL
colnames(sFRcor_MUA) <- c("Predictor","Estimate","SD","df","t","p")

temp = summary(lAmp)
sAmp            <- temp$coefficients
sAmp            <- data.frame(Predictor = row.names(sAmp), sAmp);
rownames(sAmp)  <- NULL
colnames(sAmp)  <- c("Predictor","Estimate","SD","df","t","p")

temp = summary(lCV2)
sCV2             <- temp$coefficients
sCV2            <- data.frame(Predictor = row.names(sCV2), sCV2);
rownames(sCV2)  <- NULL
colnames(sCV2)  <- c("Predictor","Estimate","SD","df","t","p")

temp = summary(lBurst)
sBurst          <- temp$coefficients
sBurst          <- data.frame(Predictor = row.names(sBurst), sBurst);
rownames(sBurst)  <- NULL
colnames(sBurst)  <- c("Predictor","Estimate","SD","df","t","p")

# Model coefficients to LaTeX table (and reorder)
sFRcor_SUA$id <- 1:nrow(sFRcor_SUA)
coef <- merge(sFRcor_SUA, sFRcor_MUA, by="Predictor")
coef <- merge(coef, sBurst, by="Predictor")
coef <- merge(coef, sAmp, by="Predictor")
coef <- merge(coef, sCV2, by="Predictor")

coef <- coef[order(coef$id), ]
coef <- coef[, c(-1, -7)] 
rownames(coef)  <- NULL

coef[,3] = round(coef[,3],digits=0)
coef[,8] = round(coef[,8],digits=0)
coef[,13] = round(coef[,13],digits=0)
coef[,18] = round(coef[,18],digits=0)

coef$Predictor  <- c("\\textit{Intercept}","S3", "S2", "S1", "Wake", "REM", "Post")
coef <- coef[, c(21,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]

# posthoc coefficients

phFRcor$id  <- 1:nrow(phFR)
# phFR$dfNeg <- NA

ph <- merge(phFRcor,phAmp, by = "Comparison")
ph <- merge(ph,phCV2,   by = "Comparison")
ph <- merge(ph,phBurst, by = "Comparison")
ph <- ph[order(ph$id), ]
ph <- ph[, -7] 
rownames(ph)  <- NULL

# Concatenate in one LaTeX table
colnames(ph)   <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")
colnames(coef) <- c("Predictor", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")

stats_window <- bind_rows(coef,ph)
colnames(stats_window) <- c("", "Coef $\\beta$","SE($\\beta$)","df", "z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}", "Coef $\\beta$","SE($\\beta$)","df","z", "\\textit{p}")
stats_window[,6]  = ifelse(stats_window[,6]  > .05, paste(round(stats_window[,6], digits=2),sep=""), ifelse(stats_window[,6]  < .0001,  "<.0001\\textsuperscript{***}",  ifelse(stats_window[,6] < .001,"<.001\\textsuperscript{**}", ifelse(stats_window[,6]  < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
stats_window[,11] = ifelse(stats_window[,11] > .05, paste(round(stats_window[,11],digits=2),sep=""), ifelse(stats_window[,11] < .0001, "<.0001\\textsuperscript{***}", ifelse(stats_window[,11] < .001,"<.001\\textsuperscript{**}", ifelse(stats_window[,11] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
stats_window[,16] = ifelse(stats_window[,16] > .05, paste(round(stats_window[,16],digits=2),sep=""), ifelse(stats_window[,16] < .0001, "<.0001\\textsuperscript{***}", ifelse(stats_window[,16] < .001,"<.001\\textsuperscript{**}", ifelse(stats_window[,16] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))
stats_window[,21] = ifelse(stats_window[,21] > .05, paste(round(stats_window[,21],digits=2),sep=""), ifelse(stats_window[,21] < .0001, "<.0001\\textsuperscript{***}", ifelse(stats_window[,21] < .001,"<.001\\textsuperscript{**}", ifelse(stats_window[,21] < .01, "<.01\\textsuperscript{*}", "<.05\\textsuperscript{.}"))))

options(knitr.kable.NA = '')

kbl(stats_window, "latex", booktabs = T, linesep = "", label = 'unit_window_stats',
    escape = FALSE, digits = 2,
    caption = "Effect of sleep stage on firing rates during IED")%>%
  kable_styling(latex_options = c("HOLD_position"))%>%
  pack_rows("Sleep stages", 1, 7) %>% # latex_gap_space = "2em"
  pack_rows("Post-hoc comparisons", 8, 28) %>%
  add_header_above(c(" ", "Firing Rate" = 5, "Amplitude" = 5, "CV2" = 5, "Bursts" = 5)) %>%
  kable_styling(latex_options = c("scale_down"))%>%
  save_kable("D:/Dropbox/Apps/Overleaf/Hspike/tables/unit_window_stats.tex")


set_theme(
  base = theme_article(),
  # panel.bordercol = NA
  axis.textsize = 0.9
)

FRcor_SUA_plot_model <- plot_model(
  lFRcor_SUA,  
  title = "",
  colors = "bw",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  font_size(labels.y = 9, labels.x = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3"))
 
FRcor_MUA_plot_model <- plot_model(
  lFRcor_MUA,  
  title = "",
  colors = "bw",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5) +
  # axis.labels = c("REM","Post","Wake","S1","S2","S3")) +  
  font_size(labels.y = 9, labels.x = 9) +
  # ylim(-3, 0) + 
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3"))

CV2_plot_model <- plot_model(
  lCV2,  
  title = "",
  colors = "bw",
  # axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5,
  axis.labels = c("REM","Post","Wake","S1","S2","S3")) +  
  font_size(labels.y = 9, labels.x = 9) +
  ylim(-0.01, 0.1) + 
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3"))

CV2_burst_plot_model <- plot_model(
  lCV2_burst,  
  title = "",
  colors = "bw",
  # axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5,
  axis.labels = c("REM","Post","Wake","S1","S2","S3")) +  
  font_size(labels.y = 9, labels.x = 9) +
  ylim(-0.05, 0.05) +
  # scale_y_continuous(breaks=c(-0.05, 0, 0.05)) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3"))

Amp_plot_model <- plot_model(
  lAmp,  
  title = "",
  colors = "bw",
  # axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5,
  axis.labels = c("REM","Post","Wake","S1","S2","S3")) +  
  font_size(labels.y = 9, labels.x = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3"))

Burst_plot_model <- plot_model(
  lBurst,  
  title = "",
  colors = "bw",
  # axis.labels = "",
  axis.title = "",
  show.values = TRUE,
  show.p = TRUE,
  decimals = 4,
  digits = 2,
  value.offset = 0.5,
  value.size = 2.5,
  axis.labels = c("REM","Post","Wake","S1","S2","S3")) +  
  font_size(labels.y = 9, labels.x = 9) +
  scale_x_discrete(expand=expansion(mult=c(0.1,0.2)), 
                   labels=c("Post","REM","Wake","S1","S2","S3"))

Burst_plot_model$data$title = "Model fixed effects"
Burst_plot_model <- Burst_plot_model + facet_wrap(~title, scales="free_y")
FR_plot_model$data$title = "Model fixed effects"
FR_plot_model <- FR_plot_model + facet_wrap(~title, scales="free_y")
FRcor_SUA_plot_model$data$title = "Model fixed effects"
FRcor_SUA_plot_model <- FRcor_SUA_plot_model + facet_wrap(~title, scales="free_y")
FRcor_MUA_plot_model$data$title = "Model fixed effects"
FRcor_MUA_plot_model <- FRcor_MUA_plot_model + facet_wrap(~title, scales="free_y")
Amp_plot_model$data$title = "Model fixed effects"
Amp_plot_model <- Amp_plot_model + facet_wrap(~title, scales="free_y")
CV2_plot_model$data$title = "Model fixed effects"
CV2_plot_model <- CV2_plot_model + facet_wrap(~title, scales="free_y")
CV2_burst_plot_model$data$title = "Model fixed effects"
CV2_burst_plot_model <- CV2_burst_plot_model + facet_wrap(~title, scales="free_y")

ggarrange(FRcor_SUA_plot, FRcor_SUA_plot_model, FRcor_MUA_plot, FRcor_MUA_plot_model, Burst_plot, Burst_plot_model,
          ncol = 2, nrow = 3,
          vjust = 1.5, hjust = -1, 
          labels = c("A","B","C","D","E","F","G","H"), 
          legend = "right", 
          common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/window_FR_boxplots.pdf")


ggarrange(Amp_plot, Amp_plot_model, CV2_plot, CV2_plot_model,  CV2_burst_plot, CV2_burst_plot_model,
          ncol = 2, nrow = 3,
          vjust = 1.5, hjust = -1, 
          labels = c("A","B","C","D","E","F","G","H"), 
          legend = "right", 
          common.legend = TRUE, 
          font.label = list(size = 14, color = "black", face = "bold")) %>%
  ggexport(filename = "D:/Dropbox/Apps/Overleaf/Hspike/images/window_CV2_boxplots.pdf")

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
# 
# ## POLAR ##
# 
# library("ggplot2")
# library("dplyr")
# 
# ## plot power values PER HOUR
# 
# data$bin    <- as.integer(cut(data$minute, seq(0, 24*60, by = 60)))
# data_binned <- setNames(aggregate(data$alpha, c(list(data$patient), list(data$bin)), mean), c("patient", "bin", "alpha"))
# data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Nalpha = alpha/max(alpha))) 
# 
# ggplot(data=data_binned, aes(x = bin, y = alpha, fill=patient, col=patient)) +
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   geom_bar(stat="identity", width = 1) +
#   geom_vline(xintercept = seq(0.5, 21.5, by = 3), colour = "grey90")  + 
#   theme_bw() +
#   theme(panel.border = element_blank(),
#         legend.key = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_blank(),
#         panel.grid  = element_blank(),
#   ) + 
#   scale_x_continuous(breaks=seq(0.5, 21.5, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
#   coord_polar(theta = "x", start = 0)
# 
# 
# ## plot power values (use bindivision for hour divisions)
# bindivision = 1
# data$bin    <- as.integer(cut(data$minute, seq(0, 24*60, by = 60/bindivision)))
# y1 <- setNames(aggregate(data$alpha, c(list(data$patient), list(data$bin)), median), c("patient", "bin", "alpha"))
# y2 <- setNames(aggregate(data$theta, c(list(data$patient), list(data$bin)), median), c("patient", "bin", "theta"))
# y3 <- setNames(aggregate(data$beta,  c(list(data$patient), list(data$bin)), median), c("patient", "bin", "beta"))
# y4 <- setNames(aggregate(data$delta, c(list(data$patient), list(data$bin)), median), c("patient", "bin", "delta"))
# y1sd <- setNames(aggregate(data$alpha, c(list(data$patient), list(data$bin)), sd), c("patient", "bin", "SDalpha"))
# y2sd <- setNames(aggregate(data$theta, c(list(data$patient), list(data$bin)), sd), c("patient", "bin", "SDtheta"))
# y3sd <- setNames(aggregate(data$beta,  c(list(data$patient), list(data$bin)), sd), c("patient", "bin", "SDbeta"))
# y4sd <- setNames(aggregate(data$delta, c(list(data$patient), list(data$bin)), sd), c("patient", "bin", "SDdelta"))
# data_binned <- merge(y1,y2)
# data_binned <- merge(data_binned,y3)
# data_binned <- merge(data_binned,y4)
# data_binned <- merge(data_binned,y1sd)
# data_binned <- merge(data_binned,y2sd)
# data_binned <- merge(data_binned,y3sd)
# data_binned <- merge(data_binned,y4sd)
# data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Nalpha = (alpha-min(alpha)) / max(alpha-min(alpha))))
# data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Ntheta = (theta-min(theta)) / max(theta-min(theta))))  
# data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Ndelta = (delta-min(delta)) / max(delta-min(delta))))
# data_binned <- as.data.frame(data_binned %>% group_by(patient) %>% mutate(Nbeta  = (beta-min(beta))   / max(beta-min(beta))))
# data_binned_melt <-melt(data_binned[c("Ndelta","Ntheta","Nalpha","Nbeta","patient","bin")], id=c("patient","bin"))
# data_binned_melt$variable <- factor(data_binned_melt$variable, labels = c("\u03B4 power (1-4Hz)","\u03B8 power (5-7Hz)","\u03B1 power (8-14Hz)","\u03B2 power (15-25Hz)"))
# data_binned_melt$patient  <- factor(data_binned_melt$patient, levels = c(8:1))
# 
# pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/polar_power_combined.pdf")
# 
# ggplot(data=data_binned_melt, aes(x = bin, y = value, fill=patient, col=patient)) +
#   ggtitle("Normalized circadian LFP power") +
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   geom_vline(xintercept = seq(0.5, 21.5 * bindivision, by = 3 * bindivision), colour = "grey90")  + 
#   geom_hline(yintercept = seq(1, 8, by = 1), colour = "grey90")  + 
#   geom_bar(stat="identity", width = 1) +
#   theme_bw() +
#   theme(panel.border = element_blank(),
#         legend.key = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_blank(),
#         panel.grid  = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         strip.background = element_blank(),
#   ) + 
#   scale_x_continuous(breaks=seq(0.5, 21.5*bindivision, by = 3 * bindivision), 
#                      labels = c("0" = "00:00", "3" = "3:00", "6" = "6:00", "9" = "9:00", 
#                                 "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
#   coord_polar(theta = "x", start = 0) +
#   ylim(-1, 7) +
#   facet_wrap(~variable)
# 
# dev.off()
# 
# library(circular)
# 
# data$radian = data$minute / (24*60) * pi * 2 
# 
# stats_circ <- median.circular(data$radian)
# 
# 
# # load data: Patients x Units x time window
# data_IED <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/hypdata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
# 
# # prepare data
# data_IED$patient      <- factor(data_IED$patient, levels = c(8:1))
# data_IED$hour         <- data_IED$minute / (60)
# data_IED$part         <- factor(data_IED$part,)
# 
# # load data: Patients x Units x time window
# data_seizures <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/seizuredata_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
# 
# # prepare data
# data_seizures$patient  <- factor(data_seizures$patient, levels = c(8:1))
# data_seizures$hour     <- data_seizures$minute / (60)
# 
# # graphics.off()
# 
# # display.brewer.all(colorblindFriendly = TRUE)
# 
# pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/polar_density_combined.pdf")
# 
# ggplot(data=data_IED, aes(x=hour, fill=patient, col=patient)) +
#   ggtitle("Circadian interictal density & seizure occurance") +
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
#   geom_hline(yintercept = seq(-0.25+0.025, -0.05, by = 0.025), colour = "grey90") +
#   geom_histogram(position = 'stack', aes(y = stat(density)), binwidth=0.1, center = 0.05) + # make sure binwidth is multiple of 24/60
#   coord_polar(theta = "x") +
#   scale_x_continuous(breaks=seq(0, 21, by = 3), labels = c("0" = "00:00", "3" = "3:00", "6" = "6:00", "9" = "9:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
#   coord_polar(theta = "x", start = 0, direction = 1) + 
#   theme_bw() +
#   theme(panel.border = element_blank(),
#         legend.key = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_blank(),
#         panel.grid  = element_blank(),
#   ) + 
#   geom_point(size=1, data = data_seizures, aes(x = hour, y = -0.05 - as.numeric(patient)*0.025+0.025, fill = patient, col = patient)) +
#   ylim(-0.3, 0.85)
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(data=data_IED, aes(x=hour, fill=patient, col=patient)) +
#   ggtitle("Circadian interictal density & seizure occurance") +
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
#   #geom_segment(aes(x = 0, xend = 0, y = 1, yend = 0.1)) +
#   geom_hline(yintercept = seq(-0.25+0.025, -0.05, by = 0.025), colour = "grey90") +
#   geom_histogram(position = 'stack', aes(y = stat(density)), binwidth=0.1, center = 0.05) + # make sure binwidth is multiple of 24/60
#   coord_polar(theta = "x") +
#   scale_x_continuous(breaks=seq(0, 21, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
#   coord_polar(theta = "x", start = 0, direction = 1) + 
#   theme_bw() +
#   theme(panel.border = element_blank(),
#         legend.key = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_blank(),
#         panel.grid  = element_blank(),
#   ) + 
#   geom_point(size=0.5, data = data_seizures, aes(x = hour, y = -0.05 - as.numeric(patient)*0.025+0.025, fill = patient, col = patient)) +
#   ylim(-0.3, 0.85)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(circular)
# library(units)
# 
# plot.circular(data_seizures$radian, pch = 16, cex = 1, stack = TRUE,
#               axes = TRUE, start.sep=0, sep = 0.025, shrink = 1,
#               bins = 23, ticks = TRUE, tcl = 0.025, tcl.text = 0.125,
#               col = NULL, tol = 0.04, uin = NULL,
#               xlim = c(-1, 1), ylim = c(-1, 1), digits = 2, units = "rad",
#               template = NULL, zero = 0.1, rotation = NULL,
#               main = NULL, sub=NULL, xlab = "", ylab = "")
# 
# 
# plot(sin, -pi, 2*pi) # see ?plot.function
# 
# pi_rad <- as_units(pi, "radians")
# 
# 
# 
# pdf(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/hspike/R/polar_density_per patient.pdf")
# 
# ggplot(data=data, aes(x=hour, fill=patient, col=patient)) +
#   ggtitle("Circadian interictal density per patient") +
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   geom_histogram(position = 'stack', aes(y = stat(density)), binwidth=0.2, center = 0.1) + # make sure binwidth is multiple of 24/60
#   coord_polar(theta = "x") +
#   geom_vline(xintercept = seq(0, 21, by = 3), colour = "grey90") +
#   scale_x_continuous(breaks=seq(0, 21, by = 3), labels = c("0" = "00:00", "3" = "03:00", "6" = "06:00", "9" = "09:00", "12" = "12:00", "15" = "15:00", "18" = "18:00", "21" = "21:00")) +
#   coord_polar(theta = "x", start = 0, direction = 1) + 
#   theme_bw() +
#   
#   theme(panel.border = element_blank(),
#         legend.key = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_blank(),
#         panel.grid  = element_blank(),
#         strip.background = element_blank(),
#   ) + facet_wrap(~patient)
# dev.off()
# 
# 
# 
# 
# 
# 
# # get ylim and xlim from plot
# #layer_scales(p)$y$range$range
# #layer_scales(p)$x$range$range
# 
# 
# 
# 
# 
# # data$SU <- data$percRPV < 0.25
# data$SU <- factor(ifelse(data$percRPV > 0.25, "MUA", "SUA"))
# data$SU <- c("MUA","SUA","MUA","MUA","MUA","MUA","SUA","MUA","MUA","MUA","MUA","SUA","MUA","MUA","SUA","MUA","MUA","MUA","SUA","MUA","SUA","SUA","MUA","SUA","MUA")
# # data$PI <- factor(data$PI)
# data$PI <- factor(ifelse(data$template_tp < 550, "Int", "Pyr"))
# # data$PI <- revalue(data$PI, c("Int"="Interneuron", "Pyr"="Pyramidal cell"))
# 
# 
# 
# 
# data <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/aurelie/statAH.csv", sep=';', header=TRUE, dec=',', na.strings = " ")
# data$P <- data$ï..Patient
# 
# # get some p-values
# detach(package:lmerTest)
# library(lmerTest)
# l1 <- lmer(EEG ~ NSE +S100 + (1 | P) + (1 | Time), data); summary(l1)
# l1 <- lmer(EEG ~ NSE + S100 + (1 | P), data); summary(l1)
# 
# l1 <- lm(EEG ~ NSE + S100, data); summary(l1)



