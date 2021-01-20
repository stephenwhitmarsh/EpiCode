if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

setTimeLimit(100000); setSessionTimeLimit(10000)

library(ggplot2)
library("cowplot")
library("gridExtra")
library(ggpubr)
library(plyr)

# load data
data <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/pnh/spikestats.csv", sep=',', header=TRUE, dec='.', na.strings = " ")

# remove empty lines, dont know why they are there
# data = data[1:25,]

data$U <- factor(data$unit)
data$N <- factor(data$ipatient)
# data$SU <- data$percRPV < 0.25
data$SU <- factor(ifelse(data$percRPV > 0.25, "MUA", "SUA"))
data$SU <- c("MUA","SUA","MUA","MUA","MUA","MUA","SUA","MUA","MUA","MUA","MUA","SUA","MUA","MUA","SUA","MUA","MUA","MUA","SUA","MUA","SUA","SUA","MUA","SUA","MUA")
# data$PI <- factor(data$PI)
data$PI <- factor(ifelse(data$template_tp < 550, "Int", "Pyr"))
# data$PI <- revalue(data$PI, c("Int"="Interneuron", "Pyr"="Pyramidal cell"))



max_y_range <- max(data$ES_bl_freq, data$ES_bl_fano, data$ES_bl_BI, data$ES_bl_CV, na.rm = TRUE) - 
  min(data$ES_bl_freq, data$ES_bl_fano, data$ES_bl_BI, data$ES_bl_CV, na.rm = TRUE)

ds = (max(data$ES_bl_freq, na.rm = TRUE) - min(data$ES_bl_freq, na.rm = TRUE))/max_y_range/1

## PLOT ALL UNITS

m <- ggplot(data, aes(x=1, y=ES_bl_freq)) + 
#  geom_violin() + 
  theme_bw() + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "Baseline Firingrate (Hz)") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  # geom_text(label=data$unit) +
  theme(legend.title = element_blank())


ff <- ggplot(data, aes(x=1, y=ES_bl_fano) ) + 
#  geom_violin() + 
  theme_bw() + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "Fano Factor") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())


cv1 <- ggplot(data, aes(x=1, y=ES_bl_CV) ) + 
#  geom_violin() +
  theme_bw() + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "CV") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())

bi <- ggplot(data, aes(x=1, y=ES_bl_BI) ) +
#  geom_violin() + 
  theme_bw() + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
   labs(y = "Burst Index") +
   geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  # geom_text(label=data$unit) +
  theme(legend.title = element_blank())

rpv <- ggplot(data, aes(x=1, y=percRPV) ) +
#  geom_violin() + 
  theme_bw() + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "%RPV") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())

ggarrange(rpv, m, ff, cv1, bi,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 5, nrow = 2, align = "v", common.legend = TRUE, legend="bottom")
ggsave("//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/pnh/spikestats_combined2.pdf") 

## SPLIT FOR CELLTYPE

m <- ggplot(data, aes(x=1, y=ES_bl_freq)) + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") + 
 # geom_violin() + 
  theme_bw() + theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "Baseline Firingrate (Hz)") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank()) 

ff <- ggplot(data, aes(x=1, y=ES_bl_fano) ) + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") + 
 # geom_violin() + 
  theme_bw() + theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "Fano Factor") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())


cv1 <- ggplot(data, aes(x=1, y=ES_bl_CV) ) + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") +
#  geom_violin() +
  theme_bw() + theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "CV") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())

bi <- ggplot(data, aes(x=1, y=ES_bl_BI) ) + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") +
#  geom_violin() + 
  theme_bw() + theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "Burst Index") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())

rpv <- ggplot(data, aes(x=1, y=percRPV) ) + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") +
#  geom_violin() + 
  theme_bw() + theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "%RPV") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())

width <- ggplot(data, aes(x=1, y=template_width) )  + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") +
#  geom_violin() + 
  theme_bw() + theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "Half-Width (micro s)") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())

tp <- ggplot(data, aes(x=1, y=template_tp) )  + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") +
 # geom_violin() + 
  theme_bw() + theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "Trough-Peak (micro s)") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  theme(legend.title = element_blank())

e <- plot.new()


ggarrange(e, e, e, width, tp, rpv, m, ff, cv1, bi,
          labels = c("", "", "", "B", "C", "D", "E", "F", "G", "H", "I"),
          ncol = 5, nrow = 2, align = "v", common.legend = TRUE, legend="bottom" )
ggsave("//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/pnh/spikestats_split.pdf", paper = 'a4r') 


ggplot(data, aes(x=template_width, y=template_tp, color = SU, label = as.numeric(row.names(data)))) + 
  # geom_point(size = 8) +
  # theme_bw() +
  # labs(x = "Peak Width", y = "Trough-to-Peak") +
  theme(legend.position = "none", legend.title = element_blank(), panel.grid.minor = element_blank()) + geom_text()

ggsave("//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/pnh/scatter.pdf") 



