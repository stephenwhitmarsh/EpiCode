if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

setTimeLimit(100000); setSessionTimeLimit(10000)

library(ggplot2)
library("cowplot")
library("gridExtra")
library(ggpubr)
library(plyr)

# load data
data <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/pnh/spikestats_revision.csv", sep=',', header=TRUE, dec='.', na.strings = " ")

# remove empty lines, dont know why they are there
# data = data[1:25,]

data$U <- factor(data$unit)
data$N <- factor(data$nodule)
# data$SU <- data$percRPV < 0.25
# data$SU <- factor(ifelse(data$percRPV > 0.25, "MUA", "SUA"))
data$SU <- c("MUA","SUA","MUA","MUA","MUA","MUA","SUA","MUA","MUA","MUA","MUA","SUA","MUA","MUA","SUA","MUA","MUA","MUA","SUA","MUA","SUA","SUA","MUA","SUA","MUA")
# data$PI <- factor(data$PI)
data$PI <- factor(ifelse(data$template_tp < 550, "Int", "Pyr"))
# data$PI <- revalue(data$PI, c("Int"="Interneuron", "Pyr"="Pyramidal cell"))



max_y_range <- max(data$FR, data$CV, data$CV2, na.rm = TRUE) - 
  min(data$FR, data$CV, data$CV2, na.rm = TRUE)

ds = (max(data$FR, na.rm = TRUE) - min(data$FR, na.rm = TRUE))/max_y_range/1

# ds = 0.4
## PLOT ALL UNITS

m_all <- ggplot(data, aes(x=1, y=FR)) + 
#  geom_violin() + 
  ggtitle("FR (Hz)") +
  theme_bw() + theme(
    axis.text.y = element_text(colour = 'black'),
    axis.ticks.y = element_line(colour = 'black'),    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  # geom_text(label=data$unit) +
  scale_fill_manual(values=c("#FFFFFF", "#000000")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())


cv1_all <- ggplot(data, aes(x=1, y=CV) ) + 
  #  geom_violin() +
  ggtitle("CV") +
  theme_bw() + theme(
    axis.text.y = element_text(colour = 'black'),
    axis.ticks.y = element_line(colour = 'black'),    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  scale_fill_manual(values=c("#FFFFFF", "#000000")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())


cv2_all <- ggplot(data, aes(x=1, y=CV2) ) + 
#  geom_violin() + 
  ggtitle("CV2") +
  theme_bw() + theme(
    axis.text.y = element_text(colour = 'black'),
    axis.ticks.y = element_line(colour = 'black'),    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  scale_fill_manual(values=c("#FFFFFF", "#000000")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())


## SPLIT FOR CELLTYPE

m_split <- ggplot(data, aes(x=1, y=FR)) + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") + 
  ggtitle("FR (Hz)") +
  # geom_violin() + 
  theme_bw() + theme(
    axis.text.y = element_text(colour = 'black'),
    axis.ticks.y = element_line(colour = 'black'),    
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  scale_fill_manual(values=c("#FFFFFF", "#000000")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())

cv1_split <- ggplot(data, aes(x=1, y=CV) ) + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") +
  ggtitle("CV") +
  #  geom_violin() +
  theme_bw() + theme(
    axis.text.y = element_text(colour = 'black'),
    axis.ticks.y = element_line(colour = 'black'),    
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  scale_fill_manual(values=c("#FFFFFF", "#000000")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())

cv2_split <- ggplot(data, aes(x=1, y=CV2) ) + facet_wrap( ~ PI, scales="fixed", strip.position = "bottom") + 
  # geom_violin() + 
  ggtitle("CV2") +
  theme_bw() + theme(
    axis.text.y = element_text(colour = 'black'),
    axis.ticks.y = element_line(colour = 'black'),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x =element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
  labs(y = "") +
  geom_dotplot(aes(fill=SU), binaxis = "y", stackdir = "centerwhole", dotsize = ds, method = "histodot", stackgroups=TRUE) +
  scale_fill_manual(values=c("#FFFFFF", "#000000")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())


ggarrange(m_all, cv1_all, cv2_all, m_split, cv1_split, cv2_split,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2, align = "v", common.legend = TRUE, legend="bottom")
ggsave("//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/pnh/spikestats_combined_revision.pdf") 


ggarrange(m_split, cv1_split, cv2_split,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 2, align = "v", common.legend = TRUE, legend="none")
ggsave("//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/pnh/spikestats_splitonly_revision.pdf") 





ggplot(data, aes(x=template_width, y=template_tp, color = SU, label = as.numeric(row.names(data)))) + 
  # geom_point(size = 8) +
  # theme_bw() +
  # labs(x = "Peak Width", y = "Trough-to-Peak") +
  theme(legend.position = "none", legend.title = element_blank(), panel.grid.minor = element_blank()) + geom_text()

ggsave("//lexport/iss01.charpier/analyses/stephen.whitmarsh/images/pnh/scatter.pdf") 



