
####################
# HYPNOGRAM PLOTTING
####################

# load data
hyp <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/hypnogram_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
hyp$label[hyp$label == "PRESLEEP"] = "AWAKE"
hyp$label[hyp$label == "POSTSLEEP"] = "AWAKE"
hyp$label <- factor(hyp$label, levels = c("PHASE_1", "PHASE_2", "PHASE_3", "AWAKE", "REM"))
hyp$part = factor(hyp$part)

# try plotting with Circlize
position <- read.csv(file="//lexport/iss01.charpier/analyses/stephen.whitmarsh/data/hspike/offset_table.txt", sep=',', header=TRUE, dec='.', na.strings = " ")
position$part = factor(position$part)

track_height = 0.1 * 1.1 ^ c(1:8*3) 
plot_size    = 1 *   1.1 ^ c(1:8*3) 
dev.off()

isize = 1
for(ipatient in 1:8) {
  for(ipart in 1:3) {
    
    hyp_sel <- hyp[hyp$patient == ipatient & hyp$part == ipart, ]
    hyp_sel$part  = factor(hyp_sel$part)
    hyp_sel$dummy = factor(ipart);
    
    circos.par("track.height" = track_height[isize], 
               "start.degree" = 90+360-position_sel$offset[isize]*360, 
               "gap.after" = 360-position_sel$duration[isize]*360,
               "canvas.xlim" = c(-plot_size[isize], plot_size[isize]), 
               "canvas.ylim" = c(-plot_size[isize], plot_size[isize]))
    circos.initialize(hyp_sel$dummy, x = hyp_sel$dt)
    circos.track(hyp_sel$part, ylim = c(1, 5))
    circos.lines(hyp_sel$dt, as.numeric(hyp_sel$label, type = 's' ))
    circos.clear()
    
    par(new = TRUE) # <- magic
    isize = isize + 1;
  }
}

# plot with ggplot
hyp_inc <- hyp
hyp_inc$ID <- factor(as.numeric(hyp$patient) * 10 + as.numeric(hyp$part))

for(ipatient in 1:8) {
  
  for(ipart in 1:3) {
    hyp_inc[hyp_inc$patient == ipatient & hyp_inc$part == ipart, ]$Y = 
      hyp_inc[hyp_inc$patient == ipatient & hyp_inc$part == ipart, ]$Y + (ipatient-1) * 6 * 3 + (ipart-1) * 6 
  }
}

ggplot(hyp_inc, aes(x = hour0, y = Y, group = interaction(part, patient), color = as.factor(patient))) + 
  geom_step() +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  coord_polar(theta = "x", start = 0) +
  ylim(-50, 150)

