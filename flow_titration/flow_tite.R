library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(tidyquant)
require(cowplot)

#df <- read.csv(file = 'hr1_tite.csv', header = TRUE)
#df <- read.csv(file = 'rbd_tite.csv', header = TRUE)
df <- read.csv(file = 'ntd_tite.csv', header = TRUE)
#df <- read.csv(file = 's2sh_tite.csv', header = TRUE)

df[,5:6] <- df[,5:6]/10000

p <- ggplot(data=df, mapping = aes(x = Concentration, y = avg, fill = Mutant, shape = Mutant)) +
  geom_point(aes(col = Mutant), size = 4, alpha = 0.7) +
  geom_errorbar(aes(x = Concentration, ymin=avg-sd, ymax = avg+sd, col = Mutant), width = 0.05, size = 2, alpha = 0.7) +
  # ylim(0, 5) +
  scale_x_log10() +
  annotation_logticks(base = 10, sides = "b", outside = TRUE, short = unit(0.1, "cm"), mid = unit(0.17, "cm"), long = unit(0.25, "cm"), size = 1, color = "black") +
  coord_cartesian(clip = "off") +
  labs(x = expression(bold(paste("Concentration (µg/mL)"))), y = expression(bold(paste("MFI (\u00D7", "10"^{"4"}, ")")))) +
  ggtitle("S2M28 (NTD)") +
  theme_cowplot(12) +
  theme(plot.title=element_text(size=18+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
        axis.text=element_text(size=18,face="bold",colour = 'black'),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black', margin = margin(t=5)),
        axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
        axis.title.x=element_text(size=18,face="bold"),
        axis.title.y=element_text(size=18,face="bold"),
        axis.line = element_line(colour = "black", size = 1),
        axis.ticks.y.left = element_line(colour = "black", size = 1),
        axis.ticks.length.y.left =unit(.25, "cm"),
        legend.title=element_text(size=15+2,face="bold"),
        legend.text=element_text(size=15,face="bold"))

#ggsave("hr1_tite.pdf", width = 5, height = 4, plot = p)
#ggsave("rbd_tite.pdf", width = 5, height = 4, plot = p)
ggsave("ntd_tite.pdf", width = 5, height = 4, plot = p)
#ggsave("s2sh_tite.pdf", width = 5, height = 4, plot = p)