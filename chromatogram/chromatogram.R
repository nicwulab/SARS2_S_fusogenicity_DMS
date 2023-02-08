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

df <- read.csv(file = 'chroma.csv', header = TRUE)

chroma <- ggplot() +
  geom_line(data = df, aes(x = Vol_D994Q, y = mAU_D994Q, color = "D994Q"), size = 1.5, alpha = 0.8) +
  geom_line(data = df, aes(x = Vol_Q1005R, y = mAU_Q1005R, color = "Q1005R"), size = 1.5, alpha = 0.8) +
  geom_line(data = df, aes(x = Vol_QR, y = mAU_QR, color = "QR"), size = 1.5, alpha = 0.8) +
  geom_line(data = df, aes(x = Vol_2P, y = mAU_2P, color = "2P"), size = 1.5, alpha = 0.8) +
  geom_line(data = df, aes(x = Vol_2PQ, y = mAU_2PQ, color = "2PQ"), size = 1.5, alpha = 0.8) +
  geom_line(data = df, aes(x = Vol_2PR, y = mAU_2PR, color = "2PR"), size = 1.5, alpha = 0.8) +
  geom_line(data = df, aes(x = Vol_2PQR, y = mAU_2PQR, color = "2PQR"), size = 1.5, alpha = 0.8) +
  scale_color_brewer(palette="Set2") +
  theme_cowplot(12) +
  xlab("Volume (mL)") +
  ylab("Absorance (mAU)") +
  coord_cartesian(xlim = c(57, 97)) + 
  expand_limits(x = c(55, 100)) + 
  geom_vline(xintercept = 74.7, linetype="dotted", color = "grey50", size=1.5) +
  theme(plot.title=element_text(size=26,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
        axis.text=element_text(size=26,face="bold",colour = 'black'),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
        axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
        axis.title.x=element_text(size=26,face="bold"),
        axis.title.y=element_text(size=26,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title=element_blank(),
        legend.key.width = unit(1.5, "line"))

ggsave("chroma.pdf", plot = chroma, width = 9, height = 6)