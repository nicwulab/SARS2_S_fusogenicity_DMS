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

df <- read.csv(file = 'fus3.csv', header = TRUE)
df_pre <- df[which(df$Name %in% c("WT", "T961F", "K986P", "V987P", "D994E", "D994Q", "Q1005R")), ]
df_post <- df[which(df$Name %in% c("WT", "S943H", "A944S")), ]
df_combo <- df[which(df$Name %in% c("WT", "2P", "2PQ", "2PR", "2PQR")), ]

pre <- ggplot(data=df_pre, mapping = aes(x = Name, y = Mean, fill = Name)) +
  geom_bar(stat = 'identity', aes(x = Name, y = Mean), width = 0.9) +
  geom_errorbar(aes(x = Name, ymin=MFI_rep2, ymax = MFI_rep1), width = 0.2, size = 1) +
  scale_x_discrete(limits = c("WT", "T961F", "K986P", "V987P", "D994E", "D994Q", "Q1005R")) +
  scale_fill_manual(values = c("#DFBBE7", "#7775D1", "#7775D1", "#7775D1", "#7775D1", "#7775D1","#7775D1"), breaks = c("WT", "T961F", "K986P", "V987P", "D994E", "D994Q", "Q1005R")) +
  ylim(0, 3.2) +
  xlab("") +
  ylab("Fold change in mNG2+ %") +
  theme_cowplot(12) +
  theme(plot.title=element_text(size=28+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
        axis.text=element_text(size=28,face="bold",colour = 'black'),
        axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
        axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=28,face="bold"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "none")

post <- ggplot(data=df_post, mapping = aes(x = Name, y = Mean, fill = Name)) +
  geom_bar(stat = 'identity', aes(x = Name, y = Mean), width = 0.9) +
  geom_errorbar(aes(x = Name, ymin=MFI_rep2, ymax = MFI_rep1), width = 0.2, size = 1) +
  scale_x_discrete(limits = c("WT", "S943H", "A944S")) +
  scale_fill_manual(values = c("#DFBBE7", "#BC4E47", "#BC4E47"), breaks = c("WT", "S943H", "A944S")) +
  ylim(0, 3.2) +
  xlab("") +
  ylab("Fold change in mNG2+ %") +
  theme_cowplot(12) +
  theme(plot.title=element_text(size=28+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
        axis.text=element_text(size=28,face="bold",colour = 'black'),
        axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
        axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=28,face="bold"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "none")

combo <- ggplot(data=df_combo, mapping = aes(x = Name, y = Mean, fill = Name)) +
  geom_bar(stat = 'identity', aes(x = Name, y = Mean), width = 0.9) +
  geom_errorbar(aes(x = Name, ymin=MFI_rep2, ymax = MFI_rep1), width = 0.2, size = 1) +
  scale_x_discrete(limits = c("WT", "2P", "2PQ", "2PR", "2PQR")) +
  scale_fill_manual(values = c("#DFBBE7", "#92DDAE", "#92DDAE", "#92DDAE", "#92DDAE"), breaks = c("WT", "2P", "2PQ", "2PR", "2PQR")) +
  ylim(0, 3.2) +
  xlab("") +
  ylab("Fold change in mNG2+ %") +
  theme_cowplot(12) +
  theme(plot.title=element_text(size=28+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
        axis.text=element_text(size=28,face="bold",colour = 'black'),
        axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
        axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=28,face="bold"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length=unit(.25, "cm"),
        legend.position = "none")

ggsave("fus3_pre.pdf", width = 10, height = 8, plot = pre)
ggsave("fus3_post.pdf", width = 5, height = 8, plot = post)
ggsave("fus3_combo.pdf", width = 8, height = 8, plot = combo)