#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(ggrepel)
library(sinaplot)
library(ggforce)
require(cowplot)

x_nudging_normal <- function(label){
  if (label == 'V987P'){return (-1.2)}
  if (label == 'A899P'){return (-0.2)}
  if (label == 'D994E'){return (0.5)}
  if (label == 'D994Q'){return (-1.2)}
  if (label == 'Q1005R'){return (-0.7)}
  if (label == 'D950N'){return (-0.5)}
  if (label == 'Q954H'){return (-0.8)}
  if (label == 'T1027I'){return (-0.1)}
  if (label == 'Q1005R'){return (-1)}
  else{return (0)}
  }

y_nudging_normal <- function(label){
  if (label == 'S982A'){return (-0.1)}
  if (label == 'V987P'){return (-0.1)}
  if (label == 'L981F'){return (0.5)}
  if (label == 'A899P'){return (-1)}
  if (label == 'D994E'){return (-1.5)}
  if (label == 'D994Q'){return (-1.3)}
  if (label == 'D950N'){return (0)}
  if (label == 'Q954H'){return (-0.8)}
  if (label == 'Q1005R'){return (-0.6)}
  else{return (0)}
  }

x_nudging_residual <- function(label){
  if (label == 'D950N'){return (-1)}
  if (label == 'V987P'){return (-0.5)}
  if (label == 'D994E'){return (0.5)}
  if (label == 'D994Q'){return (-1)}
  if (label == 'Q1005R'){return (-1)}
  else{return (0)}
  }

y_nudging_residual <- function(label){
  if (label == 'Q1005R'){return (-0.5)}
  if (label == 'D950N'){return (-0.5)}
  if (label == 'V987P'){return (0.5)}
  if (label == 'D994E'){return (-0.5)}
  if (label == 'D994Q'){return (-0.5)}
  else{return (0)}
  }

plot_exp_vs_fus <- function(df, df_special, df_others, df_WT, graphname, ylab){
  textsize <- 7
  if (ylab == 'Fusion score'){
    label_table <- df_special %>%
                     mutate(nudge_x=mapply(x_nudging_normal, mut)) %>%
                     mutate(nudge_y=mapply(y_nudging_normal, mut))
    }
  if (ylab == 'Adjusted fusion score'){
    label_table <- df_special %>%
                     mutate(nudge_x=mapply(x_nudging_residual, mut)) %>%
                     mutate(nudge_y=mapply(y_nudging_residual, mut))
    }
  p <-  ggplot() +
          geom_point(data=df_others, aes(x=exp_score, y=fus_score, size=log10(avg_freq), color='Others'), alpha=0.4, pch=16) +
          geom_point(data=df_WT, aes(x=exp_score, y=fus_score, size=log10(avg_freq), color='WT'), alpha=0.7, pch=16) +
          geom_point(data=df_special, aes(x=exp_score, y=fus_score, size=log10(avg_freq), color=grouping), alpha=0.7, pch=16) +
          scale_color_manual('',
                             values=c('Others'='grey',
                                      'S-Closed/HexaPro'='#92DDAE',
                                      'Non-fusogenic'='#7775D1',
                                      'Fusogenic'='#BC4E47',
                                      'Variants'=qualpal(n = 5, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex[1],
                                      'WT'='#DFBBE7'),
                             drop=FALSE) +
          geom_text_repel(data=df_special, aes(x=exp_score, y=fus_score,label=mut),
                          color="black", min.segment.length=0, segment.size=0.2, size=2, force=15, force_pull=1,
                          seed=5, nudge_x=label_table$nudge_x, nudge_y=label_table$nudge_y, max.overlaps = Inf) +
          geom_smooth(data=df, mapping = aes(x=exp_score, y=fus_score), method = "lm", se=FALSE, color="black", size = 0.2, lty=2) +
          scale_size_continuous(expression(bold(log["10"]~'(avg freq)')),
                                limits = c(-4.5, -1), range = c(0.1, 2.5), breaks = c(-4, -3, -2, -1)) +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),,
                axis.title.y=element_text(size=textsize,face="bold"),
                legend.position = "right",
                legend.title    = element_text(size=textsize,face="bold"),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ylab(ylab) +
          xlab("Expression score")
  ggsave(graphname,p,width=3.2, height=2.5, bg='white', dpi=1200)
  }

coloring <- function(mut){
  if (mut %in% c("K986P", "V987P", "A892P", "A899P", "A942P")){return ("S-Closed/HexaPro")}
  else if (mut %in% c("D994E", "D994Q", "T961F", "Q1005R")){return ("Non-fusogenic")}
  else if (mut %in% c("S943H","A944S")){return ("Fusogenic")}
  else if (mut %in% c("WT")){return ("WT")}
  else if (mut %in% c("D950N", "Q954H", "N969K", "L981F", "S982A", "T1027I")){return ("Variants")}
  else {return ('Others')}
  }

##### MAIN #####
group_level <- c('S-Closed/HexaPro','Non-fusogenic', 'Fusogenic', 'Variants')

df <- read_tsv('result/S2_HR1_DMS_scores.tsv') %>%
        mutate(grouping=factor(mapply(coloring,mut),levels=c('Others',group_level, 'WT'))) %>%
        filter(mut_class %in% c('missense', 'WT'))
Model <- lm(df$fus_score~df$exp_score)
df$residual <- Model$resid
df_special <- filter(df, grouping %in% group_level)
df_WT <- filter(df, grouping=='WT')
df_others <- filter(df, grouping=='Others')
print (paste("correlation between fusion score and expression score: ", cor(df$exp_score, df$fus_score)))
plot_exp_vs_fus(df, df_special, df_others, df_WT, 'graph/Exp_vs_fus.png', 'Fusion score')

df <- mutate(df, fus_score=residual)
df_special <- mutate(df_special, fus_score=residual)
df_WT <- mutate(df_WT, fus_score=residual)
df_others <- mutate(df_others, fus_score=residual)
plot_exp_vs_fus(df, df_special, df_others, df_WT, 'graph/Exp_vs_fus_residual.png', 'Adjusted fusion score')
