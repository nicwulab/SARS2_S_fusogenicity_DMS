#R code
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

frequency_data <- read.csv(file = 'SARS_CoV2_S_S2_sequence_conservation.csv', header = TRUE)
expression_scores_by_resi <- read_tsv('S2_HR1_DMS_scores_by_resi_exp.tsv')
fusion_scores_by_resi <- read_tsv('S2_HR1_DMS_scores_by_resi_fus.tsv')
frequency_data$mean_exp_score <- expression_scores_by_resi$mean_exp_score
frequency_data$mean_fus_score <- fusion_scores_by_resi$mean_fus_score

plot_freq_vs_score <- function(df, graphname, score_type){
  textsize <- 12
  if (score_type == "exp"){
  p <- ggplot(df,aes(x=align_freq, y=mean_exp_score)) +
    geom_point(size=1,pch=16, alpha=0.5) +
    #geom_violin(size=0.5) +
    #geom_boxplot(width=0.3, size = 0.3, color="black", outlier.shape=NA, alpha=0.5) +
    #geom_beeswarm(size=0.3, pch=16, alpha=0.5) +
    theme_cowplot(12) +
    theme(plot.title=element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=textsize-1,face="bold"),
          legend.position='right') +
    scale_x_continuous(limit=c(0,1),breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))+
    labs(x=bquote(bold(paste('Sequence conservation'))),y=bquote(bold(paste('Mean expression Score'))))
  ggsave(graphname, p, height=3, width=3)
  }
  if(score_type == "fus"){
    p <- ggplot(df,aes(x=align_freq, y=mean_fus_score)) +
      geom_point(size=1,pch=16, alpha=0.5) +
      #geom_violin(size=0.5) +
      #geom_boxplot(width=0.3, size = 0.3, color="black", outlier.shape=NA, alpha=0.5) +
      #geom_beeswarm(size=0.3, pch=16, alpha=0.5) +
      theme_cowplot(12) +
      theme(plot.title=element_blank(),
            plot.background = element_rect(fill = "white"),
            axis.title=element_text(size=textsize,face="bold"),
            axis.text=element_text(size=textsize,face="bold"),
            legend.key.size=unit(0.1,'in'),
            legend.spacing.x=unit(0.03, 'in'),
            legend.title=element_blank(),
            legend.text=element_text(size=textsize-1,face="bold"),
            legend.position='right') +
      scale_x_continuous(limit=c(0,1),breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))+
      labs(x=bquote(bold(paste('Sequence conservation'))),y=bquote(bold(paste('Mean fusion Score'))))
    ggsave(graphname, p, height=3, width=3)
  }}

plot_freq_vs_score(frequency_data,"dms_expression_vs_seq_conservation_VOC_VOI.pdf", "exp")
plot_freq_vs_score(frequency_data,"dms_fusion_vs_seq_conservation_VOC_VOI.pdf", "fus")
print(cor(frequency_data$align_freq, frequency_data$mean_exp_score,method="spearman", use="complete.obs"))
print(cor(frequency_data$align_freq, frequency_data$mean_fus_score,method="spearman", use="complete.obs"))