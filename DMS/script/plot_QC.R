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
library(viridis)
library(qualpalr)
library(sinaplot)
library(ggforce)
require(cowplot)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
  }

plot_replicate_cor <- function(df, graphname, param, rep1_label, rep2_table){
  print (paste('correlation for:', graphname, cor(df$rep1, df$rep2)))
  textsize <- 7
  df$density <- get_density(df$rep1, df$rep2, n = 100)
  p <- ggplot(df,aes(x=rep1, y=rep2, color=density)) +
         geom_hex(bins = 70) +
         scale_fill_continuous(type = "viridis") +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
         labs(x=bquote(bold(paste(.(param),.(rep1_label)))),y=bquote(bold(paste(.(param),.(rep2_table)))))
  ggsave(graphname, p, height=2, width=2.5)
  }

plot_by_class <- function(df, graphname, ylab){
  df <- df %>%
          filter(mut_class != 'WT')
  textsize <- 7
  p <- ggplot(df,aes(x=mut_class, y=score, group=mut_class)) +
         geom_violin(width=1, color="black") +
         geom_sina(pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='black', alpha=0.2) +
         geom_boxplot(width=0.3, color="black", outlier.shape=NA, alpha=0) + 
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_blank(),
               axis.title.y=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
         ylab(ylab)
  ggsave(graphname, p, height=2, width=2)
  }

t_test <- function(df, class_1, class_2){
  p_value_exp <- t.test(filter(df, mut_class==class_1)$exp_score, filter(df, mut_class==class_2)$exp_score)$p.value
  p_value_fus <- t.test(filter(df, mut_class==class_1)$fus_score, filter(df, mut_class==class_2)$fus_score)$p.value
  print (paste("Exp: p-value of diff between", class_1, 'vs', class_2, ':', p_value_exp))
  print (paste("Fus: p-value of diff between", class_1, 'vs', class_2, ':', p_value_fus))
  }

df <- read_tsv('result/S2_HR1_DMS_scores.tsv')
print (nrow(df))
df_exp12 <- df %>%
            rename(rep1=exp_score_rep1) %>%
            rename(rep2=exp_score_rep2) %>%
            rename(score=exp_score)
df_exp13 <- df %>%
            rename(rep1=exp_score_rep1) %>%
            rename(rep2=exp_score_rep3) %>%
            rename(score=exp_score)
df_exp23 <- df %>%
            rename(rep1=exp_score_rep2) %>%
            rename(rep2=exp_score_rep3) %>%
            rename(score=exp_score)
df_fus <- df %>%
            rename(rep1=fus_score_rep1) %>%
            rename(rep2=fus_score_rep2) %>%
            rename(score=fus_score)
plot_replicate_cor(df_exp12, 'graph/QC_replicate_exp12.png', "Expression score", '(replicate 1)', '(replicate 2)')
plot_replicate_cor(df_exp13, 'graph/QC_replicate_exp13.png', "Expression score", '(replicate 1)', '(replicate 3)')
plot_replicate_cor(df_exp23, 'graph/QC_replicate_exp23.png', "Expression score", '(replicate 2)', '(replicate 3)')
plot_by_class(df_exp12, 'graph/Exp_by_class.png', 'Expression score')
plot_replicate_cor(df_fus, 'graph/QC_replicate_fus.png', "Fusion score", '(replicate 1)', '(replicate 2)')
plot_by_class(df_fus, 'graph/Fus_by_class.png', 'Fusion score')
t_test(df, 'silent', 'nonsense')
t_test(df, 'silent', 'missense')
t_test(df, 'missense', 'nonsense')
#write.table(select(filter(df_exp, mut_class=='missense'), mut, Input_freq, score), 'NTD_DMS_expression_score.tsv', quote=FALSE, sep="\t", row.names=FALSE)
