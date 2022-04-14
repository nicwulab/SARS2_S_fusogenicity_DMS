library(ggplot2)
library(ggrepel)
library(ggnewscale)

# Load counts file.
dat <- as.data.frame(read.table(file = 'S2_HR1_DMS_count_trimmed_aa.tsv', sep = '\t', header = TRUE))

# Add a pseudocount of 1 to all entries prevent division by 0.
for (i in 3:19){
  dat[, i] <- dat[, i] + 1
}

# For each replicate, for each mutant, total count = count in bin 0 + count in bin 1 + count in bin 2 + count in bin 3 + count in mNG2_neg + count in mNG2_pos.
total_rep1_count <- rowSums(dat[, c(4, 5, 6, 7, 8, 9)])
total_rep2_count <- rowSums(dat[, c(10, 11, 12, 13, 14, 15)])
total_rep3_count <- rowSums(dat[, c(16, 17, 18, 19)])

# Calculate the total number of counts per bin per replicate.
total_input_count <- rep(c(sum(dat$ipt_count)), times = nrow(dat))

total_bin_0_rep1_count <- rep(c(sum(dat$bin_0_rep1_count)), times = nrow(dat))
total_bin_1_rep1_count <- rep(c(sum(dat$bin_1_rep1_count)), times = nrow(dat))
total_bin_2_rep1_count <- rep(c(sum(dat$bin_2_rep1_count)), times = nrow(dat))
total_bin_3_rep1_count <- rep(c(sum(dat$bin_3_rep1_count)), times = nrow(dat))
total_mNG2_neg_rep1_count <- rep(c(sum(dat$mNG2_neg_rep1_count)), times = nrow(dat))
total_mNG2_pos_rep1_count <- rep(c(sum(dat$mNG2_pos_rep1_count)), times = nrow(dat))

total_bin_0_rep2_count <- rep(c(sum(dat$bin_0_rep2_count)), times = nrow(dat))
total_bin_1_rep2_count <- rep(c(sum(dat$bin_1_rep2_count)), times = nrow(dat))
total_bin_2_rep2_count <- rep(c(sum(dat$bin_2_rep2_count)), times = nrow(dat))
total_bin_3_rep2_count <- rep(c(sum(dat$bin_3_rep2_count)), times = nrow(dat))
total_mNG2_neg_rep2_count <- rep(c(sum(dat$mNG2_neg_rep2_count)), times = nrow(dat))
total_mNG2_pos_rep2_count <- rep(c(sum(dat$mNG2_pos_rep2_count)), times = nrow(dat))

total_bin_0_rep3_count <- rep(c(sum(dat$bin_0_rep3_count)), times = nrow(dat))
total_bin_1_rep3_count <- rep(c(sum(dat$bin_1_rep3_count)), times = nrow(dat))
total_bin_2_rep3_count <- rep(c(sum(dat$bin_2_rep3_count)), times = nrow(dat))
total_bin_3_rep3_count <- rep(c(sum(dat$bin_3_rep3_count)), times = nrow(dat))

# Calculate total frequency for each mutant for every replicate.
dat$total_rep1_freq <- total_rep1_count/(total_bin_0_rep1_count + total_bin_1_rep1_count + total_bin_2_rep1_count + total_bin_3_rep1_count + total_mNG2_neg_rep1_count + total_mNG2_pos_rep1_count)

dat$total_rep2_freq <- total_rep2_count/(total_bin_0_rep2_count + total_bin_1_rep2_count + total_bin_2_rep2_count + total_bin_3_rep2_count + total_mNG2_neg_rep2_count + total_mNG2_pos_rep2_count)

dat$total_rep3_freq <- total_rep3_count/(total_bin_0_rep3_count + total_bin_1_rep3_count + total_bin_2_rep3_count + total_bin_3_rep3_count)

# Calculate input frequency per mutant.
dat$log10_input_freq <- log10((dat$ipt_count)/total_input_count)

# Calculate average total frequency per mutant.
dat$log10_avg_total_freq <- log10(rowMeans(dat[, c('total_rep1_freq', 'total_rep2_freq', 'total_rep3_freq')]))

# Calculate the frequency of each mutant per bin per replicate.
dat$bin_0_rep1_freq <- (dat$bin_0_rep1_count)/total_bin_0_rep1_count
dat$bin_1_rep1_freq <- (dat$bin_1_rep1_count)/total_bin_1_rep1_count
dat$bin_2_rep1_freq <- (dat$bin_2_rep1_count)/total_bin_2_rep1_count
dat$bin_3_rep1_freq <- (dat$bin_3_rep1_count)/total_bin_3_rep1_count
dat$mNG2_neg_rep1_freq <- (dat$mNG2_neg_rep1_count)/total_mNG2_neg_rep1_count
dat$mNG2_pos_rep1_freq <- (dat$mNG2_pos_rep1_count)/total_mNG2_pos_rep1_count

dat$bin_0_rep2_freq <- (dat$bin_0_rep2_count)/total_bin_0_rep2_count
dat$bin_1_rep2_freq <- (dat$bin_1_rep2_count)/total_bin_1_rep2_count
dat$bin_2_rep2_freq <- (dat$bin_2_rep2_count)/total_bin_2_rep2_count
dat$bin_3_rep2_freq <- (dat$bin_3_rep2_count)/total_bin_3_rep2_count
dat$mNG2_neg_rep2_freq <- (dat$mNG2_neg_rep2_count)/total_mNG2_neg_rep2_count
dat$mNG2_pos_rep2_freq <- (dat$mNG2_pos_rep2_count)/total_mNG2_pos_rep2_count

dat$bin_0_rep3_freq <- (dat$bin_0_rep3_count)/total_bin_0_rep3_count
dat$bin_1_rep3_freq <- (dat$bin_1_rep3_count)/total_bin_1_rep3_count
dat$bin_2_rep3_freq <- (dat$bin_2_rep3_count)/total_bin_2_rep3_count
dat$bin_3_rep3_freq <- (dat$bin_3_rep3_count)/total_bin_3_rep3_count

dat$log10_bin_0_rep1_freq <- log10((dat$bin_0_rep1_count)/total_bin_0_rep1_count)
dat$log10_bin_1_rep1_freq <- log10((dat$bin_1_rep1_count)/total_bin_1_rep1_count)
dat$log10_bin_2_rep1_freq <- log10((dat$bin_2_rep1_count)/total_bin_2_rep1_count)
dat$log10_bin_3_rep1_freq <- log10((dat$bin_3_rep1_count)/total_bin_3_rep1_count)
dat$log10_mNG2_neg_rep1_freq <- log10((dat$mNG2_neg_rep1_count)/total_mNG2_neg_rep1_count)
dat$log10_mNG2_pos_rep1_freq <- log10((dat$mNG2_pos_rep1_count)/total_mNG2_pos_rep1_count)

dat$log10_bin_0_rep2_freq <- log10((dat$bin_0_rep2_count)/total_bin_0_rep2_count)
dat$log10_bin_1_rep2_freq <- log10((dat$bin_1_rep2_count)/total_bin_1_rep2_count)
dat$log10_bin_2_rep2_freq <- log10((dat$bin_2_rep2_count)/total_bin_2_rep2_count)
dat$log10_bin_3_rep2_freq <- log10((dat$bin_3_rep2_count)/total_bin_3_rep2_count)
dat$log10_mNG2_neg_rep2_freq <- log10((dat$mNG2_neg_rep2_count)/total_mNG2_neg_rep2_count)
dat$log10_mNG2_pos_rep2_freq <- log10((dat$mNG2_pos_rep2_count)/total_mNG2_pos_rep2_count)

dat$log10_bin_0_rep3_freq <- log10((dat$bin_0_rep3_count)/total_bin_0_rep3_count)
dat$log10_bin_1_rep3_freq <- log10((dat$bin_1_rep3_count)/total_bin_1_rep3_count)
dat$log10_bin_2_rep3_freq <- log10((dat$bin_2_rep3_count)/total_bin_2_rep3_count)
dat$log10_bin_3_rep3_freq <- log10((dat$bin_3_rep3_count)/total_bin_3_rep3_count)

# Calculate weighted score for each mutant. Normalize according to WT (library).
weighted_score_rep1 <- (dat$bin_0_rep1_freq * 0.25 + dat$bin_1_rep1_freq * 0.5 + dat$bin_2_rep1_freq * 0.75 + dat$bin_3_rep1_freq * 1) / (dat$bin_0_rep1_freq + dat$bin_1_rep1_freq + dat$bin_2_rep1_freq + dat$bin_3_rep1_freq)
dat$weighted_score_rep1 <- weighted_score_rep1
WT_library_score_rep1 <- rep(c(dat[dat$mut == "WT", ]$weighted_score_rep1), times = nrow(dat))
dat$norm_weighted_score_rep1 <- weighted_score_rep1/WT_library_score_rep1

weighted_score_rep2 <- (dat$bin_0_rep2_freq * 0.25 + dat$bin_1_rep2_freq * 0.5 + dat$bin_2_rep2_freq * 0.75 + dat$bin_3_rep2_freq * 1) / (dat$bin_0_rep2_freq + dat$bin_1_rep2_freq + dat$bin_2_rep2_freq + dat$bin_3_rep2_freq)
dat$weighted_score_rep2 <- weighted_score_rep2
WT_library_score_rep2 <- rep(c(dat[dat$mut == "WT", ]$weighted_score_rep2), times = nrow(dat))
dat$norm_weighted_score_rep2 <- weighted_score_rep2/WT_library_score_rep2

weighted_score_rep3 <- (dat$bin_0_rep3_freq * 0.25 + dat$bin_1_rep3_freq * 0.5 + dat$bin_2_rep3_freq * 0.75 + dat$bin_3_rep3_freq * 1) / (dat$bin_0_rep3_freq + dat$bin_1_rep3_freq + dat$bin_2_rep3_freq + dat$bin_3_rep3_freq)
dat$weighted_score_rep3 <- weighted_score_rep3
WT_library_score_rep3 <- rep(c(dat[dat$mut == "WT", ]$weighted_score_rep3), times = nrow(dat))
dat$norm_weighted_score_rep3 <- weighted_score_rep3/WT_library_score_rep3

# Calculate average weighted score per mutant.
dat$avg_weighted_score <- rowMeans(dat[, c('norm_weighted_score_rep1', 'norm_weighted_score_rep2', 'norm_weighted_score_rep3')])

# Calculate average total frequency per mutant.
dat$avg_total_freq <- rowMeans(dat[, c('total_rep1_freq', 'total_rep2_freq', 'total_rep3_freq')])
dat$log10_avg_total_freq <- log10(dat$avg_total_freq)

# Calculate mNG2_pos_freq/mNG2_neg_freq per mutant for each replicate.
dat$ratio_rep1 <- (dat$mNG2_pos_rep1_freq)/(dat$mNG2_neg_rep1_freq)
WT_library_ratio_rep1 <- rep(c(dat[dat$mut == "WT", ]$ratio_rep1), times = nrow(dat))
dat$norm_ratio_rep1 <- dat$ratio_rep1/WT_library_ratio_rep1

dat$ratio_rep2 <- (dat$mNG2_pos_rep2_freq)/(dat$mNG2_neg_rep2_freq)
WT_library_ratio_rep2 <- rep(c(dat[dat$mut == "WT", ]$ratio_rep2), times = nrow(dat))
dat$norm_ratio_rep2 <- dat$ratio_rep2/WT_library_ratio_rep2

dat$avg_norm_ratio <- rowMeans(dat[, c('norm_ratio_rep1', 'norm_ratio_rep2')])
dat$log10_avg_norm_ratio <- log10(dat$avg_norm_ratio)

# Calculate fusion score normalized to average weighted score for each mutant.
dat$norm_fusion <- log10(dat$avg_norm_ratio)/dat$avg_weighted_score

# Only retain mutants where the total frequency in each replicate is at least 10^-4.15 is at least -4.15. Cutoff frequency is based on corr_score_norm_ratio >= 0.70 (to two decimal places).
dat <- dat[(dat$total_rep1_freq >= 10^-4.15) & (dat$total_rep2_freq >= 10^-4.15) & (dat$total_rep3_freq >= 10^-4.15), ]

# Establish a linear regression model for fusion against average weighted score.
Model <- lm(dat$log10_avg_norm_ratio~dat$avg_weighted_score)

# Calculate residuals.
dat$residual <- Model$resid

# Plot and calculate correlation between input frequency and average frequency.
p_corr_input_avg_freq <- ggplot(dat, aes(x = log10_input_freq, y = log10_avg_total_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  xlab(expression(bold(log[bold("10")](bold("input frequency"))))) +
  ylab(expression(bold(log[bold("10")](bold("total frequency"))))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_input_avg_freq <- cor(dat$log10_input_freq, dat$log10_avg_total_freq, method = 'pearson')

ggsave("corr_input_avg_freq.pdf", plot = p_corr_input_avg_freq, width = 6, height = 6, units = "in", dpi = 1200)

# Plot and calculate correlation coefficients of each bin.
p_corr_bin_0_rep12 <- ggplot(dat, aes(x = log10_bin_0_rep1_freq, y = log10_bin_0_rep2_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 2)") + 
  ggtitle("Bin 0") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_0_rep12 <- cor(dat$log10_bin_0_rep1_freq, dat$log10_bin_0_rep2_freq, method = 'pearson')

ggsave("corr_bin_0_rep12.pdf", plot = p_corr_bin_0_rep12, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_0_rep13 <- ggplot(dat, aes(x = log10_bin_0_rep1_freq, y = log10_bin_0_rep3_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 3)") + 
  ggtitle("Bin 0") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_0_rep13 <- cor(dat$log10_bin_0_rep1_freq, dat$log10_bin_0_rep3_freq, method = 'pearson')

ggsave("corr_bin_0_rep13.pdf", plot = p_corr_bin_0_rep13, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_0_rep23 <- ggplot(dat, aes(x = log10_bin_0_rep2_freq, y = log10_bin_0_rep3_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 2)", y = "Frequency (replicate 3)") + 
  ggtitle("Bin 0") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_0_rep23 <- cor(dat$log10_bin_0_rep2_freq, dat$log10_bin_0_rep3_freq, method = 'pearson')

ggsave("corr_bin_0_rep23.pdf", plot = p_corr_bin_0_rep23, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_1_rep12 <- ggplot(dat, aes(x = log10_bin_1_rep1_freq, y = log10_bin_1_rep2_freq)) +
  geom_point(size = 0.75, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 2)") + 
  ggtitle("Bin 1") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_1_rep12 <- cor(dat$log10_bin_1_rep1_freq, dat$log10_bin_1_rep2_freq, method = 'pearson')

ggsave("corr_bin_1_rep12.pdf", plot = p_corr_bin_1_rep12, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_1_rep13 <- ggplot(dat, aes(x = log10_bin_1_rep1_freq, y = log10_bin_1_rep3_freq)) +
  geom_point(size = 0.75, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 3)") + 
  ggtitle("Bin 1") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_1_rep13 <- cor(dat$log10_bin_1_rep1_freq, dat$log10_bin_1_rep3_freq, method = 'pearson')

ggsave("corr_bin_1_rep13.pdf", plot = p_corr_bin_1_rep13, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_1_rep23 <- ggplot(dat, aes(x = log10_bin_1_rep2_freq, y = log10_bin_1_rep3_freq)) +
  geom_point(size = 0.75, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 2)", y = "Frequency (replicate 3)") + 
  ggtitle("Bin 1") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_1_rep23 <- cor(dat$log10_bin_1_rep2_freq, dat$log10_bin_1_rep3_freq, method = 'pearson')

ggsave("corr_bin_1_rep23.pdf", plot = p_corr_bin_1_rep23, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_2_rep12 <- ggplot(dat, aes(x = log10_bin_2_rep1_freq, y = log10_bin_2_rep2_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 2)") + 
  ggtitle("Bin 2") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_2_rep12 <- cor(dat$log10_bin_2_rep1_freq, dat$log10_bin_2_rep2_freq, method = 'pearson')

ggsave("corr_bin_2_rep12.pdf", plot = p_corr_bin_2_rep12, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_2_rep13 <- ggplot(dat, aes(x = log10_bin_2_rep1_freq, y = log10_bin_2_rep2_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Rrequency (replicate 1)", y = "Frequency (replicate 3)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_2_rep13 <- cor(dat$log10_bin_2_rep1_freq, dat$log10_bin_2_rep3_freq, method = 'pearson')

ggsave("corr_bin_2_rep13.pdf", plot = p_corr_bin_2_rep13, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_2_rep23 <- ggplot(dat, aes(x = log10_bin_2_rep2_freq, y = log10_bin_2_rep3_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 2)", y = "Frequency (replicate 3)") + 
  ggtitle("Bin 2") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_2_rep23 <- cor(dat$log10_bin_2_rep2_freq, dat$log10_bin_2_rep3_freq, method = 'pearson')

ggsave("corr_bin_2_rep23.pdf", plot = p_corr_bin_2_rep23, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_3_rep12 <- ggplot(dat, aes(x = log10_bin_3_rep1_freq, y = log10_bin_3_rep2_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 2)") + 
  ggtitle("Bin 3") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_3_rep12 <- cor(dat$log10_bin_3_rep1_freq, dat$log10_bin_3_rep2_freq, method = 'pearson')

ggsave("corr_bin_3_rep12.pdf", plot = p_corr_bin_3_rep12, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_3_rep13 <- ggplot(dat, aes(x = log10_bin_3_rep1_freq, y = log10_bin_3_rep3_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 3)") + 
  ggtitle("Bin 3") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_3_rep13 <- cor(dat$log10_bin_3_rep1_freq, dat$log10_bin_3_rep3_freq, method = 'pearson')

ggsave("corr_bin_3_rep13.pdf", plot = p_corr_bin_3_rep13, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_bin_3_rep23 <- ggplot(dat, aes(x = log10_bin_3_rep2_freq, y = log10_bin_3_rep3_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 2)", y = "Frequency (replicate 3)") + 
  ggtitle("Bin 3") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_bin_3_rep23 <- cor(dat$log10_bin_3_rep2_freq, dat$log10_bin_3_rep3_freq, method = 'pearson')

ggsave("corr_bin_3_rep23.pdf", plot = p_corr_bin_3_rep23, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_mNG2_neg_rep12 <- ggplot(dat, aes(x = log10_mNG2_neg_rep1_freq, y = log10_mNG2_neg_rep2_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 2)") + 
  ggtitle("No fusion") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_mNG2_neg_rep12 <- cor(dat$log10_mNG2_neg_rep1_freq, dat$log10_mNG2_neg_rep2_freq, method = 'pearson')

ggsave("corr_mNG2_neg_rep12.pdf", plot = p_corr_mNG2_neg_rep12, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_mNG2_pos_rep12 <- ggplot(dat, aes(x = log10_mNG2_pos_rep1_freq, y = log10_mNG2_pos_rep2_freq)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Frequency (replicate 1)", y = "Frequency (replicate 2)") + 
  ggtitle("Fusion") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_mNG2_pos_rep12 <- cor(dat$log10_mNG2_pos_rep1_freq, dat$log10_mNG2_pos_rep2_freq, method = 'pearson')

ggsave("corr_mNG2_pos_rep12.pdf", plot = p_corr_mNG2_pos_rep12, width = 6, height = 6, units = "in", dpi = 1200)

# Plot and calculate correlation coefficients of scores of each bin.
p_corr_norm_score_rep12 <- ggplot(dat, aes(x = norm_weighted_score_rep1, y = norm_weighted_score_rep2)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Score (replicate 1)", y = "Score (replicate 2)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_norm_score_rep12 <- cor(dat$norm_weighted_score_rep1, dat$norm_weighted_score_rep2, method = 'pearson')

ggsave("corr_norm_score_rep12.pdf", plot = p_corr_norm_score_rep12, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_norm_score_rep13 <- ggplot(dat, aes(x = norm_weighted_score_rep1, y = norm_weighted_score_rep3)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Score (replicate 1)", y = "Score (replicate 3)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_norm_score_rep13 <- cor(dat$norm_weighted_score_rep1, dat$norm_weighted_score_rep3, method = 'pearson')

ggsave("corr_norm_score_rep13.pdf", plot = p_corr_norm_score_rep13, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_norm_score_rep23 <- ggplot(dat, aes(x = norm_weighted_score_rep2, y = norm_weighted_score_rep3)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Score (replicate 2)", y = "Score (replicate 3)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_norm_score_rep23 <- cor(dat$norm_weighted_score_rep2, dat$norm_weighted_score_rep3, method = 'pearson')

ggsave("corr_norm_score_rep23.pdf", plot = p_corr_norm_score_rep23, width = 6, height = 6, units = "in", dpi = 1200)

p_corr_norm_fusion_rep12 <- ggplot(dat, aes(x = norm_ratio_rep1, y = norm_ratio_rep2)) +
  geom_point(size = 0.5, alpha = 0.75, color = "black") + 
  labs(x = "Ratio (replicate 1)", y = "Ratio (replicate 2)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_norm_fusion_rep12 <- cor(dat$norm_ratio_rep1, dat$norm_ratio_rep2, method = 'pearson')

ggsave("corr_norm_fusion_rep12.pdf", plot = p_corr_norm_fusion_rep12, width = 6, height = 6, units = "in", dpi = 1200)

# Plot correlation between average weighted score and normalized frequency.
p_corr_score_freq <- ggplot(dat, aes(x = avg_weighted_score, y = log10_avg_total_freq)) +
  geom_point(size = 0.5, alpha = 0.75, colour = "black") +
  xlab("Average weighted score") + 
  ylab(expression(bold(log[bold("10")](bold("input frequency"))))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_score_freq <- cor(dat$avg_weighted_score, dat$log10_avg_total_freq, method = 'pearson')

# Save 'p_corr_weighted_score' to file
ggsave("corr_score_freq.pdf", plot = p_corr_score_freq, width = 6, height = 6, units = "in", dpi = 1200)

write.table(dat, file='SARS2_S2_HR1_DMS_freq.tsv', quote = FALSE, sep = '\t', col.names = NA)

score_dat <- data.frame(dat$mut, dat$mut_class, dat$avg_weighted_score, dat$log10_avg_total_freq, dat$log10_avg_norm_ratio, dat$norm_fusion, dat$residual)

colnames(score_dat) <- c("mut", "mut_class", "avg_weighted_score", "log10_avg_total_freq", "log10_avg_norm_ratio", "norm_fusion", "residual")

write.table(score_dat, file='SARS2_S2_HR1_DMS_score.tsv', quote = FALSE, sep = '\t', col.names = NA)

# For labelling specific data points in a correlation plot of weighted scores across both replicates
labelling <- rep(c(""), times = nrow(score_dat))
score_dat$labelling <- labelling

score_dat$labelling[score_dat$mut == "K986P"] <- "K986P"
score_dat$labelling[score_dat$mut == "V987P"] <- "V987P"
score_dat$labelling[score_dat$mut == "A892P"] <- "A892P"
score_dat$labelling[score_dat$mut == "A899P"] <- "A899P"
score_dat$labelling[score_dat$mut == "A942P"] <- "A942P"
score_dat$labelling[score_dat$mut == "D994E"] <- "D994E"
score_dat$labelling[score_dat$mut == "D994Q"] <- "D994Q"
score_dat$labelling[score_dat$mut == "T961F"] <- "T961F"
score_dat$labelling[score_dat$mut == "Q1005R"] <- "Q1005R"
score_dat$labelling[score_dat$mut == "S943H"] <- "S943H"
score_dat$labelling[score_dat$mut == "A944S"] <- "A944S"
score_dat$labelling[score_dat$mut == "WT"] <- "WT"
score_dat$labelling[score_dat$mut == "D950N"] <- "D950N"
score_dat$labelling[score_dat$mut == "Q954H"] <- "Q954H"
score_dat$labelling[score_dat$mut == "N969K"] <- "N969K"
score_dat$labelling[score_dat$mut == "L981F"] <- "L981F"
score_dat$labelling[score_dat$mut == "S982A"] <- "S982A"
score_dat$labelling[score_dat$mut == "T1027I"] <- "T1027I"

coloring <- rep(c(""), times = nrow(score_dat))
score_dat$coloring <- coloring

score_dat$coloring[score_dat$mut == "K986P"] <- "violetred"
score_dat$coloring[score_dat$mut == "V987P"] <- "violetred"
score_dat$coloring[score_dat$mut == "A892P"] <- "violetred"
score_dat$coloring[score_dat$mut == "A899P"] <- "violetred"
score_dat$coloring[score_dat$mut == "A942P"] <- "violetred"
score_dat$coloring[score_dat$mut == "D994E"] <- "lightseagreen"
score_dat$coloring[score_dat$mut == "D994Q"] <- "lightseagreen"
score_dat$coloring[score_dat$mut == "T961F"] <- "lightseagreen"
score_dat$coloring[score_dat$mut == "Q1005R"] <- "lightseagreen"
score_dat$coloring[score_dat$mut == "S943H"] <- "slateblue"
score_dat$coloring[score_dat$mut == "A944S"] <- "slateblue"
score_dat$coloring[score_dat$mut == "WT"] <- "darkred"
score_dat$coloring[score_dat$mut == "D950N"] <- "lightsalmon1"
score_dat$coloring[score_dat$mut == "Q954H"] <- "lightsalmon1"
score_dat$coloring[score_dat$mut == "N969K"] <- "lightsalmon1"
score_dat$coloring[score_dat$mut == "L981F"] <- "lightsalmon1"
score_dat$coloring[score_dat$mut == "S982A"] <- "lightsalmon1"
score_dat$coloring[score_dat$mut == "T1027I"] <- "lightsalmon1"

# For coloring each specific data points in a correlation plot of weighted scores across both replicates
K986P <- score_dat[score_dat$mut == "K986P", ]
V987P <- score_dat[score_dat$mut == "V987P", ]
A892P <- score_dat[score_dat$mut == "A892P", ]
A899P <- score_dat[score_dat$mut == "A899P", ]
A942P <- score_dat[score_dat$mut == "A942P", ]
D994E <- score_dat[score_dat$mut == "D994E", ]
D994Q <- score_dat[score_dat$mut == "D994Q", ]
T961F <- score_dat[score_dat$mut == "T961F", ]
Q1005R <- score_dat[score_dat$mut == "Q1005R", ]
S943H <- score_dat[score_dat$mut == "S943H", ]
A944S <- score_dat[score_dat$mut == "A944S", ]
WT <- score_dat[dat$mut == "WT", ]
D950N <- score_dat[dat$mut == "D950N", ]
Q954H <- score_dat[dat$mut == "Q954H", ]
N969K <- score_dat[dat$mut == "N969K", ]
L981F <- score_dat[dat$mut == "L981F", ]
S982A <- score_dat[dat$mut == "S982A", ]
T1027I <- score_dat[dat$mut == "T1027I", ]

mut_labels <- rbind(K986P, V987P, A892P, A899P, A942P, D994E, D994Q, T961F, Q1005R, S943H, A944S, WT, D950N, Q954H, N969K, L981F, S982A, T1027I)

stopMut <- subset(score_dat, mut_class == "nonsense")
missenseMut <- subset(score_dat, mut_class == "missense")
silentMut <- subset(score_dat, (mut_class == "silent" | mut_class == "WT"))

# Plot ratio vs score, and calculate correlation coefficient between ratio and score.
p_corr_score_ratio <- ggplot() +
  geom_point(stopMut, mapping = aes(x = avg_weighted_score, y = log10_avg_norm_ratio, color = log10_avg_total_freq, size = log10_avg_total_freq), alpha = 0.5) +
  scale_color_gradient2(midpoint = -3.5, low = "khaki", mid = "khaki2",
                        high = "khaki4", space ="Lab") +
  new_scale_color() +
  geom_point(silentMut, mapping = aes(x = avg_weighted_score, y = log10_avg_norm_ratio, color = log10_avg_total_freq, size = log10_avg_total_freq), alpha = 0.5) +
  scale_color_gradient2(midpoint = -3.5, low = "springgreen", mid = "springgreen2",
                        high = "springgreen4", space ="Lab") +
  geom_point(data = mut_labels, mapping = aes(x = avg_weighted_score, y = log10_avg_norm_ratio), color = mut_labels$coloring, alpha = 0.5) +
  new_scale_color() +
  geom_point(missenseMut, mapping = aes(x = avg_weighted_score, y = log10_avg_norm_ratio, color = log10_avg_total_freq, size = log10_avg_total_freq), alpha = 0.5) +
  scale_color_gradient2(midpoint = -2.5, low = "deepskyblue", mid = "deepskyblue2",
                        high = "deepskyblue4", space ="Lab") +
  geom_point(data = mut_labels, mapping = aes(x = avg_weighted_score, y = log10_avg_norm_ratio), color = mut_labels$coloring, alpha = 1) +
  geom_label_repel(data = score_dat, mapping = aes(x = avg_weighted_score, y = log10_avg_norm_ratio, label = labelling), fontface = "bold", color = score_dat$coloring, size = 3.5, box.padding = unit(1.2, 'lines'), point.padding = unit(0.7, 'lines'), segment.color = score_dat$coloring, segment.size = 0.6, arrow = arrow(length = unit(0.01, 'npc')), force = 2, max.iter = 1e9, max.overlaps = 1e5) +
  geom_smooth(data = dat, mapping = aes(x = avg_weighted_score, y = log10_avg_norm_ratio), method = "lm", color = "darkorange", size = 1, fill = "orange") +
  xlab("average expression score") + 
  ylab(expression( bold(paste("log"["10"], bgroup("(", frac( "positive","negative"), ")" ) )))) +
  scale_x_continuous(limits = c(0.5, 2), breaks = seq(0.5, 2, 0.5)) +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5)) +
  scale_size_continuous(limits = c(-4.5, -1), range = c(0.1, 2.5), breaks = c(-4, -3, -2, -1)) +
  guides(color = guide_colorbar(reverse = TRUE),
         size = guide_legend(reverse = FALSE)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10, face = "bold"))

corr_score_ratio_pearson <- cor(dat$avg_weighted_score, dat$log10_avg_norm_ratio, method = 'pearson')

# Save 'p_corr_score_ratio' to file
ggsave("corr_score_ratio.pdf", plot = p_corr_score_ratio, width = 9, height = 6, units = "in", dpi = 1200)

p_corr_score_norm_ratio <- ggplot(score_dat, aes(x = avg_weighted_score, y = norm_fusion, color = log10_avg_total_freq, size = log10_avg_total_freq)) +
  geom_point(alpha = 0.7) +
  geom_point(data = mut_labels, aes(x = avg_weighted_score, y = norm_fusion), color = mut_labels$coloring, alpha = 0.7) +
  geom_smooth(data = dat, mapping = aes(x = avg_weighted_score, y = norm_fusion), method = "lm", color = "darkgreen", size = 1, fill = "forestgreen") +
  geom_label_repel(data = score_dat, mapping = aes(x = avg_weighted_score, y = norm_fusion, label = labelling), fontface = "bold", color = score_dat$coloring, size = 6, box.padding = unit(1.2, 'lines'), point.padding = unit(0.4, 'lines'), segment.color = score_dat$coloring, segment.size = 0.6, arrow = arrow(length = unit(0.01, 'npc')), force = 6, max.iter = 1e9, max.overlaps = 5e5) +
  xlab("Average expression score") + 
  ylab("Normalized fusion score") +
  scale_color_gradientn(colours=c('grey85','grey45','grey5'), 
                        limits=c(-4.2,-0.8),
                        space = "Lab") +
  scale_x_continuous(limits = c(0.5, 2.0), breaks = seq(0.5, 2, 0.5)) +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5)) +
  scale_size_continuous(limits = c(-4.5, -1), range = c(0.05, 2.5), breaks = c(-4, -3, -2, -1)) +
  guides(color = guide_colorbar(reverse = TRUE),
         size = guide_legend(reverse = FALSE)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 16, face = "bold"))

ggsave("corr_score_norm_ratio.pdf", plot = p_corr_score_norm_ratio, width = 12, height = 9, units = "in", dpi = 1200)

corr_score_norm_ratio <- cor(dat$avg_weighted_score, dat$norm_fusion, method = 'pearson')

p_corr_score_residual <- ggplot(score_dat, aes(x = avg_weighted_score, y = residual, color = log10_avg_total_freq, size = log10_avg_total_freq)) +
  geom_point(alpha = 0.7) +
  geom_point(data = mut_labels, aes(x = avg_weighted_score, y = residual), color = mut_labels$coloring, alpha = 0.7) +
  geom_smooth(data = dat, mapping = aes(x = avg_weighted_score, y = residual), method = "lm", color = "darkgreen", size = 1, fill = "forestgreen") +
  geom_label_repel(data = score_dat, mapping = aes(x = avg_weighted_score, y = residual, label = labelling), fontface = "bold", color = score_dat$coloring, size = 6, box.padding = unit(1.2, 'lines'), point.padding = unit(0.4, 'lines'), segment.color = score_dat$coloring, segment.size = 0.6, arrow = arrow(length = unit(0.01, 'npc')), force = 6, max.iter = 1e9, max.overlaps = 5e5) +
  xlab("Average expression score") + 
  ylab("Residual") +
  scale_color_gradientn(colours=c('grey85','grey45','grey5'), 
                        limits=c(-4.2,-0.8),
                        space = "Lab") +
  scale_x_continuous(limits = c(0.5, 2.0), breaks = seq(0.5, 2, 0.5)) +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5)) +
  scale_size_continuous(limits = c(-4.5, -1), range = c(0.05, 2.5), breaks = c(-4, -3, -2, -1)) +
  guides(color = guide_colorbar(reverse = TRUE),
         size = guide_legend(reverse = FALSE)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        axis.ticks = element_line(colour = "black", size = 0.75, linetype = "solid"),
        axis.ticks.length = unit(0.18, "cm"),
        axis.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 16, face = "bold"))

ggsave("corr_score_residual.pdf", plot = p_corr_score_residual, width = 12, height = 9, units = "in", dpi = 1200)

corr_score_residual <- cor(dat$avg_weighted_score, dat$residual, method = 'pearson')