library(ggplot2)

WT_100nM <- as.data.frame(read.delim("D1Results.txt", sep="\t"))
WT_50nM <- as.data.frame(read.delim("C1Results.txt", sep="\t"))
WT_25nM <- as.data.frame(read.delim("B1Results.txt", sep="\t"))
WT_13nM <- as.data.frame(read.delim("A1Results.txt", sep="\t"))

mut_100nM <- as.data.frame(read.delim("H1Results.txt", sep="\t"))
mut_50nM <- as.data.frame(read.delim("G1Results.txt", sep="\t"))
mut_25nM <- as.data.frame(read.delim("F1Results.txt", sep="\t"))
mut_13nM <- as.data.frame(read.delim("E1Results.txt", sep="\t"))

WT_100nM <- WT_100nM[-1:-7,]
WT_100nM[] <- lapply(WT_100nM, function(x) as.numeric(as.character(x)))
WT_100nM[,1] <- WT_100nM[,1] - WT_100nM[1,1]
rownames(WT_100nM) <- 1:nrow(WT_100nM)
colnames(WT_100nM) <- c("time", "raw_WT_100", "fit_WT_100")
WT_50nM <- WT_50nM[-1:-7,-1]
WT_50nM[] <- lapply(WT_50nM, function(x) as.numeric(as.character(x)))
rownames(WT_50nM) <- 1:nrow(WT_50nM)
colnames(WT_50nM) <- c("raw_WT_50", "fit_WT_50")
WT_25nM <- WT_25nM[-1:-7,-1]
WT_25nM[] <- lapply(WT_25nM, function(x) as.numeric(as.character(x)))
rownames(WT_25nM) <- 1:nrow(WT_25nM)
colnames(WT_25nM) <- c("raw_WT_25", "fit_WT_25")
WT_13nM <- WT_13nM[-1:-7,-1]
WT_13nM[] <- lapply(WT_13nM, function(x) as.numeric(as.character(x)))
rownames(WT_13nM) <- 1:nrow(WT_13nM)
colnames(WT_13nM) <- c("raw_WT_13", "fit_WT_13")

mut_100nM <- mut_100nM[-1:-7,-1]
mut_100nM[] <- lapply(mut_100nM, function(x) as.numeric(as.character(x)))
rownames(mut_100nM) <- 1:nrow(mut_100nM)
colnames(mut_100nM) <- c("raw_mut_100", "fit_mut_100")
mut_50nM <- mut_50nM[-1:-7,-1]
mut_50nM[] <- lapply(mut_50nM, function(x) as.numeric(as.character(x)))
rownames(mut_50nM) <- 1:nrow(mut_50nM)
colnames(mut_50nM) <- c("raw_mut_50", "fit_mut_50")
mut_25nM <- mut_25nM[-1:-7,-1]
mut_25nM[] <- lapply(mut_25nM, function(x) as.numeric(as.character(x)))
rownames(mut_25nM) <- 1:nrow(mut_25nM)
colnames(mut_25nM) <- c("raw_mut_25", "fit_mut_25")
mut_13nM <- mut_13nM[-1:-7,-1]
mut_13nM[] <- lapply(mut_13nM, function(x) as.numeric(as.character(x)))
rownames(mut_13nM) <- 1:nrow(mut_13nM)
colnames(mut_13nM) <- c("raw_mut_13", "fit_mut_13")

df <- cbind(WT_100nM, WT_50nM, WT_25nM, WT_13nM, mut_100nM, mut_50nM, mut_25nM, mut_13nM)

write.table(df, file = "WT_mut_BLI.txt", sep = "\t", row.names = TRUE, col.names = NA)

WT_BLI <- ggplot(df, aes(time)) +                    
  geom_line(aes(y=raw_WT_100), colour="grey70", size = 1.5) +
  geom_line(aes(y=raw_WT_50), colour="grey70", size = 1.5) +
  geom_line(aes(y=raw_WT_25), colour="grey70", size = 1.5) +
  geom_line(aes(y=raw_WT_13), colour="grey70", size = 1.5) +
  geom_line(aes(y=fit_WT_100), colour="orange", size = 1.5) +
  geom_line(aes(y=fit_WT_50), colour="orange", size = 1.5) +
  geom_line(aes(y=fit_WT_25), colour="orange", size = 1.5) +
  geom_line(aes(y=fit_WT_13), colour="orange", size = 1.5) +
  geom_vline(xintercept = 120, linetype="dotted", color = "darkgreen", size=1.5) +
  scale_x_continuous(breaks = seq(0, 420, by = 60)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.025)) +
  expand_limits(y = c(0, 0.1)) + 
  xlab("Time (s)") +
  ylab("Response (nm)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),    
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        strip.text.x = element_text(size = 10, face="bold", colour = "black"),
        strip.text.y = element_text(size = 10, face="bold", colour = "black"),
        axis.text.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=12, face="bold", colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title = element_text(size=22, face="bold", colour = "black"))

ggsave("WT_BLI.pdf", plot = WT_BLI, width = 4, height = 3)

mut_BLI <- ggplot(df, aes(time)) +                    
  geom_line(aes(y=raw_mut_100), colour="grey70", size = 1.5) +
  geom_line(aes(y=raw_mut_50), colour="grey70", size = 1.5) +
  geom_line(aes(y=raw_mut_25), colour="grey70", size = 1.5) +
  geom_line(aes(y=raw_mut_13), colour="grey70", size = 1.5) +
  geom_line(aes(y=fit_mut_100), colour="darkblue", size = 1.5) +
  geom_line(aes(y=fit_mut_50), colour="darkblue", size = 1.5) +
  geom_line(aes(y=fit_mut_25), colour="darkblue", size = 1.5) +
  geom_line(aes(y=fit_mut_13), colour="darkblue", size = 1.5) +
  geom_vline(xintercept = 120, linetype="dotted", color = "darkgreen", size=1.5) +
  scale_x_continuous(breaks = seq(0, 420, by = 60)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.025)) +
  expand_limits(y = c(0, 0.1)) + 
  xlab("Time (s)") +
  ylab("Response (nm)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),    
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        strip.text.x = element_text(size = 10, face="bold", colour = "black"),
        strip.text.y = element_text(size = 10, face="bold", colour = "black"),
        axis.text.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=12, face="bold", colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title = element_text(size=22, face="bold", colour = "black"))

ggsave("mut_BLI.pdf", plot = mut_BLI, width = 4, height = 3)
