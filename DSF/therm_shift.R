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

df <- read.csv(file = 'therm_shift.csv', header = TRUE)
df2 <- melt(df, id = "Temperature")

therm_shift <- ggplot(df2, aes(x = Temperature, y = value, color = variable)) +
  geom_line(size = 1.5) +
  theme_cowplot(12) +
  xlab("Temperature (°C)") +
  ylab("-d(RFU)/dT") +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  expand_limits(x = c(5, 100)) + 
  scale_colour_discrete(labels = c("2P", "2PQ")) +
  geom_vline(xintercept = 46.5, linetype="dotted", color = "grey50", size=1.5) +
  geom_vline(xintercept = 62, linetype="dotted", color = "darkblue", size=1.5) +
  theme(plot.title=element_text(size=22+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
        axis.text=element_text(size=26,face="bold",colour = 'black'),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
        axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
        axis.title.x=element_text(size=26,face="bold"),
        axis.title.y=element_text(size=26,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        legend.title=element_blank(),
        legend.key.width = unit(1.5, "line"))

ggsave("therm_shift.pdf", plot = therm_shift, width = 10, height = 6)