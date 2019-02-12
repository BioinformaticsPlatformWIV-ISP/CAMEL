#!/usr/bin/env Rscript
library("ggplot2")

# Parse options
options <- commandArgs(trailingOnly = TRUE)
print("Options:")
print(options)

# Plot the coverage
df = data.frame(index=rep(1, 3), cat=factor(c('fail', 'warn', 'good'), levels=c('good', 'warn', 'fail')), cov=c(10, 10, 30))
g <- ggplot(data=df, aes(x=index, y=cov, fill=cat)) +
  geom_bar(stat='identity', size=0.5, colour="black") +
  scale_fill_manual(values=c("#6EBE50", "#FFA839", "#FF8080")) +
  scale_y_continuous(breaks=c(0, 10, 20, 30, 40, 50), labels=c('0X', '10X', '20X', '30X', '40x', '>50X')) +
  coord_flip() +
  geom_hline(yintercept=as.double(options[2]), size=2) +
  ylab("Coverage") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  guides(fill=FALSE)
ggsave(options[1], g, width=5, height=0.70)
