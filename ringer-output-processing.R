library(reshape)
library(stringr)
library(ggplot2)

#Read CSV
ringer = read.csv('~/Dropbox/Tests/ringer_tests/final_ringer.csv', header=FALSE)
ringer <- melt(ringer, id=c('V1','V2','V3','V4'))

#renames columns
colnames(ringer) <- c("residue", "map", "bond", "angle.offset", "angle.rel", "value") 


ringer$angle.rel <- 5*(as.numeric(str_replace_all(ringer$angle.rel, 'V', '')) - 5)
ringer$angle <- (ringer$angle.offset + ringer$angle.rel) %% 360

ringer <- ringer[order(ringer$residue, ringer$map, ringer$bond, ringer$angle),]

(p <- ggplot(subset(ringer, map=="2mFo-DFc" & bond=="chi1")) + 
   geom_line(aes(x=angle, y=value, color=residue)) +
   theme(legend.position="none") +
   scale_x_continuous(breaks=c(0,60,120,180,240,300,360)) + 
   facet_wrap( ~ residue)
)

