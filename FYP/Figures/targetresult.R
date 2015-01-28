
library(ggplot2)

x = c(rep(0,26),rep(1,26))
y = c(rnorm(13,mean=28,sd=3.5),rnorm(13,mean=24,sd=3),
      rnorm(13,mean=40,sd=4),rnorm(13,mean=23,sd=5))
z = c(rep(0,13),rep(1,13),rep(0,13),rep(1,13))

d = data.frame(x=x,y=y,z=z)
d$x = factor(d$x,levels=c(0,1),labels=c("Contrast","Noise"))
d$z = factor(d$z,levels=c(0,1),labels=c("Focal","Distributed"))


ggplot(data=d,aes(x=z,y=y,color=x,fill=x)) +
  geom_point(stat="identity",size=3, 
             position=position_jitter(width=.05,height=0)) +
  geom_smooth(aes(group=x),method="lm") +
  xlab("Attention Condition") +
  ylab("Performance (% Correct / Total Presented)") +
  ggtitle("Performance for Unattended Targets by Task") +
  theme_classic()

