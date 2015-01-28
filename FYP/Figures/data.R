source("http://danbirman.com/birman_reset.R")

reps = 15

v1_c = rnorm(reps,1,.2)
v1_n = rnorm(reps,.1,.2)
v2_c = rnorm(reps,.9,.15)
v2_n = rnorm(reps,.1,.15)
v3_c = rnorm(reps,.7,.15)
v3_n = rnorm(reps,.23,.15)
v4_c = rnorm(reps,.4,.1)
v4_n = rnorm(reps,.3,.1)
ffa_c = rnorm(reps,.03,.1)
ffa_n = rnorm(reps,1.1,.1)
b = c(v1_c,v1_n,v2_c,v2_n,v3_c,v3_n,
          v4_c,v4_n,ffa_c,ffa_n)

cn = rep(c(rep(0,reps),rep(1,reps)),5)
greps = reps*2
group = c(rep(0,greps),rep(1,greps),rep(2,greps),rep(3,greps),
          rep(4,greps))

d = data.frame(task=cn,region=group,BOLD=b)
d$task = factor(d$task,levels=c(0,1),
                labels=c("Contrast","Noise"))
d$region = factor(d$region,levels=c(0,1,2,3,4),labels=
                    c("V1","V2","V3","V4","FFA"))

library(plyr)
dse = ddply(d,c("task","region"),summarise,m=mean(BOLD),se=sd(BOLD)/sqrt(5))

#library(ggplot2)
ggplot() +
#  geom_pointrange(data=dse,aes(x=region,y=m,color=task,
#                             ymin=m-se,ymax=m+se),size=1,shape="---") +
#  geom_smooth(data=d,aes(region,BOLD,color=task,fill=task,group=task),
#              method="lm",alpha=.05,formula=y~poly(x,4)) +
  geom_line(data=dse,aes(region,m,color=task,fill=task,group=task),size=1.5) +
  ylab("BOLD Signal Change due to Pedestal Increase (%)") +
  xlab("ROI") +
  theme_bw()

ggsave("C:\\Users\\Dan\\Documents\\Box Sync\\dan-hg-repo\\FYP\\Figures\\ROI_fig.pdf")
