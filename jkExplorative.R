library(plotly)
library(ggplot2)
library(dplyr)

growData = read.csv("../Jain Krishna/Analysis/growAna.csv") %>% select(-one_of("X"))

summary(growData[,-1:-2])
pairs(growData[,-1:-2])
corgrow=cor(growData[,-1:-2])
View(corgrow)

ind=which(growData$av_core_path_mult<0.001)

plot(growData$lambda1[ind],growData$av_core_path_mult[ind])


(ggplot(growData, aes(x = lambda1)) +
  stat_density(position="identity",geom="line")) %>% ggplotly()

plot_ly(growData) %>% add_trace(x=~lambda1,type="histogram")

plot_ly(growData) %>% add_trace(x=~core_edges,type="histogram")


runData = read.csv("runAna2.csv")
runData[is.na(runData)] = 10000
runData = runData %>%
  inner_join((growData %>% select(c("network","nw_m","lambda1"))),
             by=c("parent"="network"))

run_ana = runData %>% filter(nw_m==0.2,threshold==75)

runsagg = runData %>% group_by(parent,threshold) %>%
  summarize_at((-1:-4),mean)
runsagg_sd = runData %>% group_by(parent,threshold) %>%
  summarize_at((-1:-4),sd)

pre_ana= runsagg %>% 
  inner_join(runsagg_sd,by=c("parent","threshold"),suffix = c(".mean", ".sd")) %>%
  inner_join(growData,by=c("parent"="network"))

ana=pre_ana %>%
  filter(threshold==50,nw_m==0.25)

plot_ly(ana) %>% add_trace(x=~core_nodes,type="histogram",histnorm = "probability")

(ggplot(ana, aes(x = lifetime.mean)) +
    stat_density(position="identity",geom="line")) %>% ggplotly()

cor(ana$lambda1,ana$lifetime.mean)

View(cor(ana[,unlist(lapply(ana, function(x){is.numeric(x) & !all(duplicated(x)[-1])}))],method="spearman"))

lm1 = lm(formula = lifetime ~ lambda1 + core_between_central_max + core_comple + core_edges, ana)
summary(lm1)

cor(ana$core_between_central_max,ana$core_edges)
avlog=log(runsagg$av_core_path_mult_complete)
cor(runsagg$av_core_cyclelen,runsagg$lambda1,method="spearman")


plot_ly(ana) %>%
  add_trace(x=~lambda1,y=~lifetime.mean,type="scatter",color=~parent,mode="markers",
            error_y = list(type = "data", array = ~lifetime.sd))

plot_ly(ana) %>%
  add_trace(x=~lifetime.mean,y=(ana$lifetime.sd/ana$lifetime.mean),type="scatter",mode="markers")

plot_ly(ana) %>%
  add_trace(x=~lifetime.mean,y=~recov_time.mean,type="scatter",mode="markers")

plot_ly(run_ana) %>%
  add_trace(x=~lambda1,y=~lifetime,type="scatter",mode="markers")


mean(runsagg_sd$lifetime)
runData %>% inner_join(growData,by=c("parent"="network")) %>% filter(nw_m==0.25) %>% 
  plot_ly() %>% add_trace(x=~lifetime,type="histogram")
(runsagg %>% filter(nw_m==0.25))[277,]
which((runsagg %>% filter(nw_m==0.25))$lifetime==max((runsagg %>% filter(nw_m==0.25))$lifetime))
runsagg[1921,"parent"]

mean(ana$lifetime.mean)
mean(ana$lifetime.sd)
sd(ana$lifetime.mean)



pairs(runsagg$lifetime)
