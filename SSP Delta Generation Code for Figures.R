#SSP DELTA DATA GENERATION
library(varbvs)
library(dplyr)
library(ggplot2)
library("OneR")
library("Hmisc")
library("ggbeeswarm")
library("cowplot")
library('tibble')

#SSP Model
SSP_Process <- function(sda = 1.861, gamma = 0.018, c = 100, runs = 100000, cut = 1000){
  resultdf <- data.frame()
  
  
  for(i in 1:runs){
    dt <- 0.10                                                                                  #Difference in time
    t <- seq(from = 0, to = 1000, by = 0.01)                                                    #Set up up time sequence
    sd_of_t <- sda-gamma*t                                                                      #Calculate shrinking of selective attention
    sd_of_t <- replace(sd_of_t, sd_of_t<0, 0)                                                   #Set spread of attention to 0 to avoid negative sd
    pouter<-1-pnorm(1.5, mean = 0, sd = sd_of_t)                                                #Calculate spread for outer flankers
    pinner<-pnorm(1.5, mean = 0, sd = sd_of_t) - pnorm(.5, mean = 0, sd = sd_of_t)              #Calculate spread for inner flankers
    ptarget <- pnorm(.5, mean = 0, sd = sd_of_t)- pnorm(-.5, mean = 0, sd = sd_of_t)            #Calculate spread for target
    
    bias <- 1                                                                                   #Bias for Congruent Condition
    v_t <- bias*2*pouter+bias*2*pinner+ptarget + c*sqrt(dt)*rnorm(length(t), mean = 0, sd = 1)  #Slope at each time point
    totalC <-cumsum(v_t)                                                                        #overall time-course
    CCrossC <- min(which(totalC > cut))/1000                                                        #Check Correct boundary
    ICrossC <- min(which(totalC < -cut))/1000                                                       #check Incorrect boundary
    
    
    bias <- -1                                                                                   #Bias for Congruent Condition
    v_tIN <- bias*2*pouter+bias*2*pinner+ptarget +c*sqrt(dt)*rnorm(length(t), mean = 0, sd = 1)  #Slope at each time point
    totalIN <-cumsum(v_tIN)                                                                      #overall time-course
    CCrossIN <- min(which(totalIN > cut))/1000                                                       #Check Correct boundary
    ICrossIN <- min(which(totalIN < -cut))/1000                                                      #check Incorrect boundary
    
    
    bias <- 0                                                                                   #Bias for Congruent Condition
    v_tN <- bias*2*pouter+bias*2*pinner+ptarget + c*sqrt(dt)*rnorm(length(t), mean = 0, sd = 1) #Slope at each time point
    totalN <-cumsum(v_tN)                                                                       #overall time-course
    CCrossN <-min(which(totalN > cut))/1000                                                         #Check Correct boundary
    ICrossN <- min(which(totalN < -cut))/1000                                                       #check Incorrect boundary
    
    
    resultdf <- rbind(resultdf,(c(CCrossC, CCrossN, CCrossIN)))
  }
  colnames(resultdf) <- c('Congruent','Neutral','Incongruent')
  bill <- resultdf %>% filter( (Congruent != Inf) &(Neutral != Inf)&(Incongruent != Inf))%>%
    transmute(Congruent = mean(Congruent),
              Neutral = mean(Neutral),
              Incongruent = mean(Incongruent))
  
  
  
  return(resultdf)
}
bS<-SSP_Process(c=100, cut = 2000)
ba<-SSP_Process(c=100, cut = 4000)

df<-ba
dfc <- df$Congruent/100
dfn <- df$Neutral/100
dfi <- df$Incongruent/100

dfdp <- mkt3(dfc,dfn, dfi)

out <- tapply(dfdp$obs, list(dfdp$gr), bin, nbins = 10, method = "content")
out1 <- sapply(split(dfdp$obs[dfdp$gr=="Group1"], out$Group1), "mean")
out2 <- sapply(split(dfdp$obs[dfdp$gr=="Group2"], out$Group2), "mean")
out3 <- sapply(split(dfdp$obs[dfdp$gr=="Group3"], out$Group3), "mean")
NCmeana <- rowMeans(cbind(out2,out1))
INmeana <- rowMeans(cbind(out3,out2))
DPmeana <- rowMeans(cbind(out3,out1))
dfNCa <- mkt1(out2 - out1)
dfINa <- mkt1(out3 - out2)
dfDPa <- mkt1(out3 - out1)
dfNCa$q <- NCmeana
dfINa$q <- INmeana
dfDPa$q <- DPmeana
colors <-c("Facilitation" = "Blue", "Inhibition" = "Red", "Standard" = "Black")


sspA <- ggplot(dfNCa, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  labs(color = "Legend")+
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(aes(y= dfINa$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfINa$obs))+
  geom_line(aes(y= dfDPa$obs, color = 'Standard') )+
  geom_point(aes(y= dfDPa$obs))+
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20)) +
  scale_y_continuous(breaks=seq(0,0.2,0.02), labels = NULL) +
  scale_x_continuous(labels = NULL) +
  coord_cartesian(ylim = c(0, 0.14)) +
  theme( axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 11), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = colors)+
  
  labs(title = "SSP Accuracy Focus", x = "t [s]", y = "\u0394 RT")


df<-bS
dfc <- df$Congruent/100
dfn <- df$Neutral/100
dfi <- df$Incongruent/100

dfdp <- mkt3(dfc,dfn, dfi)

out <- tapply(dfdp$obs, list(dfdp$gr), bin, nbins = 10, method = "content")
out1 <- sapply(split(dfdp$obs[dfdp$gr=="Group1"], out$Group1), "mean")
out2 <- sapply(split(dfdp$obs[dfdp$gr=="Group2"], out$Group2), "mean")
out3 <- sapply(split(dfdp$obs[dfdp$gr=="Group3"], out$Group3), "mean")
NCmeans <- rowMeans(cbind(out2,out1))
INmeans <- rowMeans(cbind(out3,out2))
DPmeans <- rowMeans(cbind(out3,out1))
dfNCs <- mkt1(out2 - out1)
dfINs <- mkt1(out3 - out2)
dfDPs <- mkt1(out3 - out1)
dfNCs$q <- NCmeans
dfINs$q <- INmeans
dfDPs$q <- DPmeans

colors <-c("Facilitation" = "Blue", "Inhibition" = "Red", "Standard" = "Black")

sspS <- ggplot(dfNCs, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  labs(color = "Legend")+
  theme_bw() +
  theme(legend.position = "none") +
  
  geom_line(aes(y= dfINs$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfINs$obs))+
  geom_line(aes(y= dfDPs$obs, color = 'Standard') )+
  geom_point(aes(y= dfDPs$obs))+
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20)) +
  scale_y_continuous(breaks=seq(0,0.2,0.02), labels = NULL) +
  scale_x_continuous(labels = NULL) +
  coord_cartesian(ylim = c(0, 0.14)) +
  theme( axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 11), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = colors)+
  
  labs(title = "SSP Speed Focus", x = "t [s]", y = "\u0394 RT")

