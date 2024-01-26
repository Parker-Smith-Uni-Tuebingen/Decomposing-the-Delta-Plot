library(varbvs)
library(dplyr)
library(Rmisc)
library(ggplot2)            # plotting & data
library(tidyverse)          # data re-shaping
library("Hmisc")
library("ggbeeswarm")
library("cowplot")
library('tibble')


mkt1 <- function(x, name = c("obs", "gr")){
  gr <- factor(c(rep('Group1', length(x))))
  df <- data.frame(x, gr)
  names(df) <- name
  df <- tibble::as_tibble(df)
  df
}
mkt2 <- function(x, y, gr_names = "gr", obs_names = "obs",
                 group_labels = c("Group1", "Group2")){
  gr <- factor(c(rep(group_labels[1], length(x)),
                 rep(group_labels[2], length(y)))) # group labels
  obs <- c(x, y) # observations
  df <- tibble::tibble(gr, obs) # make tibble
  names(df) <- c(gr_names, obs_names)
  df
}
mkt3 <- function(x, y,z, gr_names = "gr", obs_names = "obs",
                 group_labels = c("Group1", "Group2", 'Group3')){
  gr <- factor(c(rep(group_labels[1], length(x)),
                 rep(group_labels[2], length(y)),
                 rep(group_labels[3], length(z)))) # group labels
  obs <- c(x, y,z) # observations
  df <- tibble::tibble(gr, obs) # make tibble
  names(df) <- c(gr_names, obs_names)
  df
}

DMCG<- function(A1 = 20, tau = 50, a = 2, k = 0.8, b = 50, sigma=4,  lambdad = 1/50, ad = 66,mu_c =  0.4, sides = F, bottom = F, runs = 100000){
  resultdf<-data.frame()
  for(i in 1:runs){
    # # Show both automatic and controlled signal 
    
    # DMC N
    # # Model parameters
    A<- 0;                                                                        # Set A=0 for Neutral trials
    dt<- 0.1; tmax<-1000;                                                          # Set change in time and max time
    
    # Compute time-dependent drift rate mu(t)
    t <- seq(from = dt, to = tmax, by = dt)                                       # Set time
    mu<-A*exp(-t/tau)*(exp(1)*t/(a-1)/tau)^(a-1)*((a-1)/t-1/tau)+ mu_c            # calculate slope at each time point
    # Simulate time-dependent Wiener process X(t) for a single trial
    dX<-mu*dt + sigma*sqrt(dt)*randn(1,length(t));                                # calculate change in x positon                                   
    XN<-cumsum(dX);                                                               # cumulative sum corresponds to X(t)
    
    #Calculate Cross of Boundaries - Neutral
    crossN <-min(which(XN > b))                                                   #Check Correct boundary
    IcrossN <- min(which(XN < -b))                                                #Check Inorrect boundary
    
    # DMC C
    # # Model parameters
    A<- A1;                                                                       # Set A=20 for Congruent trials
    
    # Compute time-dependent drift rate mu(t)
    t <- seq(from = dt, to = tmax, by = dt)                                       # Set time
    mu<-A*exp(-t/tau)*(exp(1)*t/(a-1)/tau)^(a-1)*((a-1)/t-1/tau) +  mu_c          # calculate slope at each time point
    # Simulate time-dependent Wiener process X(t) for a single trial
    dX<-mu*dt + sigma*sqrt(dt)*randn(1,length(t));                                # calculate change in x positon                           
    XC<-cumsum(dX);                                                               # cumulative sum corresponds to X(t)
    
    #Calculate Cross of Boundaries - Congruent
    crossC <-min(which(XC > b))                                                   #Check Correct boundary
    IcrossC <- min(which(XC < -b))                                                #Check Inorrect boundary 
    
    # DMC I
    # # Model parameters
    A<- -A1;                                                                      # Set A=-20 for Incongruent trials
    
    # Compute time-dependent drift rate mu(t)
    t <- seq(from = dt, to = tmax, by = dt)                                       # Set time
    mu<-A*exp(-t/tau)*(exp(1)*t/(a-1)/tau)^(a-1)*((a-1)/t-1/tau) + mu_c           # calculate slope at each time point
    # Simulate time-dependent Wiener process X(t) for a single trial
    dX<-mu*dt + sigma*sqrt(dt)*randn(1,length(t));                                # calculate change in x positon  
    XI<-cumsum(dX);                                                               # cumulative sum corresponds to X(t)
    
    crossI <-min(which(XI > b))                                                   #Check Correct boundary
    IcrossI <- min(which(XI < -b))                                                #Check Inorrect boundary 
    
    
    N_C <- crossN-crossC                                                          #Calculate Neutral - Congruent Differnece
    IN_N <- crossI - crossN                                                       #Calculate Incongruent - Neutral Differnece
    result <-N_C/IN_N                                                             #Calculate R
    
    resultdf <- rbind(resultdf,(c(crossC, crossN, crossI)))
    
    
    if (is.na(result)) {
      result<- 123456789
    }
    else if(result == Inf){                                                            #Check for any time courses that didn't cross the boundary
      result<- 123456789
    }
    else if(result == -Inf){                                                            #Check for any time courses that didn't cross the boundary
      result<- 123456789
    }
    else if((IcrossI <= crossI) | (IcrossN <= crossN) | (IcrossC <= crossC)){          #Check that a time course didn't hit the incorrect boundary first
      result<- 123456789
    }
  }
  
  colnames(resultdf) <- c('Congruent','Neutral','Incongruent')
  bill <- resultdf %>% filter( (Congruent != Inf) &(Neutral != Inf)&(Incongruent != Inf))%>%
    transmute(Congruent = mean(Congruent),
              Neutral = mean(Neutral),
              Incongruent = mean(Incongruent))
  
  return(resultdf)
  
}
#Flanker
aF<-DMCG(tau = 118.26, a = 2.15, A1=19.20 ,k = 0.8, b = 31.30, sigma=4,  lambdad = 1/50, ad = 66,mu_c =  0.69)
sF<-DMCG(tau = 118.26, a = 2.15, A1=19.20 ,k = 0.8, b = 71.30, sigma=4,  lambdad = 1/50, ad = 66,mu_c =  0.69)

df<-sF
dfc <- df$Congruent/10000
dfn <- df$Neutral/10000
dfi <- df$Incongruent/10000

dfdp <- mkt3(dfc,dfn, dfi)

out <- tapply(dfdp$obs, list(dfdp$gr), bin, nbins = 10, method = "content")
out1 <- sapply(split(dfdp$obs[dfdp$gr=="Group1"], out$Group1), "mean")
out2 <- sapply(split(dfdp$obs[dfdp$gr=="Group2"], out$Group2), "mean")
out3 <- sapply(split(dfdp$obs[dfdp$gr=="Group3"], out$Group3), "mean")
NCmean <- rowMeans(cbind(out2,out1))
INmean <- rowMeans(cbind(out3,out2))
DPmean <- rowMeans(cbind(out3,out1))
dfNC <- mkt1(out2 - out1)
dfIN <- mkt1(out3 - out2)
dfsF <- mkt1(out3 - out1)
dfNC$q <- NCmean
dfIN$q <- INmean
dfsF$q <- DPmean

colors <-c("Facilitation" = "Blue", "Inhibition" = "Red", "Standard" = "Black")
sF <- ggplot(dfNC, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  labs(color = "Legend")+
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(aes(y= dfIN$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfIN$obs))+
  geom_line(aes(y= dfsF$obs, color = 'Standard') )+
  geom_point(aes(y= dfsF$obs))+
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20)) +
  scale_y_continuous(breaks=seq(0,0.2,0.02), labels = NULL) +
  scale_x_continuous(labels = NULL) +
  coord_cartesian(ylim = c(0, 0.08)) +
  # theme(legend.position = c(.75, .2))+
  theme( legend.text = element_text(face = "bold", size = 10), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 11), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = colors)+
  labs( title = "DMC Speed Focus", x = "t [s]", y = "\u0394 RT")


df<-aF
dfc <- df$Congruent/10000
dfn <- df$Neutral/10000
dfi <- df$Incongruent/10000

dfdp <- mkt3(dfc,dfn, dfi)

out <- tapply(dfdp$obs, list(dfdp$gr), bin, nbins = 10, method = "content")
out1 <- sapply(split(dfdp$obs[dfdp$gr=="Group1"], out$Group1), "mean")
out2 <- sapply(split(dfdp$obs[dfdp$gr=="Group2"], out$Group2), "mean")
out3 <- sapply(split(dfdp$obs[dfdp$gr=="Group3"], out$Group3), "mean")
NCmean <- rowMeans(cbind(out2,out1))
INmean <- rowMeans(cbind(out3,out2))
DPmean <- rowMeans(cbind(out3,out1))
dfNC <- mkt1(out2 - out1)
dfIN <- mkt1(out3 - out2)
dfaF <- mkt1(out3 - out1)
dfNC$q <- NCmean
dfIN$q <- INmean
dfaF$q <- DPmean

aF <- ggplot(dfNC, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  labs(color = "Legend")+
  theme_bw() +
  geom_line(aes(y= dfIN$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfIN$obs))+
  geom_line(aes(y= dfaF$obs, color = 'Standard') )+
  geom_point(aes(y= dfaF$obs))+
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20)) +
  scale_y_continuous(breaks=seq(0,0.2,0.02), labels = NULL) +
  scale_x_continuous(labels = NULL) +
  coord_cartesian(ylim = c(0, 0.08)) +
  theme(legend.position = c(.2, .75))+
  theme( legend.text = element_text(face = "bold", size = 10), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 11), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = colors)+
  labs( title = "DMC Accuracy Focus", x = "t [s]", y = "\u0394 RT")


#Simon
aS <- DMCG(tau = 34.94, a = 2.80, A1=15.99 ,k = 0.8, b = 74.56, sigma=4,  lambdad = 1/50, ad = 66,mu_c =  0.69)
sS <- DMCG(tau = 34.94, a = 2.80, A1=15.99 ,k = 0.8, b = 34.56, sigma=4,  lambdad = 1/50, ad = 66,mu_c =  0.69)
df<-sS
dfc <- df$Congruent/10000
dfn <- df$Neutral/10000
dfi <- df$Incongruent/10000

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
dfsS <- mkt1(out3 - out1)
dfNCs$q <- NCmeans
dfINs$q <- INmeans
dfsS$q <- DPmeans


s <- ggplot(dfNCs, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  labs(color = "Legend")+
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(aes(y= dfINs$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfINs$obs))+
  geom_line(aes(y= dfsS$obs, color = 'Standard') )+
  geom_point(aes(y= dfsS$obs))+
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20)) +
  scale_y_continuous(breaks=seq(0,0.2,0.02), labels = NULL) +
  scale_x_continuous(labels = NULL) +
  coord_cartesian(ylim = c(0, 0.06)) +
  theme( legend.text = element_text( size = 10),
         axis.title.y = element_text(size = 14), 
         axis.title.x = element_text(size = 14),
         axis.text.x = element_text(size = 11), 
         axis.line = element_line(colour = "black"))+
  scale_color_manual(values = colors)+
  labs( title = "DMC Speed Focus", x = "t [s]", y = "\u0394 RT")


df<-aS
dfc <- df$Congruent/10000
dfn <- df$Neutral/10000
dfi <- df$Incongruent/10000

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
dfaS <- mkt1(out3 - out1)
dfNCa$q <- NCmeana
dfINa$q <- INmeana
dfaS$q <- DPmeana
a <- ggplot(dfNCa, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  labs(color = "Legend")+
  theme_bw() +
  geom_line(aes(y= dfINa$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfINa$obs))+
  geom_line(aes(y= dfaS$obs, color = 'Standard') )+
  geom_point(aes(y= dfaS$obs))+
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20)) +
  scale_y_continuous(breaks=seq(0,0.2,0.02), labels = NULL) +
  scale_x_continuous(labels = NULL) +
  coord_cartesian(ylim = c(0, 0.06)) +
  theme(legend.position = c(.8, .8))+
  theme( legend.text = element_text( size = 10),axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 11), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = colors)+
  labs( title = "DMC Accuracy Focus", x = "t [s]", y = "\u0394 RT")
