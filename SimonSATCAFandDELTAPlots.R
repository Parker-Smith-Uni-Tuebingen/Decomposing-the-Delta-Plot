#Alt Analysis Template 5.0 prototype Simon

#Import in the appropriate packages
# packages <- c("tapply", "ggplot2", "plyr", "purrr", "tidyverse", "magrittr", "tidyr", "rstatix", "Rmisc" )

# install.packages(setdiff(packages, rownames(installed.packages())))

library(Rmisc)
library(ggplot2)            # plotting & data
library(tidyverse)          # data re-shaping
library(tidyr)  
library(magrittr)           # pipe operator
library(purrr)      
library(rstatix)
library(tidyverse)
library(purrr)
library(plyr)
library(dplyr)
library(ez)
library(schoRsch)
library(gridExtra)
library(Rmisc)
library(ggplot2)        # plotting & data
library(tidyverse)          # data re-shaping
library(tidyr)  
library(magrittr)       # pipe operator
library(purrr)      
library(rstatix)
library(tidyverse)
library(purrr)
library(plyr)
library(dplyr)
library(ez)
library(schoRsch)
library("OneR")
# library("rogme")
library("Hmisc")
# library("retimes")
library("ggbeeswarm")
library("cowplot")
library('tibble')

#Acquire Data
setwd("~/Phase II Data Collection/Phase II Simon Experiment Data Set 1")
filesAll <- list.files(pattern = "\\.csv$")
totaldat <- list()

for (i in seq_along(filesAll)) {
  totaldat[[i]] <- read.csv(file = filesAll[i])
}

#Set up containers
{
  LSimonC <- data.frame()
  SSimon <- data.frame()
  ASimon <- data.frame()
  BSimon <- data.frame()
  dat<- totaldat[[1]]
  dataAll <- data.frame()
  dataSAll <- data.frame()
  dataAAll <- data.frame()
  dataBAll <- data.frame()
  ANOVADatA <- data.frame()
  ANOVADatB <- data.frame()
  ANOVADatS <- data.frame()
  Aprop<- data.frame()
  Bprop<- data.frame()
  Sprop<- data.frame()
  SSimon <- data.frame()
  ASimon <- data.frame()
  BSimon <- data.frame()
  TRFdat<- data.frame()
  InterestingL<-data.frame()
  InterestingA<-data.frame()
  InterestingS<-data.frame()
  InterestingB<-data.frame()
}

#Analysis of Data
for (i in 1:length(totaldat)){
  dat<- totaldat[[i]]
  temp<- mutate(dat,
                Subject = i)
  
  #Filter out outliers and non-responses
  temp <- filter(temp, temp$Type == "Sym")
  temp <- filter(temp, temp$kr_trial_l.rt != "NA")
  temp <- filter(temp, temp$kr_trial_l.rt > 0.2)
  temp <- filter(temp, temp$kr_trial_l.rt < 2)
  
  #Select information items of interest
  temp <- select(temp, Type, Subject, ConditionSym, kr_trial_l.rt, kr_trial_l.corr, InstrCond)
  
  #Make container for all joint information
  dataAll<- rbind(dataAll, temp)
  
  #Create datasets for each speed focus
  dS <- temp
  dS <- filter(temp, temp$InstrCond == "speed")
  dataSAll<- rbind(dataSAll, dS)
  
  dA <- temp
  dA <- filter(temp, temp$InstrCond == "accuracy")
  dataAAll<- rbind(dataAAll, dA)
  
  dB <- temp
  dB <- filter(temp, temp$InstrCond == "both")
  dataBAll<- rbind(dataBAll, dB)
  
  #Calculate Percent correct and mean RT's for each dataset 
  reS.plyr<- dS  %>%
    group_by(ConditionSym) %>%
    mutate(PC = mean(kr_trial_l.corr))%>%
    ungroup() %>%
    filter(kr_trial_l.corr != 0)%>%
    group_by(ConditionSym) %>%
    mutate(RT_Mean = mean(kr_trial_l.rt))%>%
    ungroup()
  
  reA.plyr<- dA  %>%
    group_by(ConditionSym) %>%
    mutate(PC = mean(kr_trial_l.corr))%>%
    ungroup() %>%
    filter(kr_trial_l.corr != 0)%>%
    group_by(ConditionSym) %>%
    mutate(RT_Mean = mean(kr_trial_l.rt))%>%
    ungroup()
  
  reB.plyr<- dB  %>%
    group_by(ConditionSym) %>%
    mutate(PC = mean(kr_trial_l.corr))%>%
    ungroup() %>%
    filter(kr_trial_l.corr != 0)%>%
    group_by(ConditionSym) %>%
    mutate(RT_Mean = mean(kr_trial_l.rt))%>%
    ungroup()
  
  #Reorganize the datasets, then find mean RTs for each Congruency Condition
  reA.plyr <- reA.plyr %>%
    transmute(Subject = Subject,
              Conditions = ConditionSym,
              PC = PC,
              RT_Mean = RT_Mean,
              Instruction = InstrCond)%>%
    distinct()
  ANOVADatA <- rbind(ANOVADatA, reA.plyr)
  reA.plyr <- select(reA.plyr, Subject, Instruction, Conditions, RT_Mean)
  reA.plyr<-reA.plyr %>%
    mutate(Conditions = factor(Conditions, levels = unique(Conditions))) %>%
    spread(Conditions, RT_Mean)
  
  reS.plyr <- reS.plyr %>%
    transmute(Subject = Subject,
              Conditions = ConditionSym,
              PC = PC,
              RT_Mean = RT_Mean,
              Instruction = InstrCond)%>%
    distinct()
  ANOVADatS <- rbind(ANOVADatS, reS.plyr)
  reS.plyr <- select(reS.plyr, Subject, Instruction, Conditions, RT_Mean)
  reS.plyr<-reS.plyr %>%
    mutate(Conditions = factor(Conditions, levels = unique(Conditions))) %>%
    spread(Conditions, RT_Mean)
  
  reB.plyr <- reB.plyr %>%
    transmute(Subject = Subject,
              Conditions = ConditionSym,
              PC = PC,
              RT_Mean = RT_Mean,
              Instruction = InstrCond)%>%
    distinct()
  ANOVADatB <- rbind(ANOVADatB, reB.plyr)
  reB.plyr <- select(reB.plyr, Subject, Instruction, Conditions, RT_Mean)
  reB.plyr<-reB.plyr %>%
    mutate(Conditions = factor(Conditions, levels = unique(Conditions))) %>%
    spread(Conditions, RT_Mean)
  
  #Calculate R for each speed condition
  reA.ply <- transmute(reA.plyr,
                       Subject = Subject,
                       Instruction = Instruction,
                       R= (IN + C)/2-N,
                       R1= (N-C)/(IN-N))
  reS.ply <- transmute(reS.plyr,
                       Subject = Subject,
                       Instruction = Instruction,
                       R= (IN + C)/2-N,
                       R1= (N-C)/(IN-N))
  reB.ply <- transmute(reB.plyr,
                       Subject = Subject,
                       Instruction = Instruction,
                       R= (IN + C)/2-N,
                       R1= (N-C)/(IN-N))
  
  #Store the R's
  Aprop <- rbind(Aprop, reA.ply)
  Bprop <- rbind(Bprop, reB.ply)
  Sprop <- rbind(Sprop, reS.ply)
  
  #Store the resulting means and PC
  InterestingA <-rbind(InterestingA, reA.plyr)
  InterestingS <-rbind(InterestingS, reS.plyr)
  InterestingB <-rbind(InterestingB, reB.plyr)
}

#Merge together all ANOVA ready datasets
allDatResults <- data.frame()
allDatResults <-rbind(ANOVADatA,ANOVADatS)
allDatResults <-rbind(allDatResults,ANOVADatB)

#Anova results for RTs
res.anovaRT<- ezANOVA(allDatResults, dv = RT_Mean, wid = Subject, within = c(Instruction, Conditions), detailed = TRUE)
print(res.anovaRT)

#Anova results for PC
res.anovaPC<- ezANOVA(allDatResults, dv = PC, wid = Subject, within = c(Instruction, Conditions), detailed = TRUE)
print(res.anovaPC)


MofPer<- function(x){
  FirstP <- filter(x, (x$TP <= 0.2) & (x$TP > 0))
  SecondP <- filter(x, (x$TP <= 0.4) & (x$TP > 0.2))
  ThirdP <- filter(x, (x$TP <= 0.6) & (x$TP > 0.4))
  FourthP <- filter(x, (x$TP <= 0.8) & (x$TP > 0.6))
  FifthP <- filter(x, (x$TP <= 1) & (x$TP > 0.8))
  return(c(mean(FirstP$PC),mean(SecondP$PC),mean(ThirdP$PC),mean(FourthP$PC),mean(FifthP$PC)))
}


testSub1<- dataAAll

# testSub1 <- filter(testing, testing$Subject == 3)
# testSub1 <- filter(testSub1, testSub1$kr_trial_l.corr == 1)

testSub1C <- filter(testSub1, testSub1$ConditionSym == "C")
testSub1I <- filter(testSub1, testSub1$ConditionSym == "IN")
testSub1N <- filter(testSub1, testSub1$ConditionSym == "N")

TP <- seq(1,length(testSub1C$kr_trial_l.rt))
TestOrder <- order(testSub1C$kr_trial_l.rt)
TestReorderedC<- data.frame("PC" = testSub1C$kr_trial_l.corr[TestOrder], "RT" = testSub1C$kr_trial_l.rt[TestOrder])
# TestReorderedC$PC <-cumsum(TestReorderedC$PC/length(TestReorderedC$PC))
TestReorderedC$TP <- TP/length(TP)

TP <- seq(1,length(testSub1I$kr_trial_l.rt))
TestOrder <- order(testSub1I$kr_trial_l.rt)
TestReorderedI<- data.frame("PC" = testSub1I$kr_trial_l.corr[TestOrder], "RT" = testSub1I$kr_trial_l.rt[TestOrder])
# TestReorderedI$PC <-cumsum(TestReorderedI$PC/length(TestReorderedI$PC))
TestReorderedI$TP <- TP/length(TP)

TP <- seq(1,length(testSub1N$kr_trial_l.rt))
TestOrder <- order(testSub1N$kr_trial_l.rt)
TestReorderedN<- data.frame("PC" = testSub1N$kr_trial_l.corr[TestOrder], "RT" = testSub1N$kr_trial_l.rt[TestOrder])
# TestReorderedN$PC <-cumsum(TestReorderedN$PC/length(TestReorderedN$PC))
TestReorderedN$TP <- TP/length(TP)





#
stoCA<-MofPer(TestReorderedC)
stoIA<-MofPer(TestReorderedI)
stoNA<-MofPer(TestReorderedN)

stoFaA<- (stoCA- stoNA)*100
stoInA<- (stoNA- stoIA)*100
stoDPA<- (stoCA- stoIA)*100

testSub1<- dataBAll

# testSub1 <- filter(testing, testing$Subject == 3)
# testSub1 <- filter(testSub1, testSub1$kr_trial_l.corr == 1)

testSub1C <- filter(testSub1, testSub1$ConditionSym == "C")
testSub1I <- filter(testSub1, testSub1$ConditionSym == "IN")
testSub1N <- filter(testSub1, testSub1$ConditionSym == "N")

TP <- seq(1,length(testSub1C$kr_trial_l.rt))
TestOrder <- order(testSub1C$kr_trial_l.rt)
TestReorderedC<- data.frame("PC" = testSub1C$kr_trial_l.corr[TestOrder], "RT" = testSub1C$kr_trial_l.rt[TestOrder])
# TestReorderedC$PC <-cumsum(TestReorderedC$PC/length(TestReorderedC$PC))
TestReorderedC$TP <- TP/length(TP)

TP <- seq(1,length(testSub1I$kr_trial_l.rt))
TestOrder <- order(testSub1I$kr_trial_l.rt)
TestReorderedI<- data.frame("PC" = testSub1I$kr_trial_l.corr[TestOrder], "RT" = testSub1I$kr_trial_l.rt[TestOrder])
# TestReorderedI$PC <-cumsum(TestReorderedI$PC/length(TestReorderedI$PC))
TestReorderedI$TP <- TP/length(TP)

TP <- seq(1,length(testSub1N$kr_trial_l.rt))
TestOrder <- order(testSub1N$kr_trial_l.rt)
TestReorderedN<- data.frame("PC" = testSub1N$kr_trial_l.corr[TestOrder], "RT" = testSub1N$kr_trial_l.rt[TestOrder])
# TestReorderedN$PC <-cumsum(TestReorderedN$PC/length(TestReorderedN$PC))
TestReorderedN$TP <- TP/length(TP)




#
stoCB<-MofPer(TestReorderedC)
stoIB<-MofPer(TestReorderedI)
stoNB<-MofPer(TestReorderedN)

stoFaB<- (stoCB- stoNB)*100
stoInB<- (stoNB- stoIB)*100
stoDPB<- (stoCB- stoIB)*100
testSub1<- dataSAll

# testSub1 <- filter(testing, testing$Subject == 3)
# testSub1 <- filter(testSub1, testSub1$kr_trial_l.corr == 1)

testSub1C <- filter(testSub1, testSub1$ConditionSym == "C")
testSub1I <- filter(testSub1, testSub1$ConditionSym == "IN")
testSub1N <- filter(testSub1, testSub1$ConditionSym == "N")

TP <- seq(1,length(testSub1C$kr_trial_l.rt))
TestOrder <- order(testSub1C$kr_trial_l.rt)
TestReorderedC<- data.frame("PC" = testSub1C$kr_trial_l.corr[TestOrder], "RT" = testSub1C$kr_trial_l.rt[TestOrder])
# TestReorderedC$PC <-cumsum(TestReorderedC$PC/length(TestReorderedC$PC))
TestReorderedC$TP <- TP/length(TP)

TP <- seq(1,length(testSub1I$kr_trial_l.rt))
TestOrder <- order(testSub1I$kr_trial_l.rt)
TestReorderedI<- data.frame("PC" = testSub1I$kr_trial_l.corr[TestOrder], "RT" = testSub1I$kr_trial_l.rt[TestOrder])
# TestReorderedI$PC <-cumsum(TestReorderedI$PC/length(TestReorderedI$PC))
TestReorderedI$TP <- TP/length(TP)

TP <- seq(1,length(testSub1N$kr_trial_l.rt))
TestOrder <- order(testSub1N$kr_trial_l.rt)
TestReorderedN<- data.frame("PC" = testSub1N$kr_trial_l.corr[TestOrder], "RT" = testSub1N$kr_trial_l.rt[TestOrder])
# TestReorderedN$PC <-cumsum(TestReorderedN$PC/length(TestReorderedN$PC))
TestReorderedN$TP <- TP/length(TP)




#
stoCS<-MofPer(TestReorderedC)
stoIS<-MofPer(TestReorderedI)
stoNS<-MofPer(TestReorderedN)

stoFaS<- (stoCS- stoNS)*100
stoInS<- (stoNS- stoIS)*100
stoDPS<- (stoCS- stoIS)*100

x<-seq(1,5)
xl<-c("0-20", "20-40","40-60","60-80", "80-100")
df<- data.frame(stoCS, stoIS, stoNS, stoCA, stoIA, stoNA, stoCB, stoIB, stoNB, xl)
dfa <- data.frame(stoFaS, stoFaB, stoFaA, stoInS, stoInB, stoInA, xl,stoDPS,stoDPA,stoDPB)

colors <-  c( "Inhibition" = "red", "Facilitation" = "blue", "Standard" = "black")
# shapes <-  c( "Congruent" = 1, "Incongruent" = 2, "Neutral" = 0)
stuff <- c("1"="0-20", "2"="20-40","3"="40-60","4"="60-80", "5"="80-100")
stuff <- c("1"="0", "2"="20","3"="40","4"="60", "5"="80")

CAFSA <- ggplot(dfa, aes(x = xl, y = stoFaA, group = 1)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = stoFaA, group = 1,color = 'Facilitation'),  size = 0.5) +
  geom_point(aes( y= stoFaA)) +
  labs(color = "Legend")+
  theme_bw() + 
  theme(legend.position = "none") +
  geom_line(aes(y= stoInA, color = 'Inhibition') )+
  geom_point(aes(y= stoInA))+
  geom_line(aes(y= stoDPA, color = 'Standard'), linetype="longdash" )+
  geom_point(aes(y= stoDPA))+
  theme(axis.title = element_text(size = 12, face = "bold"), 
        axis.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  coord_cartesian(ylim = c(-5, 45)) +
  theme(legend.position = c(.7, .8))+
  # scale_y_continuous(breaks=seq(1,0.2,-0.2)) +
  scale_color_manual(values = c(colors))+
  # scale_x_discrete((labels=stuff))+
  # labs(x = "RT Bin (%)", y = "\u0394 PE (%)")
  labs(x = "RT Bin (%)", y = "\u0394 PE (%)")


CAFSB <- ggplot(dfa, aes(x = xl, y = stoFaB, group = 1)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = stoFaB,color = 'Facilitation', group = 1),  size = 0.5) +
  geom_point(aes( y= stoFaB)) +
  # labs(color = "Legend")+
  theme_bw() + 
  theme(legend.position = "none") +
  geom_line(aes(y= stoInB, color = 'Inhibition') )+
  geom_point(aes(y= stoInB))+
  geom_line(aes(y= stoDPB, color = 'Standard'), linetype="longdash" )+
  geom_point(aes(y= stoDPB))+
  theme(axis.title = element_text(size = 12, face = "bold"), 
        axis.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  coord_cartesian(ylim = c(-5, 45)) +
  # theme(legend.position = c(.7, .8))+
  scale_color_manual(values = c(colors))+
  # labs( x = "RT Bin (%)", y = "\u0394 PE (%)")
  labs(x = "RT Bin (%)",y="")


CAFSS <- ggplot(dfa, aes(x = xl, y = stoFaS, group = 1)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = stoFaS,color = 'Facilitation', group = 1),  size = 0.5) +
  geom_point(aes( y= stoFaS)) +
  # labs(color = "Legend")+
  theme_bw() + 
  theme(legend.position = "none") +
  geom_line(aes(y= stoInS, color = 'Inhibition') )+
  geom_point(aes(y= stoInS))+
  geom_line(aes(y= stoDPS, color = 'Standard'), linetype="longdash" )+
  geom_point(aes(y= stoDPS))+
  theme(axis.title = element_text(size = 12, face = "bold"), 
        axis.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  coord_cartesian(ylim = c(-5, 45)) +
  # theme(legend.position = c(.7, .8))+
  scale_color_manual(values = c(colors))+
  # labs( x = "RT Bin (%)", y = "\u0394 PE (%)")
  labs(x = "RT Bin (%)",y="")







library(gridExtra)

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

df<-dataAAll
df<- df %>% filter(df$kr_trial_l.corr == "1")
dfc <- df %>% filter(df$ConditionSym == "C")
dfn <- df %>% filter(df$ConditionSym == "N")
dfi <- df %>% filter(df$ConditionSym == "IN")

dfdp <- mkt3(dfc$kr_trial_l.rt,dfn$kr_trial_l.rt, dfi$kr_trial_l.rt)

out <- tapply(dfdp$obs, list(dfdp$gr), bin, nbins = 10, method = "content")
out1A <- sapply(split(dfdp$obs[dfdp$gr=="Group1"], out$Group1), "mean")
out2A <- sapply(split(dfdp$obs[dfdp$gr=="Group2"], out$Group2), "mean")
out3A <- sapply(split(dfdp$obs[dfdp$gr=="Group3"], out$Group3), "mean")
NCmean <- rowMeans(cbind(out2A,out1A))
INmean <- rowMeans(cbind(out3A,out2A))
DPmean <- rowMeans(cbind(out3A,out1A))
dfNCA <- mkt1(out2A - out1A)
dfINA <- mkt1(out3A - out2A)
dfDPA <-  mkt1(out3A - out1A)
dfNCA$q <- NCmean
dfINA$q <- INmean
dfDPA$q <- DPmean


dfdisA<- data.frame(out1A, out2A, out3A)

colors <-c("Facilitation" = "Blue", "Inhibition" = "Red", "Standard" = "Black")

fiSA <- ggplot(dfNCA, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  # scale_colour_manual(values=c("grey40", "orange1")) +
  # labs(color = "Legend")+
  theme_bw() + 
  theme(legend.position = "none") +
  geom_line(aes(y= dfINA$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfINA$obs))+
  geom_line(aes(y= dfDPA$obs, color = 'Standard'), linetype="longdash" )+
  geom_point(aes(y= dfDPA$obs))+
  theme(axis.title = element_text(size = 11, face = "bold"), 
        axis.text = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_y_continuous(breaks=seq(0,0.2,0.02)) +
  coord_cartesian(ylim = c(-0.01, 0.09)) +
  # theme(legend.position = c(.20, .8))+
  scale_color_manual(values = colors)+
  labs(title = "Simon Task - Accuracy", x = "t (s)", y = "\u0394 RT (s)")


df<-dataBAll
df<- df %>% filter(df$kr_trial_l.corr == "1")
dfc <- df %>% filter(df$ConditionSym == "C")
dfn <- df %>% filter(df$ConditionSym == "N")
dfi <- df %>% filter(df$ConditionSym == "IN")

dfdp <- mkt3(dfc$kr_trial_l.rt,dfn$kr_trial_l.rt, dfi$kr_trial_l.rt)

out <- tapply(dfdp$obs, list(dfdp$gr), bin, nbins = 10, method = "content")
out1B <- sapply(split(dfdp$obs[dfdp$gr=="Group1"], out$Group1), "mean")
out2B <- sapply(split(dfdp$obs[dfdp$gr=="Group2"], out$Group2), "mean")
out3B <- sapply(split(dfdp$obs[dfdp$gr=="Group3"], out$Group3), "mean")
NCmean <- rowMeans(cbind(out2B,out1B))
INmean <- rowMeans(cbind(out3B,out2B))
DPmean <- rowMeans(cbind(out3B,out1B))
dfNCB <- mkt1(out2B - out1B)
dfINB <- mkt1(out3B - out2B)
dfDPB <-  mkt1(out3B - out1B)
dfNCB$q <- NCmean
dfINB$q <- INmean
dfDPB$q <- DPmean

dfdisB<- data.frame(out1B, out2B, out3B)

colors <-c("Facilitation" = "Blue", "Inhibition" = "Red", "Standard" = "Black")
fiSB <- ggplot(dfNCB, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  # scale_colour_manual(values=c("grey40", "orange1")) +
  # labs(color = "Legend")+
  theme_bw() + 
  theme(legend.position = "none") +
  geom_line(aes(y= dfINB$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfINB$obs))+
  geom_line(aes(y= dfDPB$obs, color = 'Standard'), linetype="longdash" )+
  geom_point(aes(y= dfDPB$obs))+
  theme(axis.title = element_text(size = 11, face = "bold"), 
        axis.text = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank()) +
  scale_y_continuous(breaks=seq(0,0.2,0.02)) +
  coord_cartesian(ylim = c(-0.01, 0.09)) +
  # theme(legend.position = c(.25, .8))+
  scale_color_manual(values = colors)+
  labs(title = "Simon Task - Balanced", x = "t (s)",y="")
  # labs(title = "Simon Task - Balanced", x = "t (s)", y = "\u0394 RT (s)")





df<-dataSAll
df<- df %>% filter(df$kr_trial_l.corr == "1")
dfc <- df %>% filter(df$ConditionSym == "C")
dfn <- df %>% filter(df$ConditionSym == "N")
dfi <- df %>% filter(df$ConditionSym == "IN")

dfdp <- mkt3(dfc$kr_trial_l.rt,dfn$kr_trial_l.rt, dfi$kr_trial_l.rt)

out <- tapply(dfdp$obs, list(dfdp$gr), bin, nbins = 10, method = "content")
out1S <- sapply(split(dfdp$obs[dfdp$gr=="Group1"], out$Group1), "mean")
out2S <- sapply(split(dfdp$obs[dfdp$gr=="Group2"], out$Group2), "mean")
out3S <- sapply(split(dfdp$obs[dfdp$gr=="Group3"], out$Group3), "mean")
NCmean <- rowMeans(cbind(out2S,out1S))
INmean <- rowMeans(cbind(out3S,out2S))
DPmean <- rowMeans(cbind(out3S,out1S))
dfNCS <- mkt1(out2S - out1S)
dfINS <- mkt1(out3S - out2S)
dfDPS <-  mkt1(out3S - out1S)
dfNCS$q <- NCmean
dfINS$q <- INmean
dfDPS$q <- DPmean


dfdisS<- data.frame(out1S, out2S, out3S)


colors <-c("Facilitation" = "Blue", "Inhibition" = "Red", "Standard" = "Black")
fiSS <- ggplot(dfNCS, aes(x = q, y = obs)) + 
  geom_abline(intercept = 0, slope = 0, "colour"="grey40", linetype="dashed") +
  geom_line(aes(y = obs,color = 'Facilitation'),  size = 0.5, ) +
  geom_point(size = 2) +
  # scale_colour_manual(values=c("grey40", "orange1")) +
  # labs(color = "Legend")+
  theme_bw() + 
  theme(legend.position = "none") +
  geom_line(aes(y= dfINS$obs, color = 'Inhibition') )+
  geom_point(aes(y= dfINS$obs))+
  geom_line(aes(y= dfDPS$obs, color = 'Standard'), linetype="longdash" )+
  geom_point(aes(y= dfDPS$obs))+
  theme(axis.title = element_text(size = 11, face = "bold"), 
        axis.text = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank()) +
  scale_y_continuous(breaks=seq(0,0.2,0.02)) +
  coord_cartesian(ylim = c(-0.01, 0.09)) +
  # theme(legend.position = c(.25, .8))+
  scale_color_manual(values = colors)+
  labs(title = "Simon Task - Speed", x = "t (s)",y="" )
# labs(title = "Simon Task - Speed", x = "t (s)", y = "\u0394 RT (s)")



#Combine Graphs of CAF and DELTA
grid.arrange(fiSA, fiSB, fiSS, CAFSA, CAFSB, CAFSS, ncol=3, nrow =2)
grid.arrange(fiSA, fiSB, fiSS, CAFSA, CAFSB, CAFSS, ncol=3, nrow =2)
