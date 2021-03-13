##Plotting and tables for the MLE confidence intervals for the best-fitting model SSE transition rate parameters

library(hisse)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)

## Relaxed Standard
Relaxed.Standard.MuH.CID8.ci <- as.data.frame(Relaxed.Standard.MuH.CID8.AdaptiveSampling$ci)
Relaxed.Standard.MuH.CID8.ci <- tibble::rownames_to_column(Relaxed.Standard.MuH.CID8.ci, "Region")

## a little help from dplyr to keep only the columns with underscores (transition rate estimates) and delete zero columns
Relaxed.Standard.MuH.CID8.ci.Rates <- Relaxed.Standard.MuH.CID8.ci %>% select(contains("_"))
Relaxed.Standard.MuH.CID8.ci.Rates <- Relaxed.Standard.MuH.CID8.ci.Rates[, colSums(Relaxed.Standard.MuH.CID8.ci.Rates != 0) > 0]

## 1
Relaxed.Standard.MuH.CID8.ci.q00.q10 <- Relaxed.Standard.MuH.CID8.ci.Rates %>% select(contains("00") & contains("10"))
## q00 > q10
Relaxed.Standard.MuH.CID8.ci.Rates.00.10 <- Relaxed.Standard.MuH.CID8.ci.q00.q10 %>% select(contains("q00"))
## q10 > q00
Relaxed.Standard.MuH.CID8.ci.Rates.10.00 <- Relaxed.Standard.MuH.CID8.ci.q00.q10 %>% select(contains("q10"))

## 2
Relaxed.Standard.MuH.CID8.ci.Rates.q01q00 <- Relaxed.Standard.MuH.CID8.ci.Rates %>% select(contains("01") & contains("00"))
## q00 > q01
Relaxed.Standard.MuH.CID8.ci.Rates.00.01 <- Relaxed.Standard.MuH.CID8.ci.Rates.q01q00 %>% select(contains("q00"))
## q01 > q00
Relaxed.Standard.MuH.CID8.ci.Rates.01.00 <- Relaxed.Standard.MuH.CID8.ci.Rates.q01q00 %>% select(contains("q01"))

## 3
Relaxed.Standard.MuH.CID8.ci.Rates.q01q11 <- Relaxed.Standard.MuH.CID8.ci.Rates %>% select(contains("01") & contains("11"))
## q01 > q11
Relaxed.Standard.MuH.CID8.ci.Rates.01.11 <- Relaxed.Standard.MuH.CID8.ci.Rates.q01q11 %>% select(contains("q01"))
## q11 > q01
Relaxed.Standard.MuH.CID8.ci.Rates.11.01 <- Relaxed.Standard.MuH.CID8.ci.Rates.q01q11 %>% select(contains("q11"))

## 4
Relaxed.Standard.MuH.CID8.ci.Rates.q10q11 <- Relaxed.Standard.MuH.CID8.ci.Rates %>% select(contains("10") & contains("11"))
## q10 > q11
Relaxed.Standard.MuH.CID8.ci.Rates.10.11 <- Relaxed.Standard.MuH.CID8.ci.Rates.q10q11 %>% select(contains("q10"))
## q11 > q10
Relaxed.Standard.MuH.CID8.ci.Rates.11.10 <- Relaxed.Standard.MuH.CID8.ci.Rates.q10q11 %>% select(contains("q11"))

## Collate all into a new dataframe for plotting

RelaxedStandard.qrates.adaptivesampling.ci.q <- cbind(Relaxed.Standard.MuH.CID8.ci.Rates.00.10, 
                                                 Relaxed.Standard.MuH.CID8.ci.Rates.10.00, 
                                                 Relaxed.Standard.MuH.CID8.ci.Rates.00.01, 
                                                 Relaxed.Standard.MuH.CID8.ci.Rates.01.00, 
                                                 Relaxed.Standard.MuH.CID8.ci.Rates.01.11, 
                                                 Relaxed.Standard.MuH.CID8.ci.Rates.11.01, 
                                                 Relaxed.Standard.MuH.CID8.ci.Rates.10.11, 
                                                 Relaxed.Standard.MuH.CID8.ci.Rates.11.10)


RelaxedStandard.qrates.adaptivesampling.ci.q <- select(RelaxedStandard.qrates.adaptivesampling.ci.q,contains("A"))

RelaxedStandard.qrates.adaptivesampling.ci.q <-as.data.frame(RelaxedStandard.qrates.adaptivesampling.ci.q, row.names = c("0%",
                                                                                                               "25%",
                                                                                                               "50%",
                                                                                                               "75%",
                                                                                                               "100%"))

RelaxedStandard.qrates.adaptivesampling.ci.q <- tibble::rownames_to_column(RelaxedStandard.qrates.adaptivesampling.ci.q, "Region")
#write.csv(RelaxedStandard.qrates.adaptivesampling.ci.q,"MuH_CID8_RelaxedStandard.qrates.adaptivesampling.ci.q.csv")
RelaxedStandard.qrates.adaptivesampling.ci.q.Melted <- melt((RelaxedStandard.qrates.adaptivesampling.ci.q))
RelaxedStandard.qrates.adaptivesampling.ci.q.Melted$Region <- factor(RelaxedStandard.qrates.adaptivesampling.ci.q.Melted$Region, levels = c("0%", "25%", "50%","75%","100%"))

Relaxed.Standard.level.order <- factor(RelaxedStandard.qrates.adaptivesampling.ci.q.Melted$variable, level = c('q00A_01A',
                                                                                                     'q11A_10A',
                                                                                                     'q10A_11A',
                                                                                                     'q00A_10A',
                                                                                                     'q01A_00A',
                                                                                                     'q01A_11A',
                                                                                                     'q11A_01A',
                                                                                                     'q10A_00A'))

RelaxedStandard.qrates.adaptivesampling.ci.q.Melted.Jitter <- ggplot(RelaxedStandard.qrates.adaptivesampling.ci.q.Melted, aes(x=Relaxed.Standard.level.order, y = value))+ geom_jitter(aes(color = Region), 
                                                                                                                                                                        position = position_jitter(0), size = 3) +
  xlab("") + ylab("Transition Rate")+ggtitle("Relaxed Standard Adaptive Sampling (MuH CID8)")+ scale_color_manual(values=c( "#EFF3FF", "#BDD7E7" ,"#6BAED6", "#3182BD", "#08519C"))+theme_classic()+coord_flip()+ scale_fill_discrete(breaks=c("0%","25%","50%","75%","100%"))


## Strict Standard

Strict.Standard.MuH.CID5.ci <- as.data.frame(Strict.Standard.MuH.CID5.AdaptiveSampling$ci)
Strict.Standard.MuH.CID5.ci <- tibble::rownames_to_column(Strict.Standard.MuH.CID5.ci, "Region")

## State-Averaged Transition Rates 
#RelaxedStandardSolution <- Relaxed.Standard.MuH.CID8$solution

## Transpose the dataframe
#RelaxedStandardSolution.Transposed <- as.data.frame(t(as.matrix(RelaxedStandardSolution)))

## a little help from dplyr to keep only the columns with underscores (transition rate estimates) and delete zero columns
Strict.Standard.MuH.CID5.ci.Rates <- Strict.Standard.MuH.CID5.ci %>% select(contains("_"))
Strict.Standard.MuH.CID5.ci.Rates <- Strict.Standard.MuH.CID5.ci.Rates[, colSums(Strict.Standard.MuH.CID5.ci.Rates != 0) > 0]

## Categorize into the eight different transitions
## 00 > 10 DONE
## 00 > 01 DONE
## 10 > 00 DONE
## 01 > 00 DONE
## 01 > 11 DONE
## 11 > 01 DONE
## 11 > 10
## 10 > 11

## 1
Strict.Standard.MuH.CID5.ci.q00.q10 <- Strict.Standard.MuH.CID5.ci.Rates %>% select(contains("00") & contains("10"))
## q00 > q10
Strict.Standard.MuH.CID5.ci.Rates.00.10 <- Strict.Standard.MuH.CID5.ci.q00.q10 %>% select(contains("q00"))
## q10 > q00
Strict.Standard.MuH.CID5.ci.Rates.10.00 <- Strict.Standard.MuH.CID5.ci.q00.q10 %>% select(contains("q10"))

## 2
Strict.Standard.MuH.CID5.ci.Rates.q01q00 <- Strict.Standard.MuH.CID5.ci.Rates %>% select(contains("01") & contains("00"))
## q00 > q01
Strict.Standard.MuH.CID5.ci.Rates.00.01 <- Strict.Standard.MuH.CID5.ci.Rates.q01q00 %>% select(contains("q00"))
## q01 > q00
Strict.Standard.MuH.CID5.ci.Rates.01.00 <- Strict.Standard.MuH.CID5.ci.Rates.q01q00 %>% select(contains("q01"))

## 3
Strict.Standard.MuH.CID5.ci.Rates.q01q11 <- Strict.Standard.MuH.CID5.ci.Rates %>% select(contains("01") & contains("11"))
## q01 > q11
Strict.Standard.MuH.CID5.ci.Rates.01.11 <- Strict.Standard.MuH.CID5.ci.Rates.q01q11 %>% select(contains("q01"))
## q11 > q01
Strict.Standard.MuH.CID5.ci.Rates.11.01 <- Strict.Standard.MuH.CID5.ci.Rates.q01q11 %>% select(contains("q11"))

## 4
Strict.Standard.MuH.CID5.ci.Rates.q10q11 <- Strict.Standard.MuH.CID5.ci.Rates %>% select(contains("10") & contains("11"))
## q10 > q11
Strict.Standard.MuH.CID5.ci.Rates.10.11 <- Strict.Standard.MuH.CID5.ci.Rates.q10q11 %>% select(contains("q10"))
## q11 > q10
Strict.Standard.MuH.CID5.ci.Rates.11.10 <- Strict.Standard.MuH.CID5.ci.Rates.q10q11 %>% select(contains("q11"))

## Collate all into a new dataframe for plotting
StrictStandard.qrates.adaptivesampling.ci <- list(Strict.Standard.MuH.CID5.ci.Rates.00.10, 
                                                  Strict.Standard.MuH.CID5.ci.Rates.10.00, 
                                                  Strict.Standard.MuH.CID5.ci.Rates.00.01, 
                                                  Strict.Standard.MuH.CID5.ci.Rates.01.00, 
                                                  Strict.Standard.MuH.CID5.ci.Rates.01.11, 
                                                  Strict.Standard.MuH.CID5.ci.Rates.11.01, 
                                                  Strict.Standard.MuH.CID5.ci.Rates.10.11, 
                                                  Strict.Standard.MuH.CID5.ci.Rates.11.10) %>%
  setNames(c("00.10",
             "10.00",
             "00.01",
             "01.00",
             "01.11",
             "11.01",
             "10.11",
             "11.10"))

StrictStandard.qrates.adaptivesampling.ci.q <- cbind(Strict.Standard.MuH.CID5.ci.Rates.00.10, 
                                                     Strict.Standard.MuH.CID5.ci.Rates.10.00, 
                                                     Strict.Standard.MuH.CID5.ci.Rates.00.01, 
                                                     Strict.Standard.MuH.CID5.ci.Rates.01.00, 
                                                     Strict.Standard.MuH.CID5.ci.Rates.01.11, 
                                                     Strict.Standard.MuH.CID5.ci.Rates.11.01, 
                                                     Strict.Standard.MuH.CID5.ci.Rates.10.11, 
                                                     Strict.Standard.MuH.CID5.ci.Rates.11.10)


StrictStandard.qrates.adaptivesampling.ci.q <- select(StrictStandard.qrates.adaptivesampling.ci.q,contains("A"))

StrictStandard.qrates.adaptivesampling.ci.q <-as.data.frame(StrictStandard.qrates.adaptivesampling.ci.q, row.names = c("0%",
                                                                                                                       "25%",
                                                                                                                       "50%",
                                                                                                                       "75%",
                                                                                                                       "100%"))

StrictStandard.qrates.adaptivesampling.ci.q <- tibble::rownames_to_column(StrictStandard.qrates.adaptivesampling.ci.q, "Region")
#write.csv(StrictStandard.qrates.adaptivesampling.ci.q,"MuhCID5.StrictStandard.qrates.adaptivesampling.ci.q.csv")

StrictStandard.qrates.adaptivesampling.ci.q.Melted <- melt((StrictStandard.qrates.adaptivesampling.ci.q))
StrictStandard.qrates.adaptivesampling.ci.q.Melted$Region <- factor(StrictStandard.qrates.adaptivesampling.ci.q.Melted$Region, levels = c("0%", "25%", "50%","75%","100%"))
Strict.Standard.level.order <- factor(StrictStandard.qrates.adaptivesampling.ci.q.Melted$variable, level = c('q00A_01A',
                                                                                                       'q11A_10A',
                                                                                                       'q10A_11A',
                                                                                                       'q00A_10A',
                                                                                                       'q01A_00A',
                                                                                                       'q01A_11A',
                                                                                                       'q11A_01A',
                                                                                                       'q10A_00A'))
StrictStandard.qrates.adaptivesampling.ci.q.Melted.Jitter <- ggplot(StrictStandard.qrates.adaptivesampling.ci.q.Melted, aes(x = Strict.Standard.level.order,y=value))+ geom_jitter(aes(color = Region), 
                                                                                                                                                                                  position = position_jitter(0), size = 3) +
  xlab("") + ylab("Transition Rate")+ggtitle("Strict Standard Adaptive Sampling (MuH CID5)")+ scale_color_manual(values=c( "#EFF3FF", "#BDD7E7" ,"#6BAED6", "#3182BD", "#08519C"))+theme_classic()+coord_flip()+ scale_fill_discrete(breaks=c("0%","25%","50%","75%","100%"))

#brewer.pal(n = 5, name = "Blues")

#dev.size("in") #[1] 10.77778  5.87500

## Relaxed Full
Relaxed.Full.MuH.CID8.AdaptiveSampling <- Relaxed.MuH.CID8.CI
Relaxed.Full.MuH.CID8.ci <- as.data.frame(Relaxed.Full.MuH.CID8.AdaptiveSampling$ci)
Relaxed.Full.MuH.CID8.ci <- tibble::rownames_to_column(Relaxed.Full.MuH.CID8.ci, "Region")

## a little help from dplyr to keep only the columns with underscores (transition rate estimates) and delete zero columns
Relaxed.Full.MuH.CID8.ci.Rates <- Relaxed.Full.MuH.CID8.ci %>% select(contains("_"))
Relaxed.Full.MuH.CID8.ci.Rates <- Relaxed.Full.MuH.CID8.ci.Rates[, colSums(Relaxed.Full.MuH.CID8.ci.Rates != 0) > 0]

## 1
Relaxed.Full.MuH.CID8.ci.q00.q10 <- Relaxed.Full.MuH.CID8.ci.Rates %>% select(contains("00") & contains("10"))
## q00 > q10
Relaxed.Full.MuH.CID8.ci.Rates.00.10 <- Relaxed.Full.MuH.CID8.ci.q00.q10 %>% select(contains("q00"))
## q10 > q00
Relaxed.Full.MuH.CID8.ci.Rates.10.00 <- Relaxed.Full.MuH.CID8.ci.q00.q10 %>% select(contains("q10"))

## 2
Relaxed.Full.MuH.CID8.ci.Rates.q01q00 <- Relaxed.Full.MuH.CID8.ci.Rates %>% select(contains("01") & contains("00"))
## q00 > q01
Relaxed.Full.MuH.CID8.ci.Rates.00.01 <- Relaxed.Full.MuH.CID8.ci.Rates.q01q00 %>% select(contains("q00"))
## q01 > q00
Relaxed.Full.MuH.CID8.ci.Rates.01.00 <- Relaxed.Full.MuH.CID8.ci.Rates.q01q00 %>% select(contains("q01"))

## 3
Relaxed.Full.MuH.CID8.ci.Rates.q01q11 <- Relaxed.Full.MuH.CID8.ci.Rates %>% select(contains("01") & contains("11"))
## q01 > q11
Relaxed.Full.MuH.CID8.ci.Rates.01.11 <- Relaxed.Full.MuH.CID8.ci.Rates.q01q11 %>% select(contains("q01"))
## q11 > q01
Relaxed.Full.MuH.CID8.ci.Rates.11.01 <- Relaxed.Full.MuH.CID8.ci.Rates.q01q11 %>% select(contains("q11"))

## 4
Relaxed.Full.MuH.CID8.ci.Rates.q10q11 <- Relaxed.Full.MuH.CID8.ci.Rates %>% select(contains("10") & contains("11"))
## q10 > q11
Relaxed.Full.MuH.CID8.ci.Rates.10.11 <- Relaxed.Full.MuH.CID8.ci.Rates.q10q11 %>% select(contains("q10"))
## q11 > q10
Relaxed.Full.MuH.CID8.ci.Rates.11.10 <- Relaxed.Full.MuH.CID8.ci.Rates.q10q11 %>% select(contains("q11"))

## Collate all into a new dataframe for plotting

RelaxedFull.qrates.adaptivesampling.ci.q <- cbind(Relaxed.Full.MuH.CID8.ci.Rates.00.10, 
                                                     Relaxed.Full.MuH.CID8.ci.Rates.10.00, 
                                                     Relaxed.Full.MuH.CID8.ci.Rates.00.01, 
                                                     Relaxed.Full.MuH.CID8.ci.Rates.01.00, 
                                                     Relaxed.Full.MuH.CID8.ci.Rates.01.11, 
                                                     Relaxed.Full.MuH.CID8.ci.Rates.11.01, 
                                                     Relaxed.Full.MuH.CID8.ci.Rates.10.11, 
                                                     Relaxed.Full.MuH.CID8.ci.Rates.11.10)


RelaxedFull.qrates.adaptivesampling.ci.q <- select(RelaxedFull.qrates.adaptivesampling.ci.q,contains("A"))

RelaxedFull.qrates.adaptivesampling.ci.q <-as.data.frame(RelaxedFull.qrates.adaptivesampling.ci.q, row.names = c("0%",
                                                                                                                       "25%",
                                                                                                                       "50%",
                                                                                                                       "75%",
                                                                                                                       "100%"))

RelaxedFull.qrates.adaptivesampling.ci.q <- tibble::rownames_to_column(RelaxedFull.qrates.adaptivesampling.ci.q, "Region")
#write.csv(RelaxedFull.qrates.adaptivesampling.ci.q,"MuHCID8.RelaxedFull.qrates.adaptivesampling.ci.q.csv")

RelaxedFull.qrates.adaptivesampling.ci.q.Melted <- melt((RelaxedFull.qrates.adaptivesampling.ci.q))
RelaxedFull.qrates.adaptivesampling.ci.q.Melted$Region <- factor(RelaxedFull.qrates.adaptivesampling.ci.q.Melted$Region, levels = c("0%", "25%", "50%","75%","100%"))

Relaxed.Full.level.order <- factor(RelaxedFull.qrates.adaptivesampling.ci.q.Melted$variable, level = c('q00A_01A',
                                                                                                     'q11A_10A',
                                                                                                     'q10A_11A',
                                                                                                     'q00A_10A',
                                                                                                     'q01A_00A',
                                                                                                     'q01A_11A',
                                                                                                     'q11A_01A',
                                                                                                     'q10A_00A'))


RelaxedFull.qrates.adaptivesampling.ci.q.Melted.Jitter <- ggplot(RelaxedFull.qrates.adaptivesampling.ci.q.Melted, aes(x = Relaxed.Full.level.order, y =value))+ geom_jitter(aes(color = Region), 
                                                                                                                                                                                  position = position_jitter(0), size = 3) +
  xlab("") + ylab("Transition Rate")+ggtitle("Relaxed Full Adaptive Sampling (MuH CID8)")+ scale_color_manual(values=c( "#EFF3FF", "#BDD7E7" ,"#6BAED6", "#3182BD", "#08519C"))+theme_classic()+coord_flip()+ scale_fill_discrete(breaks=c("0%","25%","50%","75%","100%"))

## Strict Full
Strict.Full.MuH.CID7.ci <- as.data.frame(Strict.Full.MuH.CID7.AdaptiveSampling$ci)
Strict.Full.MuH.CID7.ci <- tibble::rownames_to_column(Strict.Full.MuH.CID7.ci, "Region")

## a little help from dplyr to keep only the columns with underscores (transition rate estimates) and delete zero columns
Strict.Full.MuH.CID7.ci.Rates <- Strict.Full.MuH.CID7.ci %>% select(contains("_"))
Strict.Full.MuH.CID7.ci.Rates <- Strict.Full.MuH.CID7.ci.Rates[, colSums(Strict.Full.MuH.CID7.ci.Rates != 0) > 0]

## 1
Strict.Full.MuH.CID7.ci.q00.q10 <- Strict.Full.MuH.CID7.ci.Rates %>% select(contains("00") & contains("10"))
## q00 > q10
Strict.Full.MuH.CID7.ci.Rates.00.10 <- Strict.Full.MuH.CID7.ci.q00.q10 %>% select(contains("q00"))
## q10 > q00
Strict.Full.MuH.CID7.ci.Rates.10.00 <- Strict.Full.MuH.CID7.ci.q00.q10 %>% select(contains("q10"))

## 2
Strict.Full.MuH.CID7.ci.Rates.q01q00 <- Strict.Full.MuH.CID7.ci.Rates %>% select(contains("01") & contains("00"))
## q00 > q01
Strict.Full.MuH.CID7.ci.Rates.00.01 <- Strict.Full.MuH.CID7.ci.Rates.q01q00 %>% select(contains("q00"))
## q01 > q00
Strict.Full.MuH.CID7.ci.Rates.01.00 <- Strict.Full.MuH.CID7.ci.Rates.q01q00 %>% select(contains("q01"))

## 3
Strict.Full.MuH.CID7.ci.Rates.q01q11 <- Strict.Full.MuH.CID7.ci.Rates %>% select(contains("01") & contains("11"))
## q01 > q11
Strict.Full.MuH.CID7.ci.Rates.01.11 <- Strict.Full.MuH.CID7.ci.Rates.q01q11 %>% select(contains("q01"))
## q11 > q01
Strict.Full.MuH.CID7.ci.Rates.11.01 <- Strict.Full.MuH.CID7.ci.Rates.q01q11 %>% select(contains("q11"))

## 4
Strict.Full.MuH.CID7.ci.Rates.q10q11 <- Strict.Full.MuH.CID7.ci.Rates %>% select(contains("10") & contains("11"))
## q10 > q11
Strict.Full.MuH.CID7.ci.Rates.10.11 <- Strict.Full.MuH.CID7.ci.Rates.q10q11 %>% select(contains("q10"))
## q11 > q10
Strict.Full.MuH.CID7.ci.Rates.11.10 <- Strict.Full.MuH.CID7.ci.Rates.q10q11 %>% select(contains("q11"))

## Collate all into a new dataframe for plotting

StrictFull.qrates.adaptivesampling.ci.q <- cbind(Strict.Full.MuH.CID7.ci.Rates.00.10, 
                                                  Strict.Full.MuH.CID7.ci.Rates.10.00, 
                                                  Strict.Full.MuH.CID7.ci.Rates.00.01, 
                                                  Strict.Full.MuH.CID7.ci.Rates.01.00, 
                                                  Strict.Full.MuH.CID7.ci.Rates.01.11, 
                                                  Strict.Full.MuH.CID7.ci.Rates.11.01, 
                                                  Strict.Full.MuH.CID7.ci.Rates.10.11, 
                                                  Strict.Full.MuH.CID7.ci.Rates.11.10)


StrictFull.qrates.adaptivesampling.ci.q <- select(StrictFull.qrates.adaptivesampling.ci.q,contains("A"))

StrictFull.qrates.adaptivesampling.ci.q <-as.data.frame(StrictFull.qrates.adaptivesampling.ci.q, row.names = c("0%",
                                                                                                                 "25%",
                                                                                                                 "50%",
                                                                                                                 "75%",
                                                                                                                 "100%"))

StrictFull.qrates.adaptivesampling.ci.q <- tibble::rownames_to_column(StrictFull.qrates.adaptivesampling.ci.q, "Region")
#write.csv(StrictFull.qrates.adaptivesampling.ci.q,"MuHCID8.StrictFull.qrates.adaptivesampling.ci.q.csv")

StrictFull.qrates.adaptivesampling.ci.q.Melted <- melt((StrictFull.qrates.adaptivesampling.ci.q))
StrictFull.qrates.adaptivesampling.ci.q.Melted$Region <- factor(StrictFull.qrates.adaptivesampling.ci.q.Melted$Region, levels = c("0%", "25%", "50%","75%","100%"))

Strict.Full.level.order <- factor(StrictFull.qrates.adaptivesampling.ci.q.Melted$variable, level = c('q00A_01A',
                                                                                                     'q11A_10A',
                                                                                                     'q10A_11A',
                                                                                                     'q00A_10A',
                                                                                                     'q01A_00A',
                                                                                                     'q01A_11A',
                                                                                                     'q11A_01A',
                                                                                                     'q10A_00A'))

StrictFull.qrates.adaptivesampling.ci.q.Melted.Jitter <- ggplot(StrictFull.qrates.adaptivesampling.ci.q.Melted, aes(x=Strict.Full.level.order, y = value))+ geom_jitter(aes(color = Region), 
                                                                                                                                                                        position = position_jitter(0), size = 3) +
  xlab("") + ylab("Transition Rate")+ggtitle("Strict Full Adaptive Sampling (MuH CID7)")+ scale_color_manual(values=c( "#EFF3FF", "#BDD7E7" ,"#6BAED6", "#3182BD", "#08519C"))+theme_classic()+coord_flip()+ scale_fill_discrete(breaks=c("0%","25%","50%","75%","100%"))

## Plot all
AdaptiveSamplingAll.Final <- ggarrange(
  StrictFull.qrates.adaptivesampling.ci.q.Melted.Jitter, RelaxedFull.qrates.adaptivesampling.ci.q.Melted.Jitter, StrictStandard.qrates.adaptivesampling.ci.q.Melted.Jitter,RelaxedStandard.qrates.adaptivesampling.ci.q.Melted.Jitter, labels = c("A", "C", "B", "D"),
  common.legend = TRUE, legend = "top")

#For mean + sd
# + stat_summary(aes(color = value), size = 0.3,fun.data="mean_sdl", alpha=0.7,color="gray70", fun.args = list(mult=1))

