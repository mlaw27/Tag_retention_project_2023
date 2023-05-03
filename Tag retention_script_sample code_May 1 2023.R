library(car)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(survival)
library(ranger)
library (ggplot2)
library(ggfortify)
library(ggpattern)
library(glmm)
library(lme4)
library(lmerTest)
library(nlme)
library(emmeans) 
library(broom)
library(metafor)
library(gridExtra)
library(forestplot)
library(ggforestplot)
library(readxl)
library(janitor)
library(survminer)
library(psych)


#####Part 1: Experiential Series: Tag retention and mortality in Atlantic salmon#####


setwd('')

Tag_data_all <-read.csv('TAGRETENTION.1.csv') %>% as.tibble() %>%  
  rename_all(tolower) %>% print()

####Experimental Mortality analyses#####

#simplify with just the core stuff that I need
#weight gain calc^ 
#(%) = 100[(Wt-Wi)/Wi]. 


Tag_mort <- Tag_data_all %>% select(survival..days.,survival, 
                                    tag.retention..days., dummy.type, 
                                    tag.weight..g., tag.length..mm.,
                                    initial.weight..g.,final.weight..g., 
                                    initial.length..mm., final.length..mm.,sex,
                                    rejection..y.n.,tank..) %>% 
  mutate(Perc_weight_change = 100*((final.weight..g. - initial.weight..g.)/initial.weight..g.)) %>% 
  rename(Tag_type = dummy.type) %>% 
  mutate(logWi = log10(initial.weight..g.)) %>% 
  mutate(logWf = log10(final.weight..g.)) %>% 
  mutate(logLi = log10(initial.length..mm.)) %>% 
  mutate(logLf = log10(final.length..mm.)) %>% 
  mutate(DeltaW = logWf - logWi) %>% 
  mutate(DeltaL = logLf - logLi) %>%
  rename(tank = tank..) %>% 
  print()



str(Tag_mort)
Tag_mort$sex<-as.factor(Tag_mort$sex)
Tag_mort$Tag_type <-as.factor(Tag_mort$Tag_type)

#using the package 'survival' to calculates if rate of failure changes with tag 
#type and other covaraites (ie weight change, k change, sex, tag properties
#and rejection)


#resetting dat for use in the survival package

Tag_mort_Surv <- Tag_mort %>% rename(time = survival..days., 
                                     status = survival, 
                                     rejection = rejection..y.n.) %>% 
  select(Tag_type, time, status, sex, rejection, 
         Perc_weight_change, tank) %>% 
  print()
#180 8


Tag_mort_Surv$status <-str_replace(Tag_mort_Surv$status, "Y","0") 
Tag_mort_Surv$status <-str_replace(Tag_mort_Surv$status, "N","1") 

Tag_mort_Surv$status

Tag_mort_Surv$status <-as.numeric(Tag_mort_Surv$status)
str(Tag_mort_Surv$status)
#write.csv(Tag_mort_Surv, "Tag retention data_cleaned up_April 24 2023.csv")

Mort_results <-with(Tag_mort_Surv, Surv(time, status)) %>% print() 
#making survival object ^


####Mortality survival model######

                                     
#http://www.sthda.com/english/wiki/cox-proportional-hazards-model

Mort_results_Cox_simple <- coxph(Surv(time, status)~Tag_type,
                                 data =Tag_mort_Surv) 

summary(Mort_results_Cox_simple)

# Call:
#   coxph(formula = Surv(time, status) ~ Tag_type, data = Tag_mort_Surv)
# 
# n= 180, number of events= 6 
# 
#             coef    exp(coef) se(coef)       z Pr(>|z|)
# Tag_typeV7 1.1428    3.1356   1.1547 0.990    0.322
# Tag_typeV8 0.6919    1.9976   1.2248 0.565    0.572
# 
#              exp(coef)  exp(-coef)lower .95 upper .95
# Tag_typeV7     3.136     0.3189    0.3262     30.15
# Tag_typeV8     1.998     0.5006    0.1811     22.03
# 
# Concordance= 0.619  (se = 0.101 )
# Likelihood ratio test= 1.14  on 2 df,   p=0.6
# Wald test            = 1.02  on 2 df,   p=0.6
# Score (logrank) test = 1.1  on 2 df,   p=0.6


#No effect of tag type on mortility of smolts


#summary for death by group is as follows: 

A<-Tag_mort_Surv %>% group_by(Tag_type) %>%
  summarise(no_rows = length(status)) %>% print()

#overall total counts
# # A tibble: 3 x 2
# Tag_type no_rows
# <fct>      <int>
#   1 SHAM        60
# 2 V7            59
# 3 V8            61

#TOTAL = 180

#counts of number of zeros (ie survivals)
B<- Tag_mort_Surv %>% 
  group_by(Tag_type) %>% 
  summarise_each(funs(sum(.==0))) %>% print()

# Tag_type       status    
# <fct>          <int>              
# 1 SHAM           59                      
# 2 V7             56                    
# 3 V8             59                   
# 
#TOTAL = 174
#seems alright as there were 6 morts during the exp

#percent survivals:

perc_surv <- left_join(A,B) %>% mutate(Perc = (status/no_rows)*100) %>% 
  select(Tag_type, Perc) %>% print()

# # A tibble: 3 x 2
# Tag_type  Perc
# <fct>    <dbl>
#   1 SHAM      98.3
# 2 V7        94.9
# 3 V8        96.7


####Mortality survival plot####

#using the survminer package to do so:
#https://rpkgs.datanovia.com/survminer/

library(survminer)

Fit <-survfit(Surv(time, status)~Tag_type,
              data =Tag_mort_Surv) 

A<-ggsurvplot(Fit,legend.labs = c('Sham', 
                                  'V7', 'V8'),
              conf.int = F, data = Tag_mort_Surv, risk.table = T,
              risk.table.height = .3, xlab = "Time in days", ylim = c(.75, 1),
              palette = c ('red3', 'grey', 'grey32' ), 
              ggtheme = theme_classic2(base_family = "Times", base_size=14),
              linetype = "strata")
A$plot <-A$plot +scale_linetype_manual(values = c ("dotted","solid","dashed"))
print(A)

#use manual export to save image



####Time to death#####

Tag_mort_Surv %>% group_by(Tag_type) %>% 
  summarise(n=n(),
            mean = mean(time),
            sd = sd(time),
            se=sd/sqrt(n))


# # A tibble: 3 x 5
# Tag_type     n  mean    sd    se
# <fct>    <int> <dbl> <dbl> <dbl>
# 1 SHAM        60  90.8  9.30  1.2 
# 2 V7          59  87.6 19.1   2.49
# 3 V8          61  89.2 15.4   1.97





####Growth and condition indices#####


######Log w vs log L 

names(Tag_mort)


LogWvL_model<-lmer(DeltaL~DeltaW+Tag_type+sex+(1|tank), data = Tag_mort)

plot(LogWvL_model)
qqnorm(resid(LogWvL_model))
qqline(resid(LogWvL_model))


#plots look good, residuals very scattered and almost perfect qq plot

summary(LogWvL_model)




# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: DeltaL ~ DeltaW + Tag_type + sex + (1 | tank)
#    Data: Tag_mort
# 
# REML criterion at convergence: -1055.8
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.68015 -0.69254 -0.05714  0.66365  2.52236 
# 
# Random effects:
#  Groups   Name        Variance  Std.Dev.
#  tank     (Intercept) 1.702e-06 0.001305
#  Residual             9.530e-05 0.009762
# Number of obs: 172, groups:  tank, 2
# 
# Fixed effects:
#               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)   0.023679   0.003971  68.874060   5.962 9.55e-08 ***
# DeltaW        0.211017   0.011435 166.481154  18.453  < 2e-16 ***
# Tag_typeV7   -0.003824   0.001905 166.024616  -2.007   0.0464 *  
# Tag_typeV8   -0.001002   0.001826 166.864189  -0.549   0.5840    
# sexM          0.002493   0.001538 166.999898   1.621   0.1069    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#            (Intr) DeltaW Tg_tV7 Tg_tV8
# DeltaW     -0.888                     
# Tag_typeV7 -0.451  0.260              
# Tag_typeV8 -0.373  0.169  0.507       
# sexM       -0.111 -0.135 -0.005 -0.027
# 




#Pairwise comparisons 


Log_pair = emmeans(LogWvL_model, ~ Tag_type)

pairs(Log_pair)

# contrast  estimate      SE  df t.ratio p.value
# SHAM - V7  0.00382 0.00191 166   2.007  0.1137
# SHAM - V8  0.00100 0.00184 167   0.545  0.8494
# V7 - V8   -0.00282 0.00186 167  -1.514  0.2869




####SGR


#first calculate IGR

#calc'd according to crane et al. 2019 : https://onlinelibrary.wiley.com/doi/full/10.1111/raq.12396
# g = loge(w2)-loge(w1)/ t2-t1


Tag_mort.1 <- Tag_mort %>% mutate(Inst_grow = (((log(final.weight..g.))-(log(initial.weight..g.)))/survival..days.)) %>% print()



#calulate SGR 

Tag_mort.2 <- Tag_mort.1 %>%mutate(SGR = 100*((exp(Inst_grow))-1)) %>% print()

Tag_mort.2$SGR

SGR_model<-lmer(SGR~Tag_type+sex+(1|tank), data = Tag_mort.2)

plot(SGR_model)

qqnorm(resid(SGR_model))
qqline(resid(SGR_model))


#Plots looks fine


summary(SGR_model)

# 
# Linear mixed model fit by REML ['lmerModLmerTest']
# Formula: SGR ~ Tag_type + sex + (1 | tank)
# Data: Tag_mort.2
# REML criterion at convergence: -109.8278
# Random effects:
#   Groups   Name        Std.Dev.
# tank     (Intercept) 0.03291 
# Residual             0.16573 
# Number of obs: 172, groups:  tank, 2
# Fixed Effects:
#   (Intercept)   Tag_typeV7   Tag_typeV8         sexM  
# 0.77505     -0.10892     -0.06883      0.04658  
# > plot(SGR_model)
# > summary(SGR_model)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: SGR ~ Tag_type + sex + (1 | tank)
#    Data: Tag_mort.2
# 
# REML criterion at convergence: -109.8
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.25169 -0.72199  0.03397  0.65050  2.50875 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  tank     (Intercept) 0.001083 0.03291 
#  Residual             0.027468 0.16573 
# Number of obs: 172, groups:  tank, 2
# 
# Fixed effects:
#              Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)   0.77505    0.03543   3.15954  21.878 0.000148 ***
# Tag_typeV7   -0.10892    0.03124 167.06378  -3.487 0.000624 ***
# Tag_typeV8   -0.06883    0.03057 167.35106  -2.251 0.025666 *  
# sexM          0.04658    0.02590 167.63028   1.799 0.073860 .  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#            (Intr) Tg_tV7 Tg_tV8
# Tag_typeV7 -0.434              
# Tag_typeV8 -0.429  0.486       
# sexM       -0.444  0.032 -0.005
# 





#Pairwise contrasts from the above model 

marginal = emmeans(SGR_model, ~ Tag_type)

pairs(marginal)

# 
# contrast  estimate     SE  df t.ratio p.value
# SHAM - V7   0.1089 0.0312 167   3.486  0.0018
# SHAM - V8   0.0688 0.0306 167   2.246  0.0665
# V7 - V8    -0.0401 0.0315 168  -1.273  0.4121
# 
# Results are averaged over the levels of: sex 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 
# 


#####mean values - growth


#SGR

Tag_mort.2 %>%
  group_by(Tag_type) %>%
  summarize(mean(SGR, na.rm = TRUE))

# # A tibble: 3 x 2
# Tag_type `mean(SGR, na.rm = TRUE)`
# <fct>                        <dbl>
#   1 SHAM                         0.802
# 2 V7                           0.694
# 3 V8                           0.736
# 



####growth bar plots

setwd('')
windowsFonts(Times=windowsFont("TT Times New Roman")) 





###SGR

Tag_mort.2 %>% group_by(Tag_type) %>% summarise (med = median(SGR))

ggplot(Tag_mort.2, aes(x= Tag_type, y=SGR, fill = Tag_type)) +
  geom_boxplot()+ 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c('SHAM' = "white", 'V7' = "grey", 'V8'= "grey32"))+
  theme(axis.text.x = element_text(size = 14, family= 'Times'),
        axis.text.y = element_text(size = 12, family= 'Times'))+
  theme(axis.title = element_text(size = 17, family= 'Times'))+
  labs(x ="Tagging Group", y = expression('Specific Growth Rate'~('%'~Delta~g~day^-1)))+
  scale_y_continuous(breaks=seq(0,.03,.001))+
  scale_x_discrete(labels= Labels_1)+
  theme(legend.position = "None" )+
  #annotate("text", x = .6, y = 2, label = "C", fontface = 'bold', family ='Times', size =8)+
  #annotate("text", x = 2.7, y = 2, label = "SGR", fontface = 'italic', family ='Times', size =8)+
  annotate("text", x = 1, y = 1.4, label = "A", fontface = 'bold', family ='Times', size =6)+
  annotate("text", x = 2, y = 1.4, label = "B", fontface = 'bold', family ='Times', size =6)+
  annotate("text", x = 3, y = 1.4, label = "AB", fontface = 'bold', family ='Times', size =6)+
  ylim(0, 2) 
ggsave('.tiff')



####Length weight relationship plot
Tag_mort.4<-Tag_mort.2 %>% rename(Treatment = Tag_type) %>%
  print()


Tag_mort.4$Treatment <-str_replace(Tag_mort.4$Treatment, "SHAM","Sham") 



ggplot(Tag_mort.4, aes(x=DeltaW, y=DeltaL, shape=Treatment, color=Treatment)) +
  geom_point()+
  geom_smooth(method = lm,se=F,aes(group=1),color='black')+ #for plotting just an overall regression
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c('Sham' = "#D81B60", 'V7' = "#1E88E5", 'V8'= "#FFC107"))+
  scale_size_manual(values=c(5,5,5))+
  theme(axis.text.x = element_text(size = 14, family= 'Times'),
        axis.text.y = element_text(size = 12, family= 'Times'))+
  theme(axis.title = element_text(size = 17, family= 'Times'))+
  xlab(bquote(Delta~Log[10]~'Body Mass'))+
  ylab(bquote(Delta~Log[10]~'Fork Length'))+
  theme(legend.position = c(0.2, 0.85))+theme(text=element_text(family="Times"))+
  theme(legend.text=element_text(size=11))


ggsave('.tiff')




#####Experimental Tag Retention Analyses######


#####Tag retention survival model######


setwd('')


#prep data frame and clean it up a bit   


Tagging_data<- read.csv("TAGRETENTION.1.csv")%>% as.tibble() %>% 
  select(TAG.ID,SEX, SURVIVAL..days.,TAG.RETENTION..days., DUMMY.TYPE, #taking what I need
         TAG.WEIGHT..g.,TAG.VOLUME..mm.,SUTURES.REMAINING, RESPONSE.TO.TAG,RESPONSE.TO.TAG.1,
         STATE.OF.TAG.REJECTION,REJECTION, REJECTION..Y.N.) %>% 
  rename(surv_time = SURVIVAL..days., tag_type = DUMMY.TYPE, #renaming the cols
         Tag_weight = TAG.WEIGHT..g., Tag_vol = TAG.VOLUME..mm.,
         sutures_left = SUTURES.REMAINING,Response = RESPONSE.TO.TAG,
         Reject_state = STATE.OF.TAG.REJECTION,rejection_progress =REJECTION, rejection_yanaw = REJECTION..Y.N.,
         Tag_reten_days = TAG.RETENTION..days., tag_response.1 = RESPONSE.TO.TAG.1) %>% 
  rename_all(tolower) %>%  print()

Tagging_data$sex <-as.factor(Tagging_data$sex)
Tagging_data$tag_type <-as.factor(Tagging_data$tag_type)
Tagging_data$tag.id <-as.factor(Tagging_data$tag.id)

str(Tagging_data)#180

#dropping controls and any fish that didn't have a tag retention 
#also adding in a col for status 

str(Tagging_data$tag_reten_days)

Tagging_data.1 <- Tagging_data %>% drop_na(tag_reten_days) %>% print()



Tagging_data.2<-Tagging_data.1%>% 
  mutate(status = case_when(tag_reten_days == 92~0,tag_reten_days < 92~1)) %>% 
  print()


#survival curve of the number of days that days were retained in the body

#resetting data for use in the survival package
Tag_retain_Surv <- Tagging_data.2 %>% rename(time = tag_reten_days) %>% 
  print()


Tag_retain_Surv$status <-as.numeric(Tag_retain_Surv$status)
str(Tag_retain_Surv$status)


Reten_results_Cox_simple <- coxph(Surv(time, status)~tag_type,
                                  data =Tag_retain_Surv) 

summary(Reten_results_Cox_simple)



# all:
#   coxph(formula = Surv(time, status) ~ tag_type, data = Tag_retain_Surv)
# 
# n= 95, number of events= 20 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)
# tag_typeV7 -0.4526    0.6360   0.4565 -0.991    0.322
# tag_typeV8      NA        NA   0.0000     NA       NA
# 
# exp(coef) exp(-coef) lower .95 upper .95
# tag_typeV7     0.636      1.572    0.2599     1.556
# tag_typeV8        NA         NA        NA        NA
# 
# Concordance= 0.554  (se = 0.057 )
# Likelihood ratio test= 1  on 1 df,   p=0.3
# Wald test            = 0.98  on 1 df,   p=0.3
# Score (logrank) test = 1  on 1 df,   p=0.3



###Mean retention time

#all expelled fish

Tagging_data.2 %>% 
  filter(tag_reten_days!= 92) %>% 
  summarise(n=n(),
            mean = mean(tag_reten_days ),
            sd = sd(tag_reten_days ),
            se=sd/sqrt(n))
# # A tibble: 1 x 4
# n  mean    sd    se
# <int> <dbl> <dbl> <dbl>
#   1    20  33.3  4.70  1.05


#By group
Tagging_data.2 %>% 
  group_by(tag_type) %>% 
  filter(tag_reten_days!= 92) %>% 
  summarise(n=n(),
            mean = mean(tag_reten_days ),
            sd = sd(tag_reten_days ),
            se=sd/sqrt(n))

# # A tibble: 2 x 5
# tag_type     n  mean    sd    se
# <fct>    <int> <dbl> <dbl> <dbl>
#   1 V7           8  31.9  4.79  1.69
# 2 V8          12  34.2  4.59  1.33
# 
# 



####Tag Retention survival plot####



Fit_reten <-survfit(Surv(time, status)~tag_type,
                    data =Tag_retain_Surv) 

B<-ggsurvplot(Fit_reten,legend.labs = c('V7', 'V8'),
              conf.int = F, data = Tag_retain_Surv, risk.table = T,
              risk.table.height = .3, xlab = "Time in days", ylim = c(.5, 1),
              palette = c ( 'grey', 'grey32' ), 
              break.x.by = 25,
              ggtheme = theme_classic2(base_family = "Times", base_size=14),
              linetype = "strata")
B$plot <-B$plot +scale_linetype_manual(values = c ("solid","dashed"))
x <- which(B$table$scales$find("x"))
B$table$scales$scales[[x]] <- scale_x_continuous(breaks = seq(0, 100, 25))
print(B)



#####Percent retention in salmon####

#Just want to see how many rejected overall 

summary(Tagging_data$tag_type)

# SHAM   V7   V8 
# 60   59   61 
#total = 180

Tagging_data_all<-Tagging_data %>% filter(tag_type != 'SHAM') %>%
  print()

Tagging_data_all$tag_type<-fct_drop(Tagging_data_all$tag_type, only = "SHAM")

Tagging_data_all$tag_reten_days

#120 x 11  

#so getting percents of those that retained across the whole exp,
#regardless of when dropped

table(Tagging_data_all$rejection_yanaw)

# No Yes 
# 52  63 

#5 fish were NA in the study that were the tank mortalities 
#not counting that towards the study.

#overall, 54.78% had rejected their tags over the 92 days period. 



table(Tagging_data_all$rejection_progress)

#In progress          No         Yes 
#         23          52          40 

#Of the rejected fish (N = 63), 23 were in 
#an intermediate form of rejecting the tag
#This is 36.5% of fish were 'in progress' 


Tagging_data_all %>% 
  mutate_if(is.factor, as.character) %>% 
  group_by(tag_type) %>% 
  summarise(yes = sum(grepl("Yes", rejection_yanaw)),
            no = sum(grepl("No", rejection_yanaw))) %>% print()

# # A tibble: 2 x 3
#tag_type   yes    no
#<chr>    <int> <int>
# 1 V7      29    27
# 2 V8      34    25
# > 


#by tag type, v7 51.79% rejected tag, V8 57.63 % rejected tag


#Breaking it down by the three groups instead here

Tagging_data_all %>% 
  mutate_if(is.factor, as.character) %>% # your example data was sotred as factor
  group_by(tag_type) %>% 
  summarise(yes = sum(grepl("Yes", rejection_progress)),
            no = sum(grepl("No", rejection_progress)),
            In_prog = sum(grepl("In progress", rejection_progress))) %>% print()

# # A tibble: 2 x 4
# tag_type   yes    no In_prog
# <chr>    <int> <int>   <int>
#   1 V7     16    27      13
# 2 V8       24    25      10



#Of the 'yes' fish from reject_yawnaw, for V7, 55.17% were def 'yes, 44.83 were in progress
#for V8 fish, 70.59% were def 'yes', 29.41% were 'in prog' 


###chi square for salmon tag rejection #####


Chi_prog<-table(Tagging_data_all$tag_type, Tagging_data_all$rejection_progress) %>% print()

test_prog <- chisq.test(Chi_prog) %>% print()

# Pearson's Chi-squared test
# 
# data:  Chi_prog
# X-squared = 1.9913, df = 2, p-value = 0.3695


####plots for rejection chi square####

#removing out the Na's from the fish that died prior to quantifying 

Yanaw <-Tagging_data_all %>% 
  count(tag_type, rejection_yanaw) %>%
  group_by(tag_type) %>%drop_na() %>% 
  mutate(n = n/sum(n) * 100) %>% print()
Yanaw$rejection_yanaw<-fct_relevel(Yanaw$rejection_yanaw, "Yes", "No") %>% print()


ggplot(data = Yanaw) + aes(tag_type, n, fill = rejection_yanaw, label = paste0(round(n, 2), "%")) + 
  geom_col(colour="black") +
  geom_text(position=position_stack(0.5))+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0, 105, 10))+
  scale_fill_manual(values = c('No' = "white", 'Yes' = "grey45"))+
  theme(axis.text.x = element_text(size = 14, family= 'Times', face = 'bold'),
        axis.text.y = element_text(size = 12, family= 'Times'))+
  theme(axis.title = element_text(size = 17, family= 'Times'))+
  labs(x ="Tag Type", y = "Tags Rejected (%)")+
  theme(legend.position="none")+
  annotate("text", x = 1, y = 90, label = "Yes", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 1, y = 40, label = "No", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 2, y = 90, label = "Yes", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 2, y = 40, label = "No", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 2, y = 40, label = "No", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = .48, y = 110, label = "", fontface = 'bold', family ='Times', size =8)+
  annotate("text", x = .48, y = 105, label = "A", fontface = 'bold', family ='Times', size =8)+
  theme(text=element_text(family = 'Times'))

ggsave("Overall_tags_discharged_Nov 23 2021.tiff")



Rej_prog <-Tagging_data_all %>% 
  count(tag_type, rejection_progress) %>%
  group_by(tag_type) %>%drop_na() %>% 
  mutate(n = n/sum(n) * 100) %>% print()

Rej_prog$rejection_progress<-fct_relevel(Rej_prog$rejection_progress, "Yes", "In progress", "No") %>% print()
Rej_prog

ggplot(data = Rej_prog) + aes(tag_type, n, fill = rejection_progress, label = paste0(round(n, 2), "%")) + 
  geom_col(colour="black") +
  geom_text(position=position_stack(0.5))+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(text=element_text(family = 'Times'))+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0, 105, 10))+
  scale_fill_manual(values = c('No' = "white", 'In progress' = "grey", "Yes" = 'grey45'))+
  theme(axis.text.x = element_text(size = 14, family= 'Times', face = 'bold'),
        axis.text.y = element_text(size = 12, family= 'Times'))+
  theme(axis.title = element_text(size = 17, family= 'Times'))+
  labs(x ="Tag Type", y = "Tags Rejected (%)")+
  theme(legend.position="none")+
  annotate("text", x = 1, y = 95, label = "Yes", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 2, y = 95, label = "Yes", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 1, y = 65.5, label = "In progress", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 2, y = 56.5, label = "In progress", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 1, y = 30, label = "No", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = 2, y = 30, label = "No", fontface = 'italic', family ='Times', size =5)+
  annotate("text", x = .48, y = 110, label = "", fontface = 'bold', family ='Times', size =8)


ggsave("Tag Progression_discharged_Nov 23 2021.tiff")



####Response to tag#####


#Just want to get some counts/representation of what the fate internally of the tag was 

State_bytag <-Tagging_data_all %>% 
  count(tag_type, tag_response.1) %>%
  group_by(tag_type) %>%drop_na() %>% 
  filter(tag_response.1  != 'mortality') %>% 
  mutate(Numberz = n) %>% select(-n) %>% print()

view(State_bytag)

tag_type <- c("V7", 'V8') %>% as.factor()
tag_response.1<-c("adhesion around the tag", "embedded in pyloric caeca")
Numberz <-c("0","0") %>% as.numeric()



placeholders_response <-cbind.data.frame(tag_type,tag_response.1, Numberz) %>% print()

State_bytag.1 <-bind_rows(placeholders_response, State_bytag) %>% print() 
State_bytag.2<-State_bytag.1 %>% rename('Transmitter type' = "tag_type") %>% print()


State_bytag.1$tag_response.1 <- factor(State_bytag.2$tag_response.1, levels= 
                                         c('expelled','encapsulated in mesentery','encapsulated in mesentery; enveloped in viscera',
                                           'partially encapsulated','enveloped in viscera','embedded in pyloric caeca',
                                           "adhesion around the tag",'no encapsulation'))

unique(State_bytag.1$tag_response.1)



setwd('')


ggplot(State_bytag.1,aes(x=tag_response.1,y=Numberz,fill=tag_type))+
  geom_bar(stat="identity",position="dodge")+ coord_flip()+
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = c('V7' = "grey", 'V8' = "grey32"))+
  theme(axis.text.x = element_text(size = 14, family= 'Times', face = 'bold'),
        axis.text.y = element_text(size = 14, family= 'Times'))+
  theme(axis.title = element_text(size = 17, family= 'Times'))+
  labs(x ="", y = "Number of Cases")+
  expand_limits(x=c(0,8), y=c(0, 42))+
  scale_x_discrete(labels = c('Expelled','Encapasulated, mesentary','Encapsulated, mesentary; Enveloped, viscera',
                              'Partially encapsulated','Enveloped, viscera','Embedded, pyloric caeca',
                              "Adhesion around tag",  'No encapsulation'))+
  theme(legend.title=element_blank())


ggsave("Tag_response_april 12 2023.tiff", width = 12, height = 6)




#####Salmon Meristics####



#mean body mass + SE
describe(Tag_data_all$initial.weight..g.)

# vars   n   mean    sd median trimmed   mad    min    max  range skew kurtosis   se
# X1    1 173 208.67 39.61 204.96  208.78 38.64 105.14 303.91 198.78 0.03    -0.46 3.01
# 



describe(Tag_data_all$initial.length..mm.)

# vars   n   mean    sd median trimmed   mad min max range  skew kurtosis   se
# X1    1 179 269.65 18.42    270  270.59 17.79 216 308    92 -0.42    -0.02 1.38

view(Tag_data_all)



Tag_mort %>% group_by(Tag_type) %>%
  dplyr::summarize (Mean = mean( initial.weight..g., na.rm=TRUE) )

options(pillar.sigfig=6)


# # A tibble: 3 x 2
# Tag_type    Mean
# <fct>      <dbl>
#   1 SHAM     214.868
#   2 V7       212.668
#   3 V8       198.731




#####Part 2: Meta analysis of the literature######


#This set of analyses uses the metafor package 


####Load in and clean up data#####


DF_A<-read_excel("Meta_data_filtering_Jan 20 2023_currated.xlsx") %>% print()
#318 x46


#Removes out unsuitable papers
DF_A.1 <- DF_A %>% filter(Reject == "N")  %>% print()
#162 x 46


DF_A.2<-DF_A.1 %>% clean_names() %>% print()


#take out what I need and rename so stuff

DF_A.3<-DF_A.2 %>% 
  select( "paper_id_19","year_2", "authors_21",'taxa_18',"species","common_name_2",
          "water_type","treatment", "tag_type_2", 'experiment_num', "sample_size_n",           
          "num_tag_lost", "tag_reten_perc","total_mortalities","mort_perc",  
          "exp_duration","tag_length","tag_diam","other_tags_2",
          "incision_size_max_mm",
          "tag_mass_air","fish_mass_ave","fish_mass_ave_actual", 'notes') %>%                  
  rename( paper_id = 'paper_id_19', year = year_2, author = 'authors_21', 
          tag_type = 'tag_type_2',tag_loss_N = 'num_tag_lost',
          mort_N = 'total_mortalities', incision_size = 'incision_size_max_mm',
          taxa = 'taxa_18') %>% 
  print()       

#162 x 24


DF_A.3 %>% summarise(Unique_Elements = n_distinct(paper_id)) 

#40 total papers that are included in the analysis
#writing out the data file with the filtered results

write_csv(DF_A.3, 'Tag_meta filtered_Jan 25 2023.csv')


Meta_data <- read_csv('Tag_meta filtered_Jan 25 2023.csv') %>% print()
#162x24

####Meta: Tag retention#######

####Meta: tag retention data prep#####

#dropping any controls or sham fish from the df
Tagged_fish_DF <- Meta_data %>% filter(!grepl('control|sham', tag_type)) %>% print()
#110x24


sort(Tagged_fish_DF$tag_type) %>% print() #visual check shows no control fish
str(Tagged_fish_DF)

length(Tagged_fish_DF$tag_loss_N) #110

Tagged_fish_DF %>% drop_na(tag_loss_N) %>% print()
#80


#simplifying and reorganizing DF

Tagged_fish_DF.1 <-Tagged_fish_DF %>% select(-notes, -other_tags_2) %>% 
  separate(species, c("genus", "species"), " ") %>%
  mutate(tags_retained = sample_size_n-tag_loss_N) %>% mutate(tag_BW_ratio = tag_mass_air/fish_mass_ave*100 ) %>% print()       
#110x25   



#####Meta: Tag retention estimates######


mean(Tagged_fish_DF.1$tag_reten_perc, na.rm = T)
#[1] 86.69237

range(Tagged_fish_DF.1$tag_reten_perc, na.rm = T)
#[1]  20 100

#Transforms the proportion (yi) and generates the estimate's sampling variance (vi)
Tag_retention_DF<-  escalc(measure = "PFT", ni=sample_size_n, 
                           xi=tags_retained, 
                           data =Tagged_fish_DF.1, add = 0) 

dim(Tag_retention_DF) #110x17


###Meta: Tag retention - full model#####


Reten_model_full <- rma.mv(yi, vi,mods = ~tag_length+tag_diam+incision_size+exp_duration+tag_BW_ratio, 
                           random = ~1|paper_id/experiment_num, method = "REML", data = Tag_retention_DF)
Reten_model_full

# Multivariate Meta-Analysis Model (k = 42; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0571  0.2389     22     no                 paper_id 
# sigma^2.2  0.0094  0.0970     42     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 36) = 304.5762, p-val < .0001
# 
# Test of Moderators (coefficients 2:6):
#   QM(df = 5) = 15.6590, p-val = 0.0079
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub     ??? 
# intrcpt          1.3851  0.1472   9.4127  <.0001   1.0967   1.6736  *** 
#   tag_length      -0.0044  0.0053  -0.8414  0.4001  -0.0147   0.0059      
# tag_diam         0.0079  0.0189   0.4169  0.6768  -0.0291   0.0449      
# incision_size    0.0068  0.0097   0.7006  0.4836  -0.0122   0.0259      
# exp_duration    -0.0005  0.0006  -0.9124  0.3615  -0.0016   0.0006      
# tag_BW_ratio    -0.0610  0.0203  -3.0002  0.0027  -0.1009  -0.0212   ** 
#   
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1




####Meta: Individual tag retention models######


###Tag length

Reten_model_Length <- rma.mv(yi, vi,mods = ~tag_length, 
                             random = ~1|paper_id/experiment_num, method = "REML", data = Tag_retention_DF)
Reten_model_Length

# Multivariate Meta-Analysis Model (k = 71; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0263  0.1622     32     no                 paper_id 
# sigma^2.2  0.0245  0.1564     71     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 69) = 876.3955, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 1.6438, p-val = 0.1998
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub     ??? 
# intrcpt       1.3479  0.0723  18.6467  <.0001   1.2062  1.4896  *** 
#   tag_length   -0.0036  0.0028  -1.2821  0.1998  -0.0090  0.0019      
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 



###Tag diameter

Reten_model_Diam <- rma.mv(yi, vi,mods = ~tag_diam, 
                           random = ~1|paper_id/experiment_num, method = "REML", data = Tag_retention_DF)
Reten_model_Diam

# Multivariate Meta-Analysis Model (k = 68; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0258  0.1607     32     no                 paper_id 
# sigma^2.2  0.0274  0.1656     68     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 66) = 961.0724, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 0.0059, p-val = 0.9386
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub     ??? 
# intrcpt     1.2715  0.0828  15.3506  <.0001   1.1091  1.4338  *** 
#   tag_diam   -0.0008  0.0099  -0.0771  0.9386  -0.0202  0.0187      
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# > 



####Experiment duration

Reten_model_Dur <- rma.mv(yi, vi,mods = ~exp_duration, 
                          random = ~1|paper_id/experiment_num, method = "REML", data = Tag_retention_DF)
Reten_model_Dur

# 
# Multivariate Meta-Analysis Model (k = 80; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0247  0.1572     35     no                 paper_id 
# sigma^2.2  0.0245  0.1566     80     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 78) = 861.1038, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 1.0380, p-val = 0.3083
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub     ??? 
# intrcpt         1.3008  0.0468  27.7973  <.0001   1.2091  1.3925  *** 
#   exp_duration   -0.0003  0.0003  -1.0188  0.3083  -0.0010  0.0003      
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1





####BW ratio

Reten_model_BWratio <- rma.mv(yi, vi,mods = ~tag_BW_ratio,
                              random = ~1|paper_id/experiment_num, method = "REML", data = Tag_retention_DF)
Reten_model_BWratio
# 
# Multivariate Meta-Analysis Model (k = 62; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0477  0.2185     27     no                 paper_id 
# sigma^2.2  0.0139  0.1177     62     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 60) = 584.5636, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 16.6106, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub     ??? 
# intrcpt         1.3867  0.0586  23.6678  <.0001   1.2719   1.5015  *** 
#   tag_BW_ratio   -0.0694  0.0170  -4.0756  <.0001  -0.1028  -0.0360  *** 
#   
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

####Incision size

Reten_model_incision <- rma.mv(yi, vi,mods = ~incision_size, 
                               random = ~1|paper_id/experiment_num, method = "REML", data = Tag_retention_DF)
Reten_model_incision

# 
# Multivariate Meta-Analysis Model (k = 57; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0341  0.1847     29     no                 paper_id 
# sigma^2.2  0.0188  0.1373     57     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 55) = 692.8663, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 0.3998, p-val = 0.5272
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub     ??? 
# intrcpt          1.2332  0.0840  14.6795  <.0001   1.0685  1.3978  *** 
#   incision_size    0.0040  0.0064   0.6323  0.5272  -0.0084  0.0165      
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



####Meta: Retention model - FDR ####
#good resource https://stackoverflow.com/questions/61456841/apply-fdr-correction-on-a-large-number-of-outcome-variables

#First clean up and combine model info

R_a<-tidy(Reten_model_full, conf.int = T) %>% add_column(Model = 'Full') %>% print()
R_b<-tidy(Reten_model_Length, conf.int = T) %>% add_column(Model = 'Len')%>% print()
R_c<-tidy(Reten_model_Diam, conf.int = T) %>% add_column(Model = 'Diam')%>% print()
R_d<-tidy(Reten_model_incision, conf.int = T) %>% add_column(Model = 'Inc')%>% print()
R_e<-tidy(Reten_model_Dur, conf.int = T) %>% add_column(Model = 'Dur')%>% print()
R_f<-tidy(Reten_model_BWratio, conf.int = T) %>% add_column(Model = 'BWR')%>% print()

#Combine and remove intercepts 
Reten_model_all <-bind_rows(R_a, R_b, R_c, R_d, R_e, R_f) %>% 
  filter(term != 'intercept') %>% 
  group_by(term) %>% 
  mutate(padj = p.adjust(p.value, 'BH')) %>% 
  ungroup() %>% 
  print()

write_csv(Reten_model_all, '')



####Meta: Tag retention figures######

windowsFonts(Times=windowsFont("TT Times New Roman")) 
setwd('')




#first, back transform the data for individual porportions and their CI's

dat.back <- summary(Tag_retention_DF, transf=transf.ipft, ni=Tag_retention_DF$sample_size)
dat.back


####Retention vs BW ratio scatter




ggplot()+
  geom_point(data = dat.back, aes(x=tag_BW_ratio, y=tag_reten_perc, colour = tag_type, shape = tag_type), size=3)+
  geom_smooth(data = dat.back, aes(x=tag_BW_ratio, y=tag_reten_perc),method='lm', formula= y~x, color="#00798c", se=F)+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5),
        axis.text.x = element_text(size = 16, color="black"), 
        axis.text.y = element_text(size = 16,color="black"), 
        axis.title = element_text(size = 18, color="black"),
        plot.margin = unit(c(.25,.35,.25,.25), "in"),
        text=element_text(family = "Times", colour = "black"))+
  xlab("Tag mass:fish body mass (%)") + 
  ylab("Tag retention (%)")+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,110))+
  scale_colour_manual(values = c('black', 'red', '#bc5090', '#ffa600', '#0FB859'),
                      labels = c("Acoustic", "Archival", 'Acoustic, dummy', 'Radio, dummy', 'PIT'),
                      name = 'Tag type')+
  scale_shape_manual(values = c(16, 15, 16,18, 17 ),
                     labels = c("Acoustic", "Archival", 'Acoustic, dummy', 'Radio, dummy', 'PIT'),
                     name = 'Tag type')


ggsave( 'plotsReten_BW_scatter_April 25 2023.tiff',
        width = 10, height = 6, dpi = 300, units = "in") 




###Forest plot


dat.back.1 <- as.tibble(dat.back)


dat.back.1$paper_id<-as.factor(dat.back.1$paper_id)
dat.back.1$experiment_num <-as.factor(dat.back.1$experiment_num)
str(dat.back.1)


dat.back.1$ID <- seq.int(nrow(dat.back.1))


dat.back.1<-dat.back.1 %>% unite("Authoryear",author, year, sep = ' ', remove = F,  ) %>% print()
#110 x 30

dat.back.2<-dat.back.1 %>% drop_na(tag_loss_N) %>% print()
#80 x 30

dat.back.3 <-dat.back.2 %>% select(paper_id, Authoryear, yi,ci.lb,ci.ub) %>% 
  arrange(paper_id) %>% 
  drop_na() %>% print()
#80 x 4

dat.back.3$ID <- seq.int(nrow(dat.back.3))

dat.back.3$ID <-as.factor(dat.back.3$ID)
dat.back.3$Authoryear <- as.factor(dat.back.3$Authoryear)

write.csv(dat.back.3, 'Base_reten_labels.csv')


labels <-c(" ",
           '507',
           " ",
           " ",
           " ",
           '564',
           " ",
           " ",
           " ",
           '565',
           '579',
           '604',
           " ",
           '611',
           '640',
           " ",
           " ",
           '652',
           '660',
           " ",
           " ",
           " ",
           '666',
           " ",
           " ",
           '741',
           " ",
           " ",
           " ",
           '742',
           " ",
           '763',
           " ",
           '766',
           " ",
           '782',
           " ",
           '801',
           " ",
           '810',
           '820',
           " ",
           '831',
           " ",
           " ",
           '841',
           '849',
           '864',
           " ",
           " ",
           '871',
           " ",
           " ",
           " ",
           " ",
           '881',
           " ",
           '1081',
           " ",
           '1102',
           '1200',
           " ",
           '1203',
           " ",
           " ",
           " ",
           '1273',
           '1419',
           " ",
           " ",
           '1520',
           " ",
           '1531',
           " ",
           " ",
           '1543',
           " ",
           '1697',
           " ",
           '1867') %>% print()


dat.back.3$boxes <-c("white",
                     "white",
                     'grey34',
                     'grey34',
                     'grey34',
                     'grey34',
                     "white",
                     "white",
                     "white",
                     "white",
                     'grey34',
                     "white",
                     'grey34',
                     'grey34',
                     "white",
                     'grey34',
                     'grey34',
                     'grey34',
                     "white",
                     'grey34',
                     'grey34',
                     'grey34',
                     'grey34',
                     "white",
                     "white",
                     "white",
                     'grey34',
                     'grey34',
                     'grey34',
                     'grey34',
                     "white",
                     "white",
                     'grey34',
                     'grey34',
                     "white",
                     "white",
                     'grey34',
                     'grey34',
                     "white",
                     "white",
                     'grey34',
                     "white",
                     "white",
                     'grey34',
                     'grey34',
                     'grey34',
                     "white",
                     'grey34',
                     "white",
                     "white",
                     "white",
                     'grey34',
                     'grey34',
                     'grey34',
                     'grey34',
                     'grey34',
                     "white",
                     "white",
                     'grey34',
                     'grey34',
                     "white",
                     'grey34',
                     'grey34',
                     "white",
                     "white",
                     "white",
                     "white",
                     'grey34',
                     "white",
                     "white",
                     "white",
                     'grey34',
                     'grey34',
                     "white",
                     "white",
                     "white",
                     'grey34',
                     'grey34',
                     "white",
                     "white")
str(dat.back.3)
length(dat.back.3$boxes)



ggplot(data=dat.back.3,
       aes(x = ID,y = yi, ymin = ci.lb, ymax = ci.ub ))+
  geom_vline(aes(xintercept =1:80 , colour = boxes), size = 6)+
  scale_color_manual(values =c("white", "grey83")) +
  geom_pointrange()+coord_flip()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(text=element_text(family = "Times", colour = "black"))+
  xlab("") + ylab("Proportion retained (95% confidence interval)")+
  scale_x_discrete(breaks=c(dat.back.3$ID,1:76), labels=c(labels,1:76), position = "top", expand=c(.027, 0))+ #postion top is how the x axis is now on the right
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 12))+  
  theme(axis.title = element_text(size = 18), axis.title.x = element_text(vjust=-2))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(1,.5,.5,.5), "in"))+
  labs(tag = "Study") +
  theme(plot.tag.position = c(.97,1.05 ), plot.tag = element_text(size = 18))



ggsave("Tag retention point estimates_Jan 31 2023.tiff", 
       width = 8, height = 16.8, dpi = 300, units = "in", device='tiff') 







#####Meta: Mortality estimates########


#overall motility from all of the tagged fish
mean(Tagged_fish_DF.1$mort_perc, na.rm =T)
#[1] 12.39973

range(Tagged_fish_DF.1$mort_perc, na.rm =T)
#[1]  0.0 90


#Transform data and calculate variances


Tag_mort_DF<-  escalc(measure = "PFT", ni=sample_size_n, 
                      xi=mort_N, 
                      data =Tagged_fish_DF.1, add = 0)

Tag_mort_DF



###Meta: Mortality - Full model



Mort_model_full <- rma.mv(yi, vi,mods = ~tag_length+tag_diam+incision_size+exp_duration+tag_BW_ratio, 
                          random = ~1|paper_id/experiment_num, method = "REML", data = Tag_mort_DF) 

Mort_model_full

# Multivariate Meta-Analysis Model (k = 43; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0202  0.1420     21     no                 paper_id 
# sigma^2.2  0.0252  0.1587     43     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 37) = 259.0985, p-val < .0001
# 
# Test of Moderators (coefficients 2:6):
#   QM(df = 5) = 9.5139, p-val = 0.0902
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    ??? 
# intrcpt          0.1506  0.1299   1.1595  0.2463  -0.1040  0.4053     
# tag_length      -0.0088  0.0045  -1.9430  0.0520  -0.0176  0.0001   . 
# tag_diam         0.0442  0.0161   2.7488  0.0060   0.0127  0.0756  ** 
#   incision_size   -0.0015  0.0080  -0.1909  0.8486  -0.0173  0.0142     
# exp_duration    -0.0000  0.0005  -0.0155  0.9876  -0.0010  0.0010     
# tag_BW_ratio    -0.0023  0.0218  -0.1047  0.9166  -0.0451  0.0405     
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1




#K is number of studies combined in the analysis so 42 here



####Meta: Mortality - individual models######


###Tag length

Mort_model_Length <- rma.mv(yi, vi,mods = ~tag_length, 
                            random = ~1|paper_id/experiment_num, method = "REML", data = Tag_mort_DF) 
Mort_model_Length

# Multivariate Meta-Analysis Model (k = 80; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0271  0.1646     32     no                 paper_id 
# sigma^2.2  0.0105  0.1025     80     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 78) = 749.9108, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 0.0193, p-val = 0.8896
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub     ??? 
# intrcpt       0.2585  0.0567   4.5564  <.0001   0.1473  0.3697  *** 
#   tag_length   -0.0002  0.0018  -0.1388  0.8896  -0.0037  0.0032      




###Tag diameter

Mort_model_Diam <- rma.mv(yi, vi,mods = ~tag_diam, 
                          random = ~1|paper_id/experiment_num, method = "REML", data = Tag_mort_DF)
Mort_model_Diam


# Multivariate Meta-Analysis Model (k = 77; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0294  0.1715     32     no                 paper_id 
# sigma^2.2  0.0073  0.0852     77     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 75) = 700.0846, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 0.8897, p-val = 0.3456
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub    ??? 
# intrcpt     0.2036  0.0678  3.0025  0.0027   0.0707  0.3365  ** 
#   tag_diam    0.0069  0.0073  0.9432  0.3456  -0.0074  0.0213     
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 


####Experiment duration

Mort_model_Dur <- rma.mv(yi, vi,mods = ~exp_duration, 
                         random = ~1|paper_id/experiment_num, method = "REML", data = Tag_mort_DF)
Mort_model_Dur

# Multivariate Meta-Analysis Model (k = 94; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0499  0.2233     35     no                 paper_id 
# sigma^2.2  0.0087  0.0930     94     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 92) = 4064.1908, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 0.0074, p-val = 0.9316
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub     ??? 
# intrcpt         0.2845  0.0453   6.2780  <.0001   0.1957  0.3733  *** 
#   exp_duration   -0.0000  0.0002  -0.0859  0.9316  -0.0004  0.0004      
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


####BW ratio

Mort_model_BWratio <- rma.mv(yi, vi,mods = ~tag_BW_ratio,
                             random = ~1|paper_id/experiment_num, method = "REML", data = Tag_mort_DF)
Mort_model_BWratio


# Multivariate Meta-Analysis Model (k = 67; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0365  0.1909     26     no                 paper_id 
# sigma^2.2  0.0085  0.0921     67     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 65) = 441.8434, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 0.0575, p-val = 0.8105
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub     ??? 
# intrcpt         0.2529  0.0514  4.9220  <.0001   0.1522  0.3535  *** 
#   tag_BW_ratio    0.0036  0.0149  0.2397  0.8105  -0.0256  0.0327      
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 

####Incision size

Mort_model_incision <- rma.mv(yi, vi,mods = ~incision_size, 
                              random = ~1|paper_id/experiment_num, method = "REML", data = Tag_mort_DF)
Mort_model_incision

# 
# Multivariate Meta-Analysis Model (k = 60; method: REML)
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed                   factor 
# sigma^2.1  0.0596  0.2442     28     no                 paper_id 
# sigma^2.2  0.0136  0.1166     60     no  paper_id/experiment_num 
# 
# Test for Residual Heterogeneity:
#   QE(df = 58) = 3829.7383, p-val < .0001
# 
# Test of Moderators (coefficient 2):
#   QM(df = 1) = 0.0791, p-val = 0.7785
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb   ci.ub    ??? 
# intrcpt          0.3246  0.1009   3.2169  0.0013   0.1268  0.5223  ** 
#   incision_size   -0.0021  0.0074  -0.2813  0.7785  -0.0165  0.0124     
# 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 




####Mort model - FDR ####

#First clean up and combine model info


M_a<-tidy(Mort_model_full, conf.int = T) %>% add_column(Model = 'Full') %>% print()
M_b<-tidy(Mort_model_Length, conf.int = T) %>% add_column(Model = 'Len')%>% print()
M_c<-tidy(Mort_model_Diam, conf.int = T) %>% add_column(Model = 'Diam')%>% print()
M_d<-tidy(Mort_model_incision, conf.int = T) %>% add_column(Model = 'Inc')%>% print()
M_e<-tidy(Mort_model_Dur, conf.int = T) %>% add_column(Model = 'Dur')%>% print()
M_f<-tidy(Mort_model_BWratio, conf.int = T) %>% add_column(Model = 'BWR')%>% print()

#Combine and remove intercepts 
Mort_model_all <-bind_rows(M_a, M_b, M_c, M_d, M_e, M_f) %>% 
  filter(term != 'intercept') %>% 
  group_by(term) %>% 
  mutate(padj = p.adjust(p.value, 'BH')) %>% 
  ungroup() %>% 
  print()

write_csv(Mort_model_all, '')



####Meta: Mortality figures######


windowsFonts(Times=windowsFont("TT Times New Roman")) 
setwd('')


#first, back transform the data for individual proportions and their CI's

dat.back_mort <- summary(Tag_mort_DF, transf=transf.ipft, ni=Tag_mort_DF$sample_size)
dat.back_mort



###forest plot


dat.back_mort.1 <- as.tibble(dat.back_mort)


dat.back_mort.1$Paper_ID<-as.factor(dat.back_mort.1$paper_id)
str(dat.back_mort.1)
dat.back_mort.1$ID <- seq.int(nrow(dat.back_mort.1))

dat.back_mort.1<-dat.back_mort.1 %>% unite("authoryear",author, year, sep = ' ', remove = F,  ) %>% print()
#110x31

dat.back_mort.2<-dat.back_mort.1 %>% drop_na(mort_N) %>% print()
#97x31

dat.back_mort.3 <-dat.back_mort.2 %>%
  select(paper_id, authoryear, yi,ci.lb,ci.ub) %>% 
  arrange(paper_id) %>% 
  drop_na() %>% print()
#97x5

dat.back_mort.3$ID <- seq.int(nrow(dat.back_mort.3)) %>% print()




dat.back_mort.3$ID <-as.factor(dat.back_mort.3$ID)
dat.back_mort.3$authoryear <- as.factor(dat.back_mort.3$authoryear)

write.csv(dat.back_mort.3, "mort_lables_feb 1 2023.csv")



labels_mort <-c('505',
                " ",
                '507',
                " ",
                " ",
                " ",
                '564',
                " ",
                " ",
                " ",
                '565',
                '579',
                " ",
                " ",
                '591',
                '604',
                " ",
                '611',
                '640',
                " ",
                " ",
                " ",
                '666',
                " ",
                " ",
                '699',
                " ",
                " ",
                '741',
                " ",
                " ",
                " ",
                '742',
                " ",
                '763',
                " ",
                '766',
                " ",
                '782',
                " ",
                '801',
                " ",
                '810',
                '820',
                " ",
                " ",
                '831',
                " ",
                " ",
                " ",
                " ",
                " ",
                '841',
                '849',
                '864',
                " ",
                " ",
                '871',
                " ",
                " ",
                " ",
                " ",
                '881',
                " ",
                " ",
                " ",
                " ",
                '893',
                " ",
                '1102',
                '1200',
                " ",
                " ",
                " ",
                " ",
                " ",
                '1203',
                " ",
                " ",
                " ",
                " ",
                " ",
                " ",
                '1273',
                '1419',
                " ",
                " ",
                '1520',
                " ",
                '1531',
                " ",
                " ",
                '1543',
                " ",
                '1697',
                " ",
                '1867') %>% print()




dat.back_mort.3$boxes <-c("white",
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          "white",
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          'grey34',
                          'grey34',
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          "white",
                          "white",
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          'grey34',
                          'grey34',
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          "white",
                          "white",
                          "white",
                          "white",
                          'grey34',
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          "white",
                          "white",
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          'grey34',
                          "white",
                          "white",
                          "white",
                          "white",
                          "white",
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          'grey34',
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          'grey34',
                          'grey34',
                          'grey34',
                          "white",
                          "white",
                          'grey34',
                          'grey34')



ggplot(data=dat.back_mort.3,
       aes(x = ID,y = yi, ymin = ci.lb, ymax = ci.ub ))+
  geom_vline(aes(xintercept =1:97 , colour = boxes), size = 5.5)+
  scale_color_manual(values =c("white", "grey83"))+
  geom_pointrange()+coord_flip()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(text=element_text(family = "Times", colour = "black"))+
  xlab("") + ylab("Proportion mortality (95% confidence interval)")+
  scale_x_discrete(breaks=c(dat.back_mort.3$ID,1:91), labels=c(labels_mort,1:91), position = "top", expand=c(.03, 0))+ #postion top is how the x axis is now on the right
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 12))+  
  theme(axis.title = element_text(size = 18), axis.title.x = element_text(vjust=-2))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(1,.5,.5,.5), "in"))+
  labs(tag = "Study") +
  theme(plot.tag.position = c(.97,1.05 ), plot.tag = element_text(size = 18))



ggsave("Mortality point estimates_Feb 1 2023.tiff", width = 8, height = 16.8, dpi = 300, units = "in", device='tiff') 




































