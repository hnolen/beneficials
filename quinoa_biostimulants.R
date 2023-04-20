#####QDM017 analyses

##working with the *_nona.csv files in the OneDrive QDM-017 folder

##Exp design - 7 x 2 factorial x 6 replicates CRD - 
#analyze as RCBD using exp_rep as blocking factor

library(car)
library(agricolae)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(DescTools)
library(data.table)
library(stats)
library(ggsignif)

#make sure everything is numeric
spad_dat$spad<-as.numeric(spad_dat$spad)


tb_dat$tot_biomass<-as.numeric(tb_dat$tot_biomass)

ag_dat$ag_biomass<-as.numeric(ag_dat$ag_biomass)


bg_dat$bg_biomass<-as.numeric(bg_dat$bg_biomass)











###SPAD MODEL - high effect of block

spad_mod<-lm(spad ~ stress*inoculant*exp_rep, spad_dat)

Anova(spad_mod) ##stress p < 0.001

##SPAD - assumption checking
#need to check if these models fit ANOVA assumptions: shapiro wilks, levene's test (this is assessing experiments separately so there is no blocking factor)
#resids x preds plot, generate resids for SW test:
spad$resids<-residuals(spad_mod)
spad$preds<-predict(spad_mod)
spad$sq_preds<-spad$preds^2
plot(resids ~ preds, data = spad) # 

#test for normality
shapiro.test(spad$resids) # p-value = 0.2

#test for homogeneity of variance
leveneTest(spad ~ treatment, data = spad_dat) # p-value = 0.07







##Comparing microbes vs water in nutrient, drought, and no stress
mw_spad<-subset(spad_dat, inoculant != "Lalrise")

mw_spad_mod<-lm(spad ~ inoculant*stress + exp_rep, mw_spad)
Anova(mw_spad_mod) # inoculant p = 0.8, stress p < 0.001


##Comparing Lalrise vs water in nutrient, drought, and no stress
lal_spad<-subset(spad_dat, inoculant == "Lalrise")
w_spad<-subset(spad_dat, inoculant == "a1water")
lal_w_spad<-rbind(lal_spad, w_spad)

lalw_spad_mod<-lm(spad ~ inoculant*stress +exp_rep, lal_w_spad)
Anova(lalw_spad_mod) # inoculant:stress p = 0.03

lalw_spad_mod<-lm(spad ~ treatment +exp_rep, lal_w_spad)
tuk<-HSD.test(lalw_spad_mod, "treatment")

p_lalvwater<-ggplot(lal_w_spad, aes(x=stress, y = spad, fill = inoculant)) +
  geom_boxplot() +
  ylim(0,55) +
  theme_classic() +
  scale_x_discrete(labels = c("Drought stress", "No stress", "Nutrient stress")) +
  xlab("") +
  ylab("SPAD") +
  theme(axis.title.y = element_text(size = 24, color = "black")) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  stat_signif(annotation = "*", textsize = 10, xmin = 2.8, xmax = 3.19, y_position = 40) +
  theme(legend.text = element_text(size =18)) +
  theme(legend.title = element_text(size =18, face = "bold")) +
  scale_fill_discrete(labels = c("Sterile water", "Lalrise Start SC"))







#########======================================================##########

##TOTAL BIOMASS MODEL - no effect of block

#########======================================================##########

tot_mod<-lm(tot_biomass ~ stress*inoculant*exp_rep, tb_dat)
 
Anova(tot_mod) #stress p < 0.001



##TOTAL BIOMASS - assumption checking
tb_dat$resids<-residuals(tot_mod)
tb_dat$preds<-predict(tot_mod)
tb_dat$sq_preds<-tb_dat$preds^2
plot(resids ~ preds, data = tb_dat) # looks like a funnel

#shapiro wilk - Tests for normality of residuals
shapiro.test(tb_dat$resids) # p-value = 1.925e-10


#levene's test - Tests for homogeneity of variance
leveneTest(tot_biomass ~ treatment, data = tb_dat) # p-value = 0.07


##TOT BIOMASS - need to transform - square root transformation
tb_dat$trans_tot_biomass<-sqrt(0.5+tb_dat$tot_biomass)

#TOTAL BIOMASS - TRANSFORMED MODEL

trans_tbmod<-lm(trans_tot_biomass ~ stress*inoculant + exp_rep, tb_dat)
Anova(trans_tbmod) #stress p < 0.001 


##just looking at stress
tb_trans_stress_mod<-lm(tot_biomass ~ stress + exp_rep, tb_dat)
Anova(tb_trans_stress_mod)
tuk<-HSD.test(tb_trans_stress_mod, "stress")


#transformed shapiro wilk
tb_dat$trans_resids<-residuals(trans_tbmod)
tb_dat$trans_preds<-predict(trans_tbmod)
tb_dat$trans_sq_preds<-tb_dat$trans_preds^2
plot(trans_resids ~ trans_preds, data = tb_dat) # looks like a funnel

#shapiro wilk - Tests for normality of residuals
shapiro.test(tb_dat$trans_resids) # p-value = 1.05e-8

leveneTest(trans_tot_biomass ~ treatment, data = tb_dat) #p = 0.077



##Comparing microbes vs water in nutrient, drought, and no stress
mw_tot<-subset(tb_dat, inoculant != "Lalrise")

mw_tb_mod<-lm(trans_tot_biomass ~ inoculant*stress + exp_rep, mw_tot)
Anova(mw_tb_mod) # inoculant p = 0.5, stress p < 0.001


##Comparing Lalrise vs water in nutrient, drought, and no stress
lal_tb<-subset(tb_dat, inoculant == "Lalrise")
w_tb<-subset(tb_dat, inoculant == "a1water")
lal_w_tb<-rbind(lal_tb, w_tb)

lalw_tb_mod<-lm(trans_tot_biomass ~ inoculant*stress + exp_rep, lal_w_tb)
Anova(lalw_tb_mod) # inoculant p = 0.3, stress p = 0.1

ggplot(lal_w_tb, aes(x = stress, y = tot_biomass, fill = inoculant)) +
  geom_boxplot()

ggplot(mw_tot, aes(x = stress, y = tot_biomass, fill = inoculant)) +
  geom_boxplot()


lal_w_tb<- subset(lal_w_tb, tot_biomass != 0.5559) #removing outlier

tb_lalvwater<-ggplot(lal_w_tb, aes(x=stress, y = tot_biomass, fill = inoculant)) +
  geom_boxplot() +
  theme_classic() +
  scale_x_discrete(labels = c("Drought stress", "No stress", "Nutrient stress")) +
  xlab("") +
  ylab("Total biomass (g)") +
  theme(axis.title.y = element_text(size = 24, color = "black")) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =18)) +
  theme(legend.title = element_text(size =18, face = "bold")) +
  scale_fill_discrete(labels = c("Sterile water", "Lalrise Start SC"))



##putting toghther belowground and total biomass plots
ggarrange(tb_lalvwater, bg_lalvwater, nrow = 1, ncol =2, common.legend = TRUE, legend = "right")




#########======================================================##########

##ABOVEGROUND BIOMASS MODEL 


#########======================================================##########

ag_mod<-lm(ag_biomass ~ stress*inoculant*exp_rep, ag_dat)

Anova(ag_mod) #stress p< 0.001


##ABOVEGROUND BIOMASS - assumption checking
ag_dat$resids<-residuals(ag_mod)
ag_dat$preds<-predict(ag_mod)
ag_dat$sq_preds<-ag_dat$preds^2
plot(resids ~ preds, data = ag_dat) # looks like a funnel

#shapiro wilk - Tests for normality of residuals
shapiro.test(ag_dat$resids) # p-value = 0.001


#levene's test - Tests for homogeneity of variance
#library(car)
leveneTest(ag_biomass ~ treatment, data = ag_dat) # p-value = 0.02

##ABOVEGROUND BIOMASS - need to transform - square root transformation
ag_dat$trans_ag_biomass<-sqrt(0.5+ag_dat$ag_biomass)

#ABOVEGROUND BIOMASS - TRANSFORMED MODEL
trans_agmod<-lm(trans_ag_biomass ~ stress*inoculant + exp_rep, ag_dat)
Anova(trans_agmod) #stress p < 0.001 



#transformed shapiro wilk
ag_dat$trans_resids<-residuals(trans_agmod)
ag_dat$trans_preds<-predict(trans_agmod)
ag_dat$trans_sq_preds<-ag_dat$trans_preds^2
plot(trans_resids ~ trans_preds, data = ag_dat) # looks like a funnel

#shapiro wilk - Tests for normality of residuals
shapiro.test(ag_dat$trans_resids) # p-value = 0.003 

leveneTest(trans_ag_biomass ~ treatment, data = ag_dat) #p = 0.02




# Have to do a Kruskal-Wallis test since we have the issue of a significant levene's test
kruskal.test(ag_bl_adj ~ stress, ag_dat) #p < 0.001
pairwise.wilcox.test(ag_dat$ag_bl_adj, ag_dat$stress, p.adjust.method = "BH")

kruskal.test(ag_bl_adj ~ inoculant, ag_dat) #p = 0.67
pairwise.wilcox.test(ag_dat$ag_bl_adj, ag_dat$inoculant, p.adjust.method = "BH")


## Subset by stress
dr_ag<-subset(ag_dat, stress == "drought")
kruskal.test(ag_bl_adj ~ inoculant, dr_ag) #p = 0.33
pairwise.wilcox.test(dr_ag$ag_bl_adj, dr_ag$inoculant, p.adjust.method = "BH")

#subset into just lalrise and water
dr_ag_lal<-subset(dr_ag, inoculant == "Lalrise")
dr_ag_w<-subset(dr_ag, inoculant == "a1water")
dr_ag_lal_w<-rbind(dr_ag_lal, dr_ag_w)
kruskal.test(ag_bl_adj ~ inoculant, dr_ag_lal_w) #p = 0.03
pairwise.wilcox.test(dr_ag_lal_w$ag_bl_adj, dr_ag_lal_w$inoculant, p.adjust.method = "BH")

##plot
ggplot(dr_ag_lal_w, aes(x = inoculant, y = ag_biomass)) +
  geom_boxplot()


nu_ag<-subset(ag_dat, stress == "nutrient")
kruskal.test(ag_bl_adj ~ inoculant, nu_ag) #p = 0.9
pairwise.wilcox.test(nu_ag$ag_bl_adj, nu_ag$inoculant, p.adjust.method = "BH")

#subset into just lalrise and water
nu_ag_lal<-subset(nu_ag, inoculant == "Lalrise")
nu_ag_w<-subset(nu_ag, inoculant == "a1water")
nu_ag_lal_w<-rbind(nu_ag_lal, nu_ag_w)
kruskal.test(ag_bl_adj ~ inoculant, nu_ag_lal_w) #p = 0.69

no_ag<-subset(ag_dat, stress == "none")
kruskal.test(ag_bl_adj ~ inoculant, no_ag) #p = 0.29

#subset into just lalrise and water
no_ag_lal<-subset(no_ag, inoculant == "Lalrise")
no_ag_w<-subset(no_ag, inoculant == "a1water")
no_ag_lal_w<-rbind(no_ag_lal, no_ag_w)
kruskal.test(ag_bl_adj ~ inoculant, no_ag_lal_w) #p = 0.56



lal_ag<-subset(ag_dat, inoculant == "Lalrise")
w_ag<-subset(ag_dat, inoculant == "a1water")
lal_w_ag<-rbind(lal_ag, w_ag)

##plotting that there were diffs in ag_biomass in lalrise vs water under drought stress but not nutrient or none
ag_lalvwater<-ggplot(lal_w_ag, aes(x=stress, y = ag_biomass, fill = inoculant)) +
  geom_boxplot() +
  theme_classic() +
  scale_x_discrete(labels = c("Drought stress", "No stress", "Nutrient stress")) +
  xlab("") +
  ylab("Aboveground biomass (g)") +
  theme(axis.title.y = element_text(size = 24, color = "black")) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  stat_signif(annotation = "*", textsize = 10, xmin = 0.8, xmax = 1.19, y_position = 0.052) +
  theme(legend.text = element_text(size =18)) +
  theme(legend.title = element_text(size =18, face = "bold")) +
  scale_fill_discrete(labels = c("Sterile water", "Lalrise Start SC"))



ag_dr<-subset(ag_dat, stress == "drought")
ag_dr_lal<-subset(ag_dr, inoculant == "Lalrise")
median(ag_dr_lal$ag_biomass) # 0.0138
ag_dr_w<-subset(ag_dr, inoculant == "a1water")
median(ag_dr_w$ag_biomass) # 0.0095

#########======================================================##########

##BELOWGROUND BIOMASS MODEL 

#########======================================================##########

bg_mod<-lm(bg_biomass ~ stress*inoculant*exp_rep, bg_dat)
 
Anova(bg_mod) #stress p = 0.003, exp_rep p = 0.005



##BELOWGROUND BIOMASS - assumption checking
bg_dat$resids<-residuals(bg_mod)
bg_dat$preds<-predict(bg_mod)
bg_dat$sq_preds<-bg_dat$preds^2
plot(resids ~ preds, data = bg_dat) # looks like a funnel

#shapiro wilk - Tests for normality of residuals
shapiro.test(bg_dat$resids) # p-value = 5.845e-15


#levene's test - Tests for homogeneity of variance
#library(car)
leveneTest(bg_biomass ~ treatment, data = bg_dat) # p-value = 0.22







##plotting  lalrise vs water under stress
#taking out outlier
lal_w_bg<-subset(lal_w_bg, bg_biomass != 0.4827)


bg_lalvwater<-ggplot(lal_w_bg, aes(x=stress, y = bg_biomass, fill = inoculant)) +
  geom_boxplot() +
  theme_classic() +
  scale_x_discrete(labels = c("Drought stress", "No stress", "Nutrient stress")) +
  xlab("") +
  ylab("Belowground biomass (g)") +
  theme(axis.title.y = element_text(size = 24, color = "black")) +
  theme(axis.text.x = element_text(size =20, color = "black")) +
  theme(axis.text.y = element_text(size =20, color = "black")) +
  theme(legend.text = element_text(size =18)) +
  theme(legend.title = element_text(size =18, face = "bold")) +
  scale_fill_discrete(labels = c("Sterile water", "Lalrise Start SC"))










#########======================================================##########

##WILT DATA 

#########======================================================##########

wilt_dr_dat<- subset(wilt_dat, stress == "drought")

kruskal.test(wilt_cat ~ inoculant, wilt_dr_dat) #p = 0.9


#########======================================================##########

##Visualizing data

#########======================================================##########


p1<-ggplot(ag_dat, aes(x=inoculant, y = ag_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 


p2<-ggplot(bg_dat, aes(x=inoculant, y = bg_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

p3<-ggplot(tb_dat, aes(x=inoculant, y = tot_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

p4<-ggplot(spad_dat, aes(x=inoculant, y = spad, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

ptot<-ggarrange(p1, p2, p3, p4, ncol = 2, nrow =2, common.legend = TRUE)


##visualizing effects of stress

p5<-ggplot(ag_dat, aes(x=stress, y = ag_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 


p6<-ggplot(bg_dat, aes(x=stress, y = bg_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

p7<-ggplot(tb_dat, aes(x=stress, y = tot_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

p8<-ggplot(spad_dat, aes(x=stress, y = spad, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

ptot2<-ggarrange(p5, p6, p7, p8, ncol = 2, nrow =2, legend = "none")







##Full models all showed high effects of stress - subsetting data by stress to evalute effect of inoculum
##subsetting into two different datasets:
## Drought stress vs no stress
## Nutrient stress vs no stress

ag_dr<- subset(ag_dat, stress != "nutrient")
bg_dr<- subset(bg_dat, stress != "nutrient")
tb_dr<- subset(tb_dat, stress != "nutrient")
spad_dr<- subset(spad_dat, stress != "nutrient")


ag_n<- subset(ag_dat, stress != "drought")
bg_n<- subset(bg_dat, stress != "drought")
tb_n<- subset(tb_dat, stress != "drought")
spad_n<- subset(spad_dat, stress != "drought")

#spad models
dr_spad_mod<-lm(spad_bl_adj ~ stress*inoculant + exp_rep, spad_dr)
Anova(dr_spad_mod) #NS stress p = 0.55, inoculant p = 0.81 - drought did not impact spad values?

n_spad_mod<-lm(spad_bl_adj ~ stress*inoculant + exp_rep, spad_n)
Anova(n_spad_mod) #stress p < 0.001, inoculant p = 0.7 - nutrient stress did impact spad values


##subsetting just lalrise and water out of n_sub data
n_lal_dat<-subset(spad_n, inoculant == "Lalrise")
n_w_dat<-subset(spad_n, inoculant == "a1water")
n_l_w_dat<-rbind(n_lal_dat, n_w_dat)
n_l_w_dat<-na.omit(n_l_w_dat)
n_lw_spad_mod<-lm(spad_bl_adj ~ treatment + exp_rep, n_l_w_dat)
Anova(n_lw_spad_mod) #p < 0.001
tuk<-HSD.test(n_lw_spad_mod, "treatment")

##subsetting just nutrient stress from the n_l_w_dat
nut_l_w_dat<- subset(n_l_w_dat, stress == "nutrient")
nut_lw_spad_mod<-lm(spad_bl_adj ~ treatment + exp_rep, nut_l_w_dat)
Anova(nut_lw_spad_mod) #difference between sterile water v lalrise under nutrient stress - p = 0.01*

p_for_anissa<-ggplot(n_l_w_dat, aes(x=stress, y = spad, fill = inoculant)) +
  geom_boxplot() +
  ylim(0,55) +
  theme_classic() +
  scale_x_discrete(labels = c("No stress", "Nutrient stress")) +
  xlab("") +
  ylab("SPAD") +
  theme(axis.title.y = element_text(size = 20, color = "black")) +
  theme(axis.text.x = element_text(size =16, color = "black")) +
  theme(axis.text.y = element_text(size =16, color = "black")) +
  stat_signif(annotation = "*", textsize = 10, xmin = 1.8, xmax = 2.19, y_position = 40) +
  theme(legend.text = element_text(size =18)) +
  theme(legend.title = element_text(size =18, face = "bold")) +
  scale_fill_discrete(labels = c("Sterile water", "Lalrise Start SC"))



##visualizing spad values across inoculants from nutrient stress + none dataset
p17<-ggplot(n_sub, aes(x=inoculant, y = spad, fill = stress)) +
  geom_boxplot() +
  theme_classic() +
  scale_x_discrete(labels = c("Sterile water", "Lalrise", "NFB", "PSB", "NFB + PSB" )) +
  xlab("") +
  ylab("SPAD") +
  theme(axis.title.y = element_text(size = 16, color = "black")) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black"))




#total biomass models
dr_tb_mod<-lm(tot_biomass ~ stress*inoculant + exp_rep, dr_sub)
Anova(dr_tb_mod) # stress p < 0.001

n_tb_mod<-lm(tot_biomass ~ stress*inoculant + exp_rep, n_sub)
Anova(n_tb_mod) #stress p < 0.001




#aboveground biomass models
dr_ab_mod<-lm(ag_bl_adj ~ stress*inoculant + exp_rep, dr_sub)
Anova(dr_ab_mod) # stress p < 0.001

n_ab_mod<-lm(ag_bl_adj ~ stress*inoculant + exp_rep, n_sub)
Anova(n_ab_mod) #stress p < 0.001



##visualizing aboveground biomass values across inoculants from drought stress + none dataset
p18<-ggplot(dr_sub, aes(x=inoculant, y = ag_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() +
  scale_x_discrete(labels = c("Sterile water", "Lalrise", "NFB", "PSB", "NFB + PSB" )) +
  xlab("") +
  ylab("SPAD") +
  theme(axis.title.y = element_text(size = 16, color = "black")) +
  theme(axis.text.x = element_text(size =12, color = "black")) +
  theme(axis.text.y = element_text(size =12, color = "black"))





#belowground biomass models
dr_bg_mod<-lm(bg_bl_adj ~ stress*inoculant + exp_rep, dr_sub)
Anova(dr_bg_mod) # NS

n_bg_mod<-lm(bg_bl_adj ~ stress*inoculant + exp_rep, n_sub)
Anova(n_bg_mod) #stress p = 0.001


DunnettTest(dr_sub$bg_bl_adj, dr_sub$treatment) #NS


DunnettTest(n_sub$bg_bl_adj, n_sub$treatment) #NS







##subsetting by stress even further, making drought, nutrient, and no stress datasets


ag_dr<- subset(ag_dat, stress == "drought")
bg_dr<- subset(bg_dat, stress == "drought")
tb_dr<- subset(tb_dat, stress == "drought")
spad_dr<- subset(spad_dat, stress == "drought")

ag_n<- subset(ag_dat, stress == "nutrient")
bg_n<- subset(bg_dat, stress == "nutrient")
tb_n<- subset(tb_dat, stress == "nutrient")
spad_n<- subset(spad_dat, stress == "nutrient")

ag_w<- subset(ag_dat, stress == "none")
bg_w<- subset(bg_dat, stress == "none")
tb_w<- subset(tb_dat, stress == "none")
spad_w<- subset(spad_dat, stress == "none")

#drought stress models
dr_spad_mod<-lm(spad_bl_adj ~ inoculant + exp_rep, spad_dr)
Anova(dr_spad_mod) # NS p = 0.8

dr_tb_mod<-lm(tot_biomass ~ inoculant + exp_rep, tb_dr)
Anova(dr_tb_mod) # NS p = 0.3

dr_bg_mod<-lm(bg_bl_adj ~ inoculant + exp_rep, bg_dr)
Anova(dr_bg_mod) #NS p - 0.5

dr_ag_mod<-lm(ag_bl_adj ~ inoculant + exp_rep, ag_dr)
Anova(dr_ag_mod) #NS p - 0.5


#nutrient stress models
n_spad_mod<-lm(spad_bl_adj ~ inoculant + exp_rep, spad_n)
Anova(n_spad_mod) # NS p = 0.2

n_tb_mod<-lm(tot_biomass ~ inoculant + exp_rep, tb_n)
Anova(n_tb_mod) # NS p = 0.8

n_bg_mod<-lm(bg_bl_adj ~ inoculant + exp_rep, bg_n)
Anova(n_bg_mod) #NS p - 0.9

n_ag_mod<-lm(ag_bl_adj ~ inoculant + exp_rep, ag_n)
Anova(n_ag_mod) #NS p - 0.7



#no stress models
w_spad_mod<-lm(spad_bl_adj ~ inoculant + exp_rep, spad_w)
Anova(w_spad_mod) # NS p = 0.7

w_tb_mod<-lm(tot_biomass ~ inoculant + exp_rep, tb_w)
Anova(w_tb_mod) # NS p = 0.5

w_bg_mod<-lm(bg_bl_adj ~ inoculant + exp_rep, bg_w)
Anova(w_bg_mod) #NS p - 0.4

w_ag_mod<-lm(ag_bl_adj ~ inoculant + exp_rep, ag_w)
Anova(w_ag_mod) #NS p - 0.3





##visualizing
p9<-ggplot(dr_sub, aes(x=inoculant, y = spad, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

p10<-ggplot(n_sub, aes(x=inoculant, y = spad, fill = stress)) +
  geom_boxplot() +
  theme_classic() 


p11<-ggplot(dr_sub, aes(x=inoculant, y = tot_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

p12<-ggplot(n_sub, aes(x=inoculant, y = tot_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 


p13<-ggplot(dr_sub, aes(x=inoculant, y = ag_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

p14<-ggplot(n_sub, aes(x=inoculant, y = ag_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 


p15<-ggplot(dr_sub, aes(x=inoculant, y = bg_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 

p16<-ggplot(n_sub, aes(x=inoculant, y = bg_biomass, fill = stress)) +
  geom_boxplot() +
  theme_classic() 





##subsetting the data with just lalrise and water
lalrise<-subset(spad_dat, inoculant == "Lalrise")
water<-subset(spad_dat, inoculant == "a1water")

sub_dat<-rbind(lalrise, water)
sub_dat$treatment<-as.factor(sub_dat$treatment)

subs_mod<-lm(spad ~ stress*inoculant + exp_rep, data = sub_dat)

Anova(subs_mod) #stress:inoculant p = 0.033
tuk<-HSD.test(subs_mod, "treatment")


'''
                           spad groups
a1None-sterile water   39.00909      a
None-Commercial        36.84545      a
Drought-sterile water  35.78333     ab
Drought-Commercial     30.83333    abc
Nutrient-Commercial    28.62143     bc
Nutrient-sterile water 21.24167      c
'''

lsd<-LSD.test(subs_mod, "treatment")

'''
                           spad groups
a1None-sterile water   39.00909      a
None-Commercial        36.84545      a
Drought-sterile water  35.78333      a
Drought-Commercial     30.83333     ab
Nutrient-Commercial    28.62143      b
Nutrient-sterile water 21.24167      c
'''



sub_spad<-ggplot(sub_dat, aes(x=stress, y = spad, fill = inoculant)) +
  geom_boxplot() +
  ylim(0,55) +
  theme_classic() +
  scale_x_discrete(labels = c("Drought stress", "No stress", "Nutrient stress")) +
  xlab("") +
  ylab("SPAD") +
  theme(axis.title.y = element_text(size = 20, color = "black")) +
  theme(axis.text.x = element_text(size =16, color = "black")) +
  theme(axis.text.y = element_text(size =16, color = "black")) +
  theme(legend.text = element_text(size =18)) +
  theme(legend.title = element_text(size =18, face = "bold")) +
  scale_fill_discrete(labels = c("Sterile water", "Lalrise Start SC")) +
  annotate("text", x = 0.81, y = 45, label = "a", size = 6) +
  annotate("text", x = 1.182, y = 51, label = "ab", size = 6) +
  annotate("text", x = 1.81, y = 53, label = "a", size = 6) +
  annotate("text", x = 2.182, y = 51, label = "a", size = 6) +
  annotate("text", x = 2.81, y = 32, label = "c", size = 6) +
  annotate("text", x = 3.182, y = 42, label = "b", size = 6)



  
  
  
  


subt_mod<-lm(tot_biomass ~ stress*inoculant + exp_rep, sub_dat)
Anova(subt_mod) #stress p = 0.003

suba_mod<-lm(ag_biomass ~ stress*inoculant + exp_rep, sub_dat)
Anova(suba_mod) #stress p < 0.001

subb_mod<-lm(bg_biomass ~ stress*inoculant + exp_rep, sub_dat)
Anova(subb_mod) #exp_rep p = 0.02



##subsetting out lalrise treatments and just looking at bacteria vs water
no_lalrise<- subset(ben_dat, inoculant != "Lalrise")

nls_mod<-lm(spad ~ stress*inoculant + exp_rep, no_lalrise)
Anova(nls_mod)

nlt_mod<-lm(tot_biomass ~ stress*inoculant + exp_rep, no_lalrise)
Anova(nlt_mod)

nla_mod<-lm(ag_biomass ~ stress*inoculant + exp_rep, no_lalrise)
Anova(nla_mod)

nlb_mod<-lm(bg_biomass ~ stress*inoculant + exp_rep, no_lalrise)
Anova(nlb_mod)

nl<-ggplot(no_lalrise, aes(x=inoculant, y = spad, fill = stress)) +
  geom_boxplot() +
  theme_classic() 




















