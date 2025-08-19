#-------------------------------------------------------------------------------
#  January 7, 2025
#  Run by sbR -v 4.2.0 
#  Purpose: Association between metabolites and hip fracture among NHS and HPFS 
#  Programmer: Zhiyuan Wu  /udd/n2zwu/Meta_hip
#-------------------------------------------------------------------------------
rm(list = ls())
options(bitmapType='cairo')
setwd('/udd/n2zwu/Meta_hip')

library(dplyr)
library(Biobase)
library(stringr)
library(fgsea)
library(chanmetab)
library(survival)
library(tidyverse)
library(broom)

load("/udd/n2zwu/Meta_hip/input/data.RData")
# summary(d$age)
# table(d$white)
# table(d$pmhc)
# summary(d$sclerostin)

#-------------------------------------------------------------------------------
# (1). Regression for case-control analysis
#-------------------------------------------------------------------------------
module=colnames(e_imputed) # 653 metabolites

results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=clogit(as.formula(paste0("case ~ ", module[i], " + strata(matchid)")), data=d)
  summ=summary(fit)
  results_mv[i, "n"]=nrow(model.frame(fit))
  results_mv[i, "OR"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p"]=summ$coefficients[1, "Pr(>|z|)"]
  
  fit=clogit(as.formula(paste0("case ~ ", module[i], "+ age+ bmi+ vitd25 + sclerostin + alcoi + osteohx + strata(matchid)")), data=d)
  summ=summary(fit)
  results_mv[i, "OR_1"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci_1"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci_1"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p_1"]=summ$coefficients[1, "Pr(>|z|)"]
}

results_mv$p_adj_fdr=p.adjust(results_mv$p_1, method="fdr", n=length(results_mv$p_1))
results_mv_1=results_mv %>% mutate(population='whole')

# hpfs
# table(hp_d$case)
# names(hp_d)
hp_d=d %>% filter(cohort=='hpfs')
module=colnames(e_hp_imputed[,1:360])

results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=clogit(as.formula(paste0("case ~ ", module[i], "+ strata(matchid)")), data=hp_d)
  summ=summary(fit)
  results_mv[i, "n"]=nrow(model.frame(fit))
  results_mv[i, "OR"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p"]=summ$coefficients[1, "Pr(>|z|)"]
  
  fit=clogit(as.formula(paste0("case ~ ", module[i], "+ age+ bmi+ vitd25 + sclerostin + alcoi + osteohx + strata(matchid)")), data=hp_d)
  summ=summary(fit)
  results_mv[i, "OR_1"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci_1"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci_1"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p_1"]=summ$coefficients[1, "Pr(>|z|)"]
}

results_mv$p_adj_fdr=p.adjust(results_mv$p_1, method="fdr", n=length(results_mv$p_1))
results_mv_2=results_mv %>% mutate(population='hp')


# n1
n1_d=d %>% filter(cohort=='nhs1')
module=colnames(e_n1_imputed[,1:518])

results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=clogit(as.formula(paste0("case ~ ", module[i], "+ strata(matchid)")), data=n1_d)
  summ=summary(fit)
  results_mv[i, "n"]=nrow(model.frame(fit))
  results_mv[i, "OR"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p"]=summ$coefficients[1, "Pr(>|z|)"]
  
  fit=clogit(as.formula(paste0("case ~ ", module[i], " + age+ bmi+ vitd25 + sclerostin + alcoi + osteohx + strata(matchid)")), data=n1_d)
  summ=summary(fit)
  results_mv[i, "OR_1"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci_1"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci_1"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p_1"]=summ$coefficients[1, "Pr(>|z|)"]
}

results_mv$p_adj_fdr=p.adjust(results_mv$p_1, method="fdr", n=length(results_mv$p_1))
results_mv_3=results_mv %>% mutate(population='n1')


# bld_dx_hip <= 120
# table(d$matchid_above_median_index, d$osteohx)
# print(naniar::miss_var_summary(d))
# summary(d$sclerostin)

module=colnames(e_imputed)  # 653 metabolites
results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=clogit(as.formula(paste0("case ~ ", module[i], " + strata(matchid)")), 
             data=d[which(d$matchid_above_median_index=='below'),])
  summ=summary(fit)
  results_mv[i, "n"]=nrow(model.frame(fit))
  results_mv[i, "OR"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p"]=summ$coefficients[1, "Pr(>|z|)"]
  
  fit=clogit(as.formula(paste0("case ~ ", module[i], "+  strata(matchid)")), 
             data=d[which(d$matchid_above_median_index=='below'),])
  summ=summary(fit)
  results_mv[i, "OR_1"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci_1"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci_1"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p_1"]=summ$coefficients[1, "Pr(>|z|)"]
}

results_mv$p_adj_fdr=p.adjust(results_mv$p_1, method="fdr", n=length(results_mv$p_1))
results_mv_4=results_mv %>% mutate(population='below')

# bld_dx_hip > 120
module=colnames(e_imputed) # 653 metabolites

results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=clogit(as.formula(paste0("case ~ ", module[i], " + strata(matchid)")), 
             data=d[which(d$matchid_above_median_index=='above'),])
  summ=summary(fit)
  results_mv[i, "n"]=nrow(model.frame(fit))
  results_mv[i, "OR"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p"]=summ$coefficients[1, "Pr(>|z|)"]
  
  fit=clogit(as.formula(paste0("case ~ ", module[i], "+ strata(matchid)")), 
             data=d[which(d$matchid_above_median_index=='above'),])
  summ=summary(fit)
  results_mv[i, "OR_1"]=round(summ$conf.int[1,1],4)
  results_mv[i, "lci_1"]=round(summ$conf.int[1,3],4)
  results_mv[i, "uci_1"]=round(summ$conf.int[1,4],4)
  results_mv[i, "p_1"]=summ$coefficients[1, "Pr(>|z|)"]
}

results_mv$p_adj_fdr=p.adjust(results_mv$p_1, method="fdr", n=length(results_mv$p_1))
results_mv_5=results_mv %>% mutate(population='above')

# combine regression results
results_mv_main<-rbind(results_mv_1,results_mv_2,results_mv_3,results_mv_4,results_mv_5)
results_mv_main<-merge(results_mv_main, f_match, by = 'met',all.x = TRUE, sort = TRUE)
# table(results_mv_main$population, results_mv_main$p_1<0.05)


# "X100009004"  "X100009051" validated in cohort analysis
# met_sig_whole<-results_mv_main %>% filter(population=='whole' & p_1<0.05) %>% pull(met) # 75
# met_sig_hp<-results_mv_main %>% filter(population=='hp' & p_1<0.05) %>% pull(met) # 36
# met_sig_n1<-results_mv_main %>% filter(population=='n1' & p_1<0.05) %>% pull(met) # 56
# intersect(met_sig_hp, met_sig_n1)  # 6 -- "X52463" "X52629" "X53197" "X53199" "X57338" "X64896"

# met_sig_hp<-results_mv_main %>% filter(population=='below' & p_1<0.05) %>% pull(met) # 67
# met_sig_n1<-results_mv_main %>% filter(population=='above' & p_1<0.05) %>% pull(met) # 48
# intersect(met_sig_hp, met_sig_n1)  # 8 -- "X34437" "X52629" "X52630" "X53199" "X554"  "X57338" "X57721" "X64896"

# results_mv_main %>% filter(population=='n1' & met %in% c("X52463", "X52629", "X53197", "X53199", "X57338", "X64896")) %>% head()


#-------------------------------------------------------------------------------
# (2). Regression for cohort analysis: two metabolite class
#-------------------------------------------------------------------------------

#table(e_dat_network$cohort,e_dat_network$case)
#128/(128+2388): 5.1%
#531/(531+7145): 6.9%
#27/(27+3383)

# 1 -- whole
# names(e_dat_network)
# table(e_dat_network$cohort)
# table(e_dat_network$caco)
# table(e_dat_network$menopmh)
# table(e_dat_network$cohort, e_dat_network$season)
# summary(e_dat_network$time)
module=colnames(e_dat_network[,2:30])

results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i], "+strata(cohort)")), 
            data=e_dat_network[which(e_dat_network$cohort!='nhs2'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i]," + fast+ season+ ageyr + pmh + bmi + osteohx + alco + caco+ endpoint + strata(cohort)")), 
            data=e_dat_network[which(e_dat_network$cohort!='nhs2'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_1<-results_mv %>% mutate(population='cohort-whole', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                         hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                         T ~ 'insignificant'))
results_cox_1$p_adj_fdr=p.adjust(results_cox_1$p_adj, method="fdr", n=length(results_cox_1$p_adj))


# 2 -- male
# print(naniar::miss_var_summary(e_dat_network[which(e_dat_network$cohort=='hpfs'),]))
module=colnames(e_dat_network[,2:30])
elements_to_remove <- c("HMDB0007973", "HMDB0008923","HMDB0008991", "HMDB0008993")
module <- module[!module %in% elements_to_remove]

results_mv=data.frame(met=module[!module %in% elements_to_remove])
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i])), 
            data=e_dat_network[which(e_dat_network$cohort=='hpfs'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i],"+fast+ season+ageyr + bmi + alco+ osteohx + caco + endpoint")), 
            data=e_dat_network[which(e_dat_network$cohort=='hpfs'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_2<-results_mv %>% mutate(population='hpfs', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                         hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                         T ~ 'insignificant'))
results_cox_2$p_adj_fdr=p.adjust(results_cox_2$p_adj, method="fdr", n=length(results_cox_2$p_adj))


# 3 -- female
# print(naniar::miss_var_summary(e_dat_network[which(e_dat_network$cohort=='nhs1'),]))

module=colnames(e_dat_network[,2:30])
results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i])), 
            data=e_dat_network[which(e_dat_network$cohort=='nhs1'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i],"+fast+ season+ ageyr + bmi + alco+ osteohx + pmh + caco + endpoint")), 
            data=e_dat_network[which(e_dat_network$cohort=='nhs1'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_3<-results_mv %>% mutate(population='nhs1', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                        hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                        T ~ 'insignificant'))
results_cox_3$p_adj_fdr=p.adjust(results_cox_3$p_adj, method="fdr", n=length(results_cox_3$p_adj))


# 4 -- cases within 10 years
module=colnames(e_dat_network[,2:30])
results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time_10,case_10) ~ ", module[i], "+strata(cohort)")), 
            data=e_dat_network[which(e_dat_network$cohort!='nhs2'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time_10,case_10) ~ ", module[i]," + fast+ season+ ageyr + bmi + alco+ osteohx + caco+ osteohx +pmh + endpoint+strata(cohort)")), 
            data=e_dat_network[which(e_dat_network$cohort!='nhs2'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_4<-results_mv %>% mutate(population='within10', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                                hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                                T ~ 'insignificant'))
results_cox_4$p_adj_fdr=p.adjust(results_cox_4$p_adj, method="fdr", n=length(results_cox_4$p_adj))


# 5 -- male
# print(naniar::miss_var_summary(e_dat_network[which(e_dat_network$cohort=='hpfs'),]))
module=colnames(e_dat_network[,2:30])
elements_to_remove <- c("HMDB0007973", "HMDB0008923","HMDB0008991", "HMDB0008993")
module <- module[!module %in% elements_to_remove]

results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time_10,case_10) ~ ", module[i])), 
            data=e_dat_network[which(e_dat_network$cohort=='hpfs'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time_10,case_10) ~ ", module[i]," + fast+ season+ ageyr + bmi + alco+ osteohx + caco + endpoint")), 
            data=e_dat_network[which(e_dat_network$cohort=='hpfs'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_5<-results_mv %>% mutate(population='hpfs-within10', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                            hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                            T ~ 'insignificant'))
results_cox_5$p_adj_fdr=p.adjust(results_cox_5$p_adj, method="fdr", n=length(results_cox_5$p_adj))


# 6 -- female
module=colnames(e_dat_network[,2:30])
results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time_10,case_10) ~ ", module[i])), 
            data=e_dat_network[which(e_dat_network$cohort=='nhs1'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time_10,case_10) ~ ", module[i]," + fast+ season+ ageyr + bmi + alco+ osteohx + caco +pmh + endpoint")), 
            data=e_dat_network[which(e_dat_network$cohort=='nhs1'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_6<-results_mv %>% mutate(population='nhs-within10', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                            hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                            T ~ 'insignificant'))
results_cox_6$p_adj_fdr=p.adjust(results_cox_6$p_adj, method="fdr", n=length(results_cox_6$p_adj))


# 7 -- cases exceeding 10 years - whole
module=colnames(e_dat_network[,2:30])
results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case_10_above) ~ ", module[i], "+strata(cohort)")), 
            data=e_dat_network[which(e_dat_network$cohort!='nhs2'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case_10_above) ~ ", module[i]," + fast+ season+ ageyr + bmi + osteohx + alco + caco +pmh + endpoint+strata(cohort)")), 
            data=e_dat_network[which(e_dat_network$cohort!='nhs2'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_7<-results_mv %>% mutate(population='exceed10', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                            hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                            T ~ 'insignificant'))
results_cox_7$p_adj_fdr=p.adjust(results_cox_7$p_adj, method="fdr", n=length(results_cox_7$p_adj))

# 8 -- cases exceeding 10 years - HPFS
module=colnames(e_dat_network[,2:30])
elements_to_remove <- c("HMDB0007973", "HMDB0008923","HMDB0008991", "HMDB0008993")
module <- module[!module %in% elements_to_remove]

results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case_10_above) ~ ", module[i])), 
            data=e_dat_network[which(e_dat_network$cohort=='hpfs'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case_10_above) ~ ", module[i]," + fast+ season+ ageyr + bmi + osteohx+ alco + caco + endpoint")), 
            data=e_dat_network[which(e_dat_network$cohort=='hpfs'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_8<-results_mv %>% mutate(population='exceed10-hpfs', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                            hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                            T ~ 'insignificant'))
results_cox_8$p_adj_fdr=p.adjust(results_cox_8$p_adj, method="fdr", n=length(results_cox_8$p_adj))


# 9 -- cases exceeding 10 years -nhs
module=colnames(e_dat_network[,2:30])
results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case_10_above) ~ ", module[i])), 
            data=e_dat_network[which(e_dat_network$cohort=='nhs1'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case_10_above) ~ ", module[i]," + fast+ season+ ageyr+ osteohx + bmi + alco + caco +pmh + endpoint")), 
            data=e_dat_network[which(e_dat_network$cohort=='nhs1'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_9<-results_mv %>% mutate(population='exceed10-nhs', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                            hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                            T ~ 'insignificant'))
results_cox_9$p_adj_fdr=p.adjust(results_cox_9$p_adj, method="fdr", n=length(results_cox_9$p_adj))


# 10 -- among controls
# table(e_dat_network$caco)
module=colnames(e_dat_network[,2:30])
results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i], "+strata(cohort)")), 
            data=e_dat_network[which(e_dat_network$cohort!='nhs2' & e_dat_network$caco!='case'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i]," + fast+ season+ ageyr + bmi+ osteohx + alco + caco +pmh + endpoint+strata(cohort)")), 
            data=e_dat_network[which(e_dat_network$cohort!='nhs2' & e_dat_network$caco!='case'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_10<-results_mv %>% mutate(population='controls', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                            hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                            T ~ 'insignificant'))
results_cox_10$p_adj_fdr=p.adjust(results_cox_10$p_adj, method="fdr", n=length(results_cox_10$p_adj))


# 11 -- hpfs
module=colnames(e_dat_network[,2:30])
elements_to_remove <- c("HMDB0007973", "HMDB0008923","HMDB0008991", "HMDB0008993")
module <- module[!module %in% elements_to_remove]

results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i])), 
            data=e_dat_network[which(e_dat_network$cohort=='hpfs' & e_dat_network$caco!='case'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i]," + fast+ season+ ageyr + bmi + osteohx + alco + caco + endpoint")), 
            data=e_dat_network[which(e_dat_network$cohort=='hpfs' & e_dat_network$caco!='case'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_11<-results_mv %>% mutate(population='controls-hpfs', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                             hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                             T ~ 'insignificant'))
results_cox_11$p_adj_fdr=p.adjust(results_cox_11$p_adj, method="fdr", n=length(results_cox_11$p_adj))

# 12 -- nhs
# table(e_dat_network$cohort)
module=colnames(e_dat_network[,2:30])
results_mv=data.frame(met=module)
for (i in 1:length(module)) {
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i])), 
            data=e_dat_network[which(e_dat_network$cohort=='nhs1' & e_dat_network$caco!='case'),])
  summ=summary(fit)
  results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll"]=summ$conf.int[1,3]
  results_mv[i, "ul"]=summ$conf.int[1,4]
  results_mv[i, "p"]=summ$coefficients[1,5]
  results_mv[i, "n"]=summ$n
  
  fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i]," + fast+ season+ ageyr+ osteohx + bmi + alco + caco +pmh + endpoint")), 
            data=e_dat_network[which(e_dat_network$cohort=='nhs1' & e_dat_network$caco!='case'),])
  summ=summary(fit)
  results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
  results_mv[i, "ll_adj"]=summ$conf.int[1,3]
  results_mv[i, "ul_adj"]=summ$conf.int[1,4]
  results_mv[i, "p_adj"]=summ$coefficients[1,5]
}

results_cox_12<-results_mv %>% mutate(population='controls-nhs', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
                                                                             hr_adj<1 & p_adj<0.05  ~ 'negative',
                                                                             T ~ 'insignificant'))
results_cox_12$p_adj_fdr=p.adjust(results_cox_12$p_adj, method="fdr", n=length(results_cox_12$p_adj))



results_cox<-rbind(results_cox_1,results_cox_2,results_cox_3,results_cox_4,results_cox_5,results_cox_6,
                   results_cox_7,results_cox_8,results_cox_9,results_cox_10,results_cox_11,results_cox_12)
results_cox<-merge(results_cox, f_blood_match, by = 'met',all.x = TRUE, sort = TRUE)
# table(results_cox$population,results_cox$group)
# results_cox %>% filter(met=='HMDB0008047') %>% head(n=20)
# table(e_dat_network$caco)

# for compare between cox and case-control
# names(results_mv_main)
# names(results_cox)

results_mv_main_match=results_mv_main %>% filter(population=='whole') %>% 
  mutate(rr=OR_1, pp=p_1, ppfdr=p_adj_fdr, met=HMDB) %>% select(met, population, rr, pp, ppfdr)
results_cox_match=results_cox %>% filter(population %in% c('cohort-whole', 'hpfs', 'nhs1')) %>% 
  mutate(rr=hr_adj, pp=p_adj, ppfdr=p_adj_fdr) %>% select(met, population, rr, pp, ppfdr, metabolite_name, class_broad)

results_cox_match_met=results_cox_match %>% pull(met)
results_mv_main_match=results_mv_main_match %>% filter(met %in% results_cox_match_met)

# results_cox_match %>% filter(met %in% c('HMDB0008047','HMDB0009069','HMDB0008923')) %>% head()
results_case_cox<-bind_rows(results_mv_main_match, results_cox_match)
results_case_cox <- results_case_cox %>%
  group_by(met) %>%
  mutate(metabolite_name = ifelse(is.na(metabolite_name), 
                           first(na.omit(metabolite_name)), 
                           metabolite_name),
         class_broad = ifelse(is.na(class_broad), 
                                  first(na.omit(class_broad)), 
                              class_broad)) %>%
  ungroup()



#-------------------------------------------------------------------------------
# (3). Regression for cohort analysis: individual metabolite -- remove
#-------------------------------------------------------------------------------
# whole
# names(e_dat_individual)
# table(e_dat_individual$menopmh)
# print(naniar::miss_var_summary(e_dat_individual))
# module=colnames(e_dat_individual[,2:46])

# results_mv=data.frame(met=module)
# for (i in 1:length(module)) {
#   fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i])), 
#             data=e_dat_individual[which(e_dat_individual$cohort!='nhs2'),])
#   summ=summary(fit)
#   results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
#   results_mv[i, "ll"]=summ$conf.int[1,3]
#   results_mv[i, "ul"]=summ$conf.int[1,4]
#   results_mv[i, "p"]=summ$coefficients[1,5]
#   results_mv[i, "n"]=summ$n
  
#   fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i],"+fast+ season + ageyr + bmi + alco+ endpoint + caco")), 
#             data=e_dat_individual[which(e_dat_individual$cohort!='nhs2'),])
#   summ=summary(fit)
#   results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
#   results_mv[i, "ll_adj"]=summ$conf.int[1,3]
#   results_mv[i, "ul_adj"]=summ$conf.int[1,4]
#   results_mv[i, "p_adj"]=summ$coefficients[1,5]
# }

# results_cox_1<-results_mv %>% mutate(population='whole', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
#                                                                          hr_adj<1 & p_adj<0.05  ~ 'negative',
#                                                                          T ~ 'insignificant'))
# results_cox_1$p_adj_fdr=p.adjust(results_cox_1$p_adj, method="fdr", n=length(results_cox_1$p_adj))

# male
# print(naniar::miss_var_summary(e_dat_individual[which(e_dat_individual$cohort=='hpfs'),]),n=20)

# module=colnames(e_dat_individual[,2:46])
# elements_to_remove <- c("HMDB0000259", "HMDB0008923", "HMDB0010379", "HMDB0010404", "HMDB0012104", 
#                         "HMDB0000033", "HMDB0011745", "HMDB0000192","HMDB0000300",'HMDB0000448')
# module <- module[!module %in% elements_to_remove]

# results_mv=data.frame(met=module[!module %in% elements_to_remove])
# for (i in 1:length(module)) {
#   fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i])), 
#             data=e_dat_individual[which(e_dat_individual$cohort=='hpfs'),])
#   summ=summary(fit)
#   results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
#   results_mv[i, "ll"]=summ$conf.int[1,3]
#   results_mv[i, "ul"]=summ$conf.int[1,4]
#   results_mv[i, "p"]=summ$coefficients[1,5]
#   results_mv[i, "n"]=summ$n
  
#   fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i],"+fast+ season+ ageyr + bmi + alco + endpoint+ caco")), 
#             data=e_dat_individual[which(e_dat_individual$cohort=='hpfs'),])
#   summ=summary(fit)
#   results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
#   results_mv[i, "ll_adj"]=summ$conf.int[1,3]
#   results_mv[i, "ul_adj"]=summ$conf.int[1,4]
#   results_mv[i, "p_adj"]=summ$coefficients[1,5]
# }

# results_cox_2<-results_mv %>% mutate(population='hpfs', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
#                                                                         hr_adj<1 & p_adj<0.05  ~ 'negative',
#                                                                         T ~ 'insignificant'))
# results_cox_2$p_adj_fdr=p.adjust(results_cox_2$p_adj, method="fdr", n=length(results_cox_2$p_adj))

# female
# print(naniar::miss_var_summary(e_dat_individual[which(e_dat_individual$cohort=='nhs1'),]),n=20)

# module=colnames(e_dat_individual[,2:46])
# elements_to_remove <- c("HMDB0010167")
# module <- module[!module %in% elements_to_remove]

# results_mv=data.frame(met=module[!module %in% elements_to_remove])
# for (i in 1:length(module)) {
#   fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i])), 
#             data=e_dat_individual[which(e_dat_individual$cohort=='nhs1'),])
#   summ=summary(fit)
#   results_mv[i, "hr"]=round(summ$conf.int[1,1],4)
#   results_mv[i, "ll"]=summ$conf.int[1,3]
#   results_mv[i, "ul"]=summ$conf.int[1,4]
#   results_mv[i, "p"]=summ$coefficients[1,5]
#   results_mv[i, "n"]=summ$n
  
#   fit=coxph(as.formula(paste0("Surv(time,case) ~ ", module[i],"+fast+ season+ ageyr + bmi + alco+menopmh+ endpoint+ caco")), 
#             data=e_dat_individual[which(e_dat_individual$cohort=='nhs1'),])
#   summ=summary(fit)
#   results_mv[i, "hr_adj"]=round(summ$conf.int[1,1],4)
#   results_mv[i, "ll_adj"]=summ$conf.int[1,3]
#   results_mv[i, "ul_adj"]=summ$conf.int[1,4]
#   results_mv[i, "p_adj"]=summ$coefficients[1,5]
# }

# results_cox_3<-results_mv %>% mutate(population='nhs1', group=case_when(hr_adj>1 & p_adj<0.05  ~ 'positive', 
#                                                                         hr_adj<1 & p_adj<0.05  ~ 'negative',
#                                                                         T ~ 'insignificant'))
# results_cox_3$p_adj_fdr=p.adjust(results_cox_3$p_adj, method="fdr", n=length(results_cox_3$p_adj))

# results_cox_metabolite<-rbind(results_cox_1,results_cox_2,results_cox_3)
# results_cox_metabolite<-merge(results_cox_metabolite, f_blood_match, by = 'met',all.x = TRUE, sort = TRUE)

# names(results_mv_main)
# columns_to_check_sig=results_mv_main %>% filter(p_1<0.05) %>% pull(HMDB)
# table(results_cox_metabolite$population, results_cox_metabolite$group)
# head(results_cox_metabolite[which(results_cox_metabolite$p_adj<0.05 & results_cox_metabolite$met %in% columns_to_check_sig),],100)

# median(e_dat_network$time/12)
# median(e_dat_network[which(e_dat_network$cohort=='hpfs'),]$time/12)
# median(e_dat_network[which(e_dat_network$cohort=='nhs1'),]$time/12)
# median(e_dat_network[which(e_dat_network$case==1),]$time/12)

# ggplot(e_dat_network[which(e_dat_network$case==1),], aes(x=time/12, col=as.factor(cohort)),) + 
# geom_density() + geom_vline(aes(xintercept=20),linetype="dashed", size=1)


#-------------------------------------------------------------------------------
# (4). Data output
#-------------------------------------------------------------------------------
save(results_mv_main, results_cox, results_case_cox, file = "/udd/n2zwu/Meta_hip/input/data_output.RData")



