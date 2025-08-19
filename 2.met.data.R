#-------------------------------------------------------------------------------
#  Run by sbR -v 4.2.0 
#  January 7, 2025
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
library(lubridate)

#-------------------------------------------------------------------------------
# (1). Data readin for case-control analysis
#-------------------------------------------------------------------------------
# extract PM pilot metabolites ####
load("/proj/pohtrs/pohtr0a/All_Cohorts_Pilots/metabolon/P751/p751.results.RData")
result_n1_pm<-results
f_n1_pm=fData(result_n1_pm)
# rownames(f_n1_pm)
# head(f_n1_pm$comp_id)
f_n1_pm_filter<-f_n1_pm[,c("comp_id","pm_pass")]

load("/proj/pohtrs/pohtr0a/All_Cohorts_Pilots/metabolon/P787/p787.results.RData")
result_hp_pm<-results
f_hp_pm=fData(result_hp_pm)
f_hp_pm_filter<-f_hp_pm[,c("comp_id","pm_pass")]

# Metabolon Match Data ####
# metabolon_pm_filter=rbind(f_n1_pm_filter, f_hp_pm_filter)
# table(metabolon_pm_filter$pm_pass)
# metabolon_pm_filter <- metabolon_pm_filter %>% group_by(comp_id) %>% mutate(pm_pass_max=max(pm_pass)) %>% ungroup() %>% filter(!duplicated(comp_id))
# setwd("~/lowcarb_chd/file")
# write.csv(metabolon_pm_filter, file="metabolon_pm_filter.csv") 

# metabolome and covariate data ####
# HP
load("/proj/hpblds/hpbld00/endpoints/hipfrac/LABCODES/Lab60110/results/lab60110.results.RData")
hipcov_hp<-read.csv("/udd/n2zwu/Meta_hip/input/hip_hp_cov.csv", header = TRUE)
result_hp<-results

p_hp=pData(result_hp)
f_hp=fData(result_hp)
e_hp=data.frame(exprs(result_hp))
# names(f_hp)
# head(f_hp)
hp_filter_COMP_ID<-dput(f_hp_pm_filter[which(f_hp_pm_filter$pm_pass==1),]$comp_id)
f_hp <- subset(f_hp, !(comp_id %in% hp_filter_COMP_ID)) # 1275 change to 489

p_hp$id <- str_sub(p_hp$id, start = 1, end = 6) 
hipcov_hp$id <- str_sub(hipcov_hp$id, start = 1, end = 6) 
p_hp <- subset(p_hp, (p_hp$id %in% hipcov_hp$id)) # 592 samples

colnames(e_hp)<-str_sub(colnames(e_hp), start = 2, end = 7) 
rownames(e_hp)<-gsub("X", "", rownames(e_hp))

e_hp<-e_hp[rownames(e_hp) %in% f_hp$comp_id,] # select metabolite passing PM
e_hp<-e_hp[,colnames(e_hp) %in% p_hp$id] # filter out QC sample
e_hp<-data.frame(t(e_hp)) # 592::489

# keep Missing <=70%
na_counts_hp <- colMeans(is.na(e_hp))
e_hp_1 <- e_hp[, na_counts_hp <= 0.2] # 489 to 249 metabolites remain
e_hp_2 <- e_hp[, na_counts_hp <= 0.7 & na_counts_hp > 0.2] # 111 metabolites remain

# N1
load("/proj/nhblds/nhbld00/endpoints/fractures/LABCODES/lab10280/lab10280.results.RData")
hipcov_n1<-read.csv("/udd/n2zwu/Meta_hip/input/hip_n1_cov.csv", header = TRUE)
result_n1<-results

p_n1=pData(result_n1)
f_n1=fData(result_n1)
e_n1=data.frame(exprs(result_n1))

n1_filter_COMP_ID<-dput(f_n1_pm_filter[which(f_n1_pm_filter$pm_pass==1),]$comp_id)
f_n1 <- subset(f_n1, !(comp_id %in% n1_filter_COMP_ID)) # 1291 change to 710

p_n1$id <- str_sub(p_n1$id, start = 1, end = 6) 
hipcov_n1$id <- str_sub(hipcov_n1$id, start = 1, end = 6) 
p_n1 <- subset(p_n1, (p_n1$id %in% hipcov_n1$id)) # 642 samples

colnames(e_n1)<-str_sub(colnames(e_n1), start = 2, end = 7) 
rownames(e_n1)<-gsub("X", "", rownames(e_n1))

e_n1<-e_n1[rownames(e_n1) %in% f_n1$comp_id,] # select metabolite passing PM
e_n1<-e_n1[,colnames(e_n1) %in% p_n1$id] # filter out QC sample
e_n1<-data.frame(t(e_n1)) # 642::710

# Missing <=70%
na_counts_n1 <- colMeans(is.na(e_n1))
e_n1_1 <- e_n1[, na_counts_n1 <= 0.2] # 710 to 404 metabolites remain
e_n1_2 <- e_n1[, na_counts_n1 <= 0.7 & na_counts_n1 > 0.2] # 114 metabolites remain

# f_match file
f_hp<-f_hp %>% mutate(met=rownames(f_hp), 
                      SUPER_PATHWAY=super_pathway,
                      SUB_PATHWAY=sub_pathway,
                      METABOLITE=metabolite,
                      CHEMICAL_NAME=chemical_name,
                      PLOT_NAME=plot_name,
                      HMDB=hmdb,
                      KEGG=kegg,
                      PLATFORM=platform,
                      meanCV=meancv) %>% select(met,SUPER_PATHWAY,SUB_PATHWAY,METABOLITE,CHEMICAL_NAME,PLOT_NAME,
                                                HMDB,KEGG,PLATFORM,meanCV,icc)
f_n1<-f_n1 %>% mutate(met=rownames(f_n1),
                      SUPER_PATHWAY=super_pathway,
                      SUB_PATHWAY=sub_pathway,
                      METABOLITE=metabolite,
                      CHEMICAL_NAME=chemical_name,
                      PLOT_NAME=plot_name,
                      HMDB=hmdb,
                      KEGG=kegg,
                      PLATFORM=platform,
                      meanCV=meancv) %>% select(met,SUPER_PATHWAY,SUB_PATHWAY,METABOLITE,CHEMICAL_NAME,PLOT_NAME,
                                                HMDB,KEGG,PLATFORM,meanCV,icc)
f_match<-rbind(f_hp,f_n1) 
f_match<- f_match[order(f_match$met,f_match$meanCV),]
f_match<-f_match %>% filter(!duplicated(met)) # 877

# Impute using half minimum by cohort
# HPFS
e_hp_imputed<-  e_hp_1 %>% mutate(id=rownames(e_hp_1))
e_hp_imputed <- data.frame(sapply(e_hp_imputed, function(x) ifelse(is.na(x), min(x, na.rm = TRUE) / 2, x))) # met =249
e_hp_imputed <- as.data.frame(lapply(e_hp_imputed, as.numeric))
rownames(e_hp_imputed)<-e_hp_imputed$id
e_hp_imputed<-e_hp_imputed[,-250]
identical(rownames(e_hp_imputed),rownames(e_hp_2))
e_hp_imputed<-cbind(e_hp_imputed,e_hp_2) ## 592::360
e_hp_imputed <-as.data.frame(scale(log(e_hp_imputed)))
# dput(colnames(e_hp_imputed))

# NHS
e_n1_imputed<-  e_n1_1 %>% mutate(id=rownames(e_n1_1))
e_n1_imputed <- data.frame(sapply(e_n1_imputed, function(x) ifelse(is.na(x), min(x, na.rm = TRUE) / 2, x))) # met =404
e_n1_imputed <- as.data.frame(lapply(e_n1_imputed, as.numeric))
rownames(e_n1_imputed)<-e_n1_imputed$id
e_n1_imputed<-e_n1_imputed[,-405]
identical(rownames(e_n1_imputed),rownames(e_n1_2))
e_n1_imputed<-cbind(e_n1_imputed,e_n1_2) ## 642::518
e_n1_imputed <-as.data.frame(scale(log(e_n1_imputed)))
# dput(colnames(e_n1_imputed))

# Source data #
common_column <- intersect(colnames(e_hp_imputed), colnames(e_n1_imputed)) # 225 both measured metabolites
e_imputed<-bind_rows (e_hp_imputed, e_n1_imputed) # 1234::653 metabolites


e_hp_imputed$id<-rownames(e_hp_imputed)
d_hp_imputed<-merge(e_hp_imputed, hipcov_hp, by = 'id',all = FALSE, sort = TRUE)
met=dput(names(d_hp_imputed[,1:361]))

d_hp_imputed_2<-d_hp_imputed %>% select(all_of(met), "case","bxc", "creatinine", "osteocalcin", "inpth", "vitd25","sclerostin","matchid","fast", "race","smkst",
                                        "age", "agec", "bmi","physact","htnhx", "dmhx", "osteohx","chdhx","cancer", "thiaz", "bisphos08",
                                        "dcalci", "scalc","dvitdi", "svitd","caffi", "alcoi","tproti", "aproti", "vproti", "dproti", "ndproti",
                                        "hiptime", "time_bld_dx_months","dtdth")
#table(d_hp_imputed_2$case)
#summary(d_hp_imputed_2[which(d_hp_imputed_2$case==0),]$dtdth)
#summary(d_hp_imputed_2[which(d_hp_imputed_2$case==0),]$hiptime)
#summary(d_hp_imputed_2[which(d_hp_imputed_2$case==0),]$time_bld_dx_months)


e_n1_imputed$id<-rownames(e_n1_imputed)
d_n1_imputed<-merge(e_n1_imputed, hipcov_n1, by = 'id',all = FALSE, sort = TRUE) #
met=dput(names(d_n1_imputed[,1:519]))
d_n1_imputed_2<-d_n1_imputed %>% select(all_of(met),"case","bxc", "calcium", "bicarbon", "creatinine", "osteocalcin", "inpth", "vitd25", "sclerostin","matchid",
                                        "fast12", "race","age", "agec", "bmi","physact","htnhx", "dmhx", "osteohx","chdhx","cancer","fallhx","thiaz", "bisphos02",
                                        "pmhc","smkst","dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti", "aproti", "vproti", "dproti", "ndproti",
                                        "hiptime", "time_bld_dx_months","dtdth")

names(d_hp_imputed_2)[names(d_hp_imputed_2) == "bisphos08"] <- "bisphos"
names(d_n1_imputed_2)[names(d_n1_imputed_2) == "bisphos02"] <- "bisphos"
names(d_n1_imputed_2)[names(d_n1_imputed_2) == "fast12"] <- "fast"
# table(d_hp_imputed_3$fast)
# table(d_n1_imputed_3$fast)
# table(hipcov_n1$pmhc)
d_hp_imputed_3=d_hp_imputed_2 %>%  mutate(bxc = tidyr::replace_na(bxc, median(bxc, na.rm = TRUE)),
                                          sex='male',cohort='hpfs', pmhc='male',
                                          osteocalcin = tidyr::replace_na(osteocalcin, median(osteocalcin, na.rm = TRUE)),
                                          inpth = tidyr::replace_na(inpth, median(inpth, na.rm = TRUE)),
                                          vitd25 = tidyr::replace_na(vitd25, median(vitd25, na.rm = TRUE)),
                                          sclerostin = tidyr::replace_na(sclerostin, median(sclerostin, na.rm = TRUE)),
                                          creatinine  = tidyr::replace_na(creatinine , median(creatinine , na.rm = TRUE)),
                                          bmi = tidyr::replace_na(bmi, median(bmi, na.rm = TRUE)),
                                          physact = tidyr::replace_na(physact, median(physact, na.rm = TRUE)),
                                          dcalci = tidyr::replace_na(dcalci, median(dcalci, na.rm = TRUE)),
                                          scalc = tidyr::replace_na(scalc, median(scalc, na.rm = TRUE)),
                                          dvitdi = tidyr::replace_na(dvitdi, median(dvitdi, na.rm = TRUE)),
                                          svitd = tidyr::replace_na(svitd, median(svitd, na.rm = TRUE)),
                                          caffi  = tidyr::replace_na(caffi , median(caffi , na.rm = TRUE)),
                                          alcoi = tidyr::replace_na(alcoi, median(alcoi, na.rm = TRUE)),
                                          tproti = tidyr::replace_na(tproti, median(tproti, na.rm = TRUE)))

d_n1_imputed_3=d_n1_imputed_2 %>%  mutate(fallhx=case_when(fallhx==2 ~ 1, T~0), # falling history, calcium, bicarbon, pmhc only in n1
                                          sex='female',cohort='nhs1',
                                          race=case_when(race==0 ~ 1, T ~ 0),
                                          fast=case_when(fast=='fasting 8+ hrs' ~ 1, T~0),
                                          agec=case_when(agec=='<60' ~ 1,agec=='60-64' ~ 2,agec=='65-69' ~ 3,agec=='70-75' ~ 4, T~5),
                                          bxc = tidyr::replace_na(bxc, median(bxc, na.rm = TRUE)),
                                          calcium = tidyr::replace_na(calcium, median(calcium, na.rm = TRUE)),
                                          osteocalcin = tidyr::replace_na(osteocalcin, median(osteocalcin, na.rm = TRUE)),
                                          inpth = tidyr::replace_na(inpth, median(inpth, na.rm = TRUE)),
                                          vitd25 = tidyr::replace_na(vitd25, median(vitd25, na.rm = TRUE)),
                                          sclerostin = tidyr::replace_na(sclerostin, median(sclerostin, na.rm = TRUE)),
                                          creatinine  = tidyr::replace_na(creatinine , median(creatinine , na.rm = TRUE)),
                                          bmi = tidyr::replace_na(bmi, median(bmi, na.rm = TRUE)),
                                          physact = tidyr::replace_na(physact, median(physact, na.rm = TRUE)),
                                          dcalci = tidyr::replace_na(dcalci, median(dcalci, na.rm = TRUE)),
                                          scalc = tidyr::replace_na(scalc, median(scalc, na.rm = TRUE)),
                                          dvitdi = tidyr::replace_na(dvitdi, median(dvitdi, na.rm = TRUE)),
                                          svitd = tidyr::replace_na(svitd, median(svitd, na.rm = TRUE)),
                                          caffi  = tidyr::replace_na(caffi , median(caffi , na.rm = TRUE)),
                                          alcoi = tidyr::replace_na(alcoi, median(alcoi, na.rm = TRUE)),
                                          tproti = tidyr::replace_na(tproti, median(tproti, na.rm = TRUE)))
# table(is.na(d_n1_imputed_3$pmhc))


# Combine data 
d_imputed<-bind_rows(d_hp_imputed_3, d_n1_imputed_3) # Chr not accepted
# table(d_imputed$fast)

d<-d_imputed %>% group_by(cohort) %>% mutate(bmic=case_when(bmi<25 ~ 1,bmi>=30 ~ 3,T~2), 
                                             white=case_when(race==1 ~ 1, T~0), 
                                             alcoiq=case_when(alcoi==0 ~ 1,alcoi<15 ~ 2, T~3),
                                             scalcq=case_when(scalc==0 ~ 1,scalc<100 ~ 2,scalc<500 ~ 3, T~4),
                                             svitd=case_when(svitd<0 ~ 0, T~ svitd), # 2 subjects <0
                                             svitdq=case_when(svitd==0 ~ 1,svitd<200 ~ 2,svitd<400 ~ 3, T~4),
                                             dcalciq=statar::xtile(dcalci, n = 5),
                                             dvitdiq=statar::xtile(dvitdi, n = 5),
                                             tprotiq=statar::xtile(tproti, n = 5),
                                             caffiq =statar::xtile(caffi, n = 5),
                                             physactq =statar::xtile(physact, n = 5),
                                             bxcq=statar::xtile(bxc, n = 5),
                                             osteocalcinq=statar::xtile(osteocalcin, n = 5),
                                             inpthq=statar::xtile(inpth, n = 5),
                                             vitd25q=case_when(vitd25<20 ~ 1,vitd25>=30 ~ 3, T~1),
                                             sclerostinq =statar::xtile(sclerostin, n = 5),
                                             pmhc=case_when(is.na(pmhc) ~ 'male', T~pmhc),
                                             smkst=case_when(smkst=='current' ~ 'current', smkst=='past' ~ 'past', T~'never'),
                                             egfr=  case_when(sex=='female' & creatinine<=0.7 ~ 142*((creatinine/0.7)^-0.241)*(0.9938^age)*1.012, # PMID: 34554658
                                                              sex=='female' & creatinine>0.7 ~  142*((creatinine/0.7)^-1.200)*(0.9938^age)*1.012, 
                                                              sex=='male'   & creatinine<=0.9 ~ 142*((creatinine/0.9)^-0.302)*(0.9938^age), 
                                                              sex=='male'   & creatinine>0.9 ~  142*((creatinine/0.9)^-1.200)*(0.9938^age)),
                                             egfrq=case_when(egfr<60 ~ 3, egfr<90 ~ 2, T~1),) %>% ungroup()


# stratify by median time between bld and dx_hip
# d %>% select(case, matchid) %>% arrange(matchid) %>% head()
# d %>% filter(case==1) %>% select(time_bld_dx_months) %>% summary()
matchid_above_median<- d %>% filter(case==1) %>% filter(time_bld_dx_months > 120) %>% pull(matchid) # n=310
d<-d %>% mutate(matchid_above_median_index=case_when(matchid %in% matchid_above_median ~ 'above', T ~ 'below'))
# table(d$matchid_above_median_index)

d<- d %>%  mutate(id=str_sub(id, start = 1, end = 6),
                  newid=case_when(cohort=="nhs1" ~ as.numeric(id)+1000000,
                                  cohort=="nhs2" ~ as.numeric(id)+2000000,
                                  cohort=="hpfs" ~ as.numeric(id)+3000000),
                  pmh=case_when(pmhc=='current PMH use' ~ 1, T ~ 0))

# summary(d$calcium)
# names(d)
# met_pca<-f_match %>% filter(SUB_PATHWAY=='Phosphatidylcholine (PC)') %>% pull(met)
# met_pca<-met_pca[met_pca %in% colnames(e_imputed)]
# d_pca_pc<-d %>% select(all_of(met_pca),case,age,bmi,alcoi,osteohx,vitd25,sclerostin,matchid,cohort)

#met_pca<-f_match %>% filter(SUB_PATHWAY=='Phosphatidylethanolamine (PE)') %>% pull(met)
#met_pca<-met_pca[met_pca %in% colnames(e_imputed)]
#d_pca_pe<-d %>% select(all_of(met_pca),case,age,bmi,alcoi,osteohx,vitd25,sclerostin,matchid,cohort)

#-------------------------------------------------------------------------------
# (2). Data readin for cohort analysis
#-------------------------------------------------------------------------------
# 1 -- death data
data_death<-read.csv("/udd/n2zwu/Meta_hip/input/data_death.csv", header = TRUE)
# summary(data_death$dtdth)
# names(data_death)
data_death<-data_death %>% filter(!duplicated(newid))  %>% select(newid, dtdth)

# 2 -- reported cases data
hfrac.reports_n1=read.fwf('/proj/nhdats/nh_dat_cdx/endpoints/hipfrac/hfrac8220.071621.reports', 
                          widths=c(6,1,-2,1,-2,5,-2,2,-4,2,2,-2,2,-2,5,-2,2,2,-2,4,-2,6),
                          col.names=c("id","cd","cohort","icda","qyr","mdx1","ydx1","conf","nhsicda","mdx","ydx","dxmonth","rc_record_id"))
hfrac.reports_n1<-hfrac.reports_n1 %>% mutate(newid=as.numeric(id)+1000000, hfdtdx=dxmonth) %>% select(newid,hfdtdx)

hfrac.reports_n2=read.fwf('/proj/n2dats/n2_dat_cdx/endpoints/hipfrac/hfrac1319.071921.reports', 
                          widths=c(6,1,-2,1,-2,5,-2,2,-4,2,2,-2,2,-2,5,-2,2,2,-2,4,-2,6),
                          col.names=c("id","cd","cohort","icda","qyr","mdx1","ydx1","conf","nhsicda","mdx","ydx","dxmonth","rc_record_id"))
hfrac.reports_n2<-hfrac.reports_n2 %>% mutate(newid=as.numeric(id)+2000000, hfdtdx=dxmonth) %>% select(newid,hfdtdx)

hfrac.reports_hf=read.table(
  file = '/proj/nhvfxs/nhvfx0a/diet_acid_load/HPFS/hpfs.fracture.dtdx.dat',
  header = FALSE,             # Change to TRUE if the file contains a header row
  sep = " ",                  # Specify the space as the separator
  col.names = c("id", "fracdtdx", "hfdtdx", "wfdtdx"),
  strip.white = TRUE          # Strip leading/trailing whitespace from unquoted character fields
)
hfrac.reports_hf<-hfrac.reports_hf %>% mutate(newid=as.numeric(id)+3000000) %>% select(newid, hfdtdx)

hfrac.reports<-rbind(hfrac.reports_hf,hfrac.reports_n1,hfrac.reports_n2) # 14354

# Filter rows where 'ID' appears more than once
# df_duplicates <- hfrac.reports %>%
#  group_by(newid) %>%
#  filter(n() > 1)
# Print the rows with duplicate IDs
# print(df_duplicates)

hfrac.reports<-hfrac.reports %>% mutate(hfdtdx=as.numeric(hfdtdx)) %>% 
  filter(hfdtdx>0) %>% distinct(newid, .keep_all = TRUE) # keep only case 10245

#d<-merge(d, hfrac.reports, by = 'newid',all.x = TRUE, sort = TRUE)
#summary(d[which(d$case==1 & d$cohort=='nhs1'),]$hfdtdx)
#summary(d[which(d$case==1 & d$cohort=='hpfs'),]$hfdtdx)

# 3 -- osteo history
n1_osteo <- haven::read_sas("/udd/hpeha/htn/nhs1_osteo.sas7bdat")
hp_osteo <- haven::read_sas("/udd/hpeha/htn/hpfs_osteo.sas7bdat")

n1_osteo <- n1_osteo %>% mutate(# osteohx=case_when(ost82e==1|ost84e==1|ost86e==1|ost88e==1|ost90e==1 ~ 1, T ~ 0),
                                osteohx = if_else(rowSums(select(., ost82e:ost90e) == 1, na.rm = TRUE) > 0, 1, 0),
                                newid=as.numeric(id)+1000000) %>% 
  select(newid, osteohx)

hp_osteo <- hp_osteo %>% mutate(osteohx = if_else(rowSums(select(., ost86e:ost94e) == 1, na.rm = TRUE) > 0, 1, 0),
                                newid=as.numeric(id)+3000000) %>% 
  select(newid, osteohx)

# table(n1_osteo$osteohx)
# table(hp_osteo$osteohx)

osteo.reports=rbind(n1_osteo, hp_osteo)
osteo.reports<-osteo.reports %>% filter(!duplicated(newid)) # 57321


# cohort metabolome data ####
all_first <- merge_metab_data(merge_type = "union",
                              # methods = c("C8-pos", "HILIC-pos"),
                              collection_to_use = "first",
                              cohorts = c("nhs1","hpfs"),
                              transformation = "transform_ln_z_score", # default transform_z_score transform_ln_z_score transform_none
                              # impute_cutoff = 0.0,
                              combine_cohorts = TRUE)

p=pData(all_first$expr_set$all_cohorts)
f=fData(all_first$expr_set$all_cohorts)
f_blood_match<-f %>% mutate(met=rownames(f)) %>% select("met","metabolite_name",'biochemical_name',"class_broad","sub_class_metabolon")
e=data.frame(exprs(all_first$expr_set$all_cohorts))
rownames(p) <- paste0("X", rownames(p))


# network data ####
class_names<-c('Phosphatidylcholines','Phosphatidylethanolamines')
met_validate_names=rownames(f[which(f$class_broad %in% class_names),]) # 16+14=30

# met_discovery<-f_match %>% filter(SUB_PATHWAY=='Phosphatidylcholine (PC)'|SUB_PATHWAY=='Phosphatidylethanolamine (PE)') %>% pull(HMDB)
# elements_to_remove <- c(NA, "HMDB0008143,HMDB0008045") # 28 with HMDB
# met_discovery <- met_discovery[!met_discovery %in% elements_to_remove]
# intersect(met_discovery, met_validate_names) # HMDB0008047 HMDB0008924 HMDB0008923 HMDB0009069
# head(f_match[which(f_match$HMDB %in% c("HMDB0008924", "HMDB0008047", "HMDB0008923", "HMDB0009069")),])

e_dat<-e[rownames(e) %in% met_validate_names,]
e_dat<-as.data.frame(t(e_dat))

na_counts <- colMeans(is.na(e_dat))
e_dat <- e_dat[, na_counts <= 0.7] # 29 metabolites remain HMDB0008924 removed


identical(rownames(e_dat), rownames(p)) # should be TRUE
e_dat=cbind(e_dat, p) 
e_dat=e_dat %>%  mutate(id=str_sub(id, start = 1, end = 6),
                        newid=case_when(cohort=="nhs1" ~ as.numeric(id)+1000000,
                                        cohort=="nhs2" ~ as.numeric(id)+2000000,
                                        cohort=="hpfs" ~ as.numeric(id)+3000000))


e_dat<-merge(e_dat, hfrac.reports, by = 'newid',   all.x = TRUE, sort = TRUE)
e_dat<-merge(e_dat, data_death,    by = 'newid',   all.x = TRUE, sort = TRUE)
e_dat<-merge(e_dat, osteo.reports, by = 'newid',   all.x = TRUE, sort = TRUE)


e_dat<-e_dat %>% mutate(hfdtdx=as.numeric(hfdtdx), 
                        hfdtdx=case_when(is.na(hfdtdx) ~ 9999, T ~ hfdtdx)) %>% 
  filter(blddate<hfdtdx) %>% 
  mutate(case=case_when(cohort=='hpfs' & hfdtdx<=1393 ~ 1, # 2016-1
                        cohort=='nhs1' & hfdtdx<=1398 ~ 1, # 2016-6
                        T ~ 0),
         time=case_when(case==1 ~ hfdtdx-blddate, 
                        
                        cohort=='hpfs' & dtdth<1393 ~ dtdth-blddate,
                        cohort=='hpfs' ~ 1393-blddate,
                        cohort=='nhs1' & dtdth<1398 ~ dtdth-blddate,
                        cohort=='nhs1' ~ 1398-blddate),
         
         pmh=case_when(menopmh=='post on PMH' ~ 1, T ~ 0),
         White=case_when(race=='White'~'White', T ~ 'Others'),
         osteohx=case_when(osteohx==1 ~ 1, T ~ 0),
         time_bld_dx_months=case_when(case==1 ~ time, T ~ NA))
# table(e_dat$cohort)


# Convert blddate to Date format
e_dat$date <- as.Date("1900-01-01") %m+% months(e_dat$blddate)
# Extract the month from the date
e_dat$month <- as.numeric(format(e_dat$date, "%m"))
# Create the season variable
e_dat$season <- ifelse(e_dat$month %in% c(12, 1, 2), "Winter",
                 ifelse(e_dat$month %in% c(3, 4, 5), "Spring",
                        ifelse(e_dat$month %in% c(6, 7, 8), "Summer", "Autumn")))

e_dat_network <- e_dat  
# table(e_dat_network$cohort, e_dat_network$osteohx)
# print(naniar::miss_var_summary(e_dat_network))
# 301/12
# filter(rowSums(is.na(select(., all_of(columns_to_check)))) != length(columns_to_check))


# individual data ####
#match_file<-read.csv("/udd/n2zwu/Meta_hip/input/SOLNAS_SOL_batch1_metabolites_matching_to Qi_07302024.csv", header = TRUE)
#names(match_file)[3]<-'HMDB_matched'

#names(f_match)
#names(match_file)

#table(match_file$matched_broad_metabolon)

#f_match_matched<-merge(f_match, match_file[,c(3,11)], by.x = 'METABOLITE',by.y='biochemical', all.x = TRUE, sort = TRUE)
#f_match_matched<-f_match_matched %>% filter(!duplicated(met)) %>% mutate(HMDB_final=case_when(!is.na(HMDB) ~ HMDB, T ~ HMDB_matched),
#                                                                         HMDB_final=sub(",.*", "", HMDB_final))
#table(duplicated(f_match_matched$HMDB_final))
#columns_to_check <-f_match_matched %>% filter(met %in% colnames(e_imputed)) %>%
#  filter(!is.na(HMDB_final) & !duplicated(HMDB_final)) %>% pull(HMDB_final) # 440
# names(f_match)
columns_to_check <-f_match %>% filter(met %in% colnames(e_imputed)) %>% mutate(HMDB=sub(",.*", "", HMDB)) %>% 
  filter(!is.na(HMDB) & !duplicated(HMDB)) %>% pull(HMDB) # 440

e_dat<-e[rownames(e) %in% columns_to_check,]
e_dat<-as.data.frame(t(e_dat))

na_counts <- colMeans(is.na(e_dat))
e_dat <- e_dat[, na_counts <= 0.7] # 45 metabolites remain

identical(rownames(e_dat), rownames(p)) # should be TRUE
e_dat=cbind(e_dat, p) 
e_dat=e_dat %>%  mutate(id=str_sub(id, start = 1, end = 6),
                        newid=case_when(cohort=="nhs1" ~ as.numeric(id)+1000000,
                                        cohort=="nhs2" ~ as.numeric(id)+2000000,
                                        cohort=="hpfs" ~ as.numeric(id)+3000000))

e_dat<-merge(e_dat, hfrac.reports, by = 'newid',all.x = TRUE, sort = TRUE)
e_dat<-merge(e_dat, data_death, by = 'newid',all.x = TRUE, sort = TRUE)

e_dat<-e_dat %>% mutate(hfdtdx=as.numeric(hfdtdx), 
                        hfdtdx=case_when(is.na(hfdtdx) ~ 9999, T ~ hfdtdx)) %>% 
  filter(blddate<hfdtdx) %>% 
  mutate(case=case_when(cohort=='hpfs' & hfdtdx<=1393 ~ 1, 
                        cohort=='nhs1' & hfdtdx<=1398 ~ 1,
                        T ~ 0),
         time=case_when(case==1 ~ hfdtdx-blddate, 
                        
                        cohort=='hpfs' & dtdth<1393 ~ dtdth-blddate,
                        cohort=='hpfs' ~ 1393-blddate,
                        cohort=='nhs1' & dtdth<1398 ~ dtdth-blddate,
                        cohort=='nhs1' ~ 1398-blddate)) # 5331 samples

e_dat$date <- as.Date("1900-01-01") %m+% months(e_dat$blddate)
# Extract the month from the date
e_dat$month <- as.numeric(format(e_dat$date, "%m"))
# Create the season variable
e_dat$season <- ifelse(e_dat$month %in% c(12, 1, 2), "Winter",
                       ifelse(e_dat$month %in% c(3, 4, 5), "Spring",
                              ifelse(e_dat$month %in% c(6, 7, 8), "Summer", "Autumn")))
e_dat_individual <- e_dat

# e_dat_network<- e_dat_network %>% filter(cohort!='nhs2') # 10192
# e_dat_individual<- e_dat_individual %>% filter(cohort!='nhs2') # 10192



# remove overlap case
case_case<-d %>% filter(case==1) %>% pull(newid)

e_dat_network<-e_dat_network %>% mutate(overlap=case_when((newid %in% case_case) & case==1  ~ 'overlapped',
                                                          T ~ 'nonoverlapped')) %>% filter(overlap=='nonoverlapped') %>% 
  mutate(case_10=case_when(case==1 & time<=120 ~ 1, T ~ 0), time_10=case_when(time<=120 ~ time, T ~ 120), 
         case_10_above=case_when(case==1 & time>120 ~ 1, T ~ 0), time_10_above=time) 

e_dat_individual<-e_dat_individual %>% mutate(overlap=case_when((newid %in% case_case) & case==1  ~ 'overlapped',
                                                          T ~ 'nonoverlapped')) %>% filter(overlap=='nonoverlapped') %>% 
  mutate(case_10=case_when(case==1 & time<=120 ~ 1, T ~ 0), time_10=case_when(time<=120 ~ time, T ~ 120), 
         case_10_above=case_when(case==1 & time>120 ~ 1, T ~ 0), time_10_above=time) 



#-------------------------------------------------------------------------------
# (3). Data output
#-------------------------------------------------------------------------------
save(d, e_hp_imputed, e_n1_imputed, e_imputed, 
     e_dat_network, # e_dat_individual,
     f_match, f_blood_match, file = "/udd/n2zwu/Meta_hip/input/data.RData")
