
#-------------------------------------------------------------------------------
#  January 7, 2025
#  Run by sbR -v 4.2.0 
#  Purpose: Association between metabolites and hip fracture among NHS and HPFS 
#  Programmer: Zhiyuan Wu  /udd/n2zwu/Meta_hip
#-------------------------------------------------------------------------------

rm(list = ls())
options(bitmapType='cairo')
setwd("~/Meta_hip/output")

library(dplyr)
library(Biobase)
library(stringr)
library(fgsea)
library(chanmetab)
library(survival)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(ggnewscale)
library(ggtreeExtra)
library(MetBrewer)
library(ggsci)
library(viridis)
library(forcats)
library(ggrepel)
library(tableone)

load("/udd/n2zwu/Meta_hip/input/data.RData")
load("/udd/n2zwu/Meta_hip/input/data_output.RData")
#-------------------------------------------------------------------------------
# (1). Table 1
#-------------------------------------------------------------------------------
table1<- CreateTableOne(vars = c( "sex", "white", "smkst", "age", "bmi", "physact", "time_bld_dx_months",
                                  "htnhx", "dmhx", "osteohx", "chdhx",
                                  "dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti",
                                  "bxc", "osteocalcin", "inpth", "vitd25","sclerostin", "thiaz", "bisphos"),
                        factorVars = c( "white", "smkst","htnhx", "dmhx", "osteohx", "chdhx", 
                                        "thiaz", "bisphos","sex"), data = d) 
table1<-print(table1,nonnormal = c("time_bld_dx_months"), cramVars = "hepato", 
              quote = FALSE, noSpaces = TRUE, showAllLevels = FALSE)

table2<- CreateTableOne(vars = c( "sex", "white", "smkst", "age", "bmi","physact", "time_bld_dx_months",
                                  "htnhx", "dmhx", "osteohx", "chdhx", 
                                  "dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti",
                                  "bxc", "osteocalcin", "inpth", "vitd25","sclerostin", "thiaz", "bisphos"),
                        factorVars = c( "white", "smkst",'pmhc',"htnhx", "dmhx", "osteohx", "chdhx", 
                                        "thiaz", "bisphos","sex"),data = d, strata = "case")
table2<-print(table2,nonnormal = c("time_bld_dx_months"), cramVars = "hepato", 
              quote = FALSE, noSpaces = TRUE, showAllLevels = FALSE)
t1=cbind(table1, table2)
write.csv(t1, "table_1.csv")

#-------------------------------------------------------------------------------
# (2). Table S1
#-------------------------------------------------------------------------------
table1<- CreateTableOne(vars = c( "white", "smkst", "age", "bmi", "physact", "time_bld_dx_months",
                                  "htnhx", "dmhx", "osteohx", "chdhx", 'pmh',
                                  "dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti",
                                  "bxc", "osteocalcin", "inpth", "vitd25","calcium","sclerostin", "thiaz", "bisphos"),
                        factorVars = c( "white", "smkst",'pmh',"htnhx", "dmhx", "osteohx", "chdhx", "thiaz", "bisphos"),
                        data = d[which(d$cohort=='nhs1'),], strata = "case")
table1<-print(table1,nonnormal = c("time_bld_dx_months"), cramVars = "hepato", 
              quote = FALSE, noSpaces = TRUE, showAllLevels = FALSE)

table2<- CreateTableOne(vars = c( "white", "smkst", "age", "bmi", "physact", "time_bld_dx_months",
                                  "htnhx", "dmhx", "osteohx", "chdhx", 
                                  "dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti",
                                  "bxc", "osteocalcin", "inpth", "vitd25","calcium","sclerostin", "thiaz", "bisphos"),
                        factorVars = c( "white", "smkst", "htnhx", "dmhx", "osteohx", "chdhx", "thiaz", "bisphos"),
                        data = d[which(d$cohort=='hpfs'),], strata = "case")
table2<-print(table2,nonnormal = c("time_bld_dx_months"), cramVars = "hepato", 
              quote = FALSE, noSpaces = TRUE, showAllLevels = FALSE)

t2=rbind(table1, table2)
write.csv(t2, "table_s1.csv")


#-------------------------------------------------------------------------------
# (3). Table 2
#-------------------------------------------------------------------------------
# names(e_dat_network)
table<- CreateTableOne(vars = c( "White", "smoke", "ageyr", "bmi","act", "time_bld_dx_months",
                                  "alco","pmh","osteohx","case","caco",'endpoint'),
                        factorVars = c( "White", "smoke","pmh","osteohx","case","caco",'endpoint'),
                        data = e_dat_network[which(e_dat_network$cohort!='nhs2'),], strata = "cohort")
table<-print(table,nonnormal = c("time_bld_dx_months"), cramVars = "hepato", 
              quote = FALSE, noSpaces = TRUE, showAllLevels = FALSE)
write.csv(table, "table_2.csv")


#-------------------------------------------------------------------------------
# (4). Table 3
#-------------------------------------------------------------------------------
# "X52463" "X52629" "X53197" "X53199" "X57338" "X64896"
results_mv_main_met_six= results_mv_main %>% filter(met %in% c("X52463", "X52629", "X53197", "X53199", "X57338", "X64896")) %>% 
  mutate(ORCI=paste0(round(OR_1,2),"(",round(lci_1,2),"-",round(uci_1,2),")"),Pvalue=round(p_1,3),Pfdr=round(p_adj_fdr,3)) %>% 
  select(population, SUB_PATHWAY, METABOLITE, HMDB, ORCI, Pvalue, Pfdr)

write.csv(results_mv_main_met_six, "table_3.csv")





