
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
# (1). Figure 1
#-------------------------------------------------------------------------------
# 1a-c
results_mv_map<-results_mv_main %>% mutate(group=case_when(OR_1>1 & p_1<0.05  ~ 'positive',
                                                           OR_1<1 & p_1<0.05  ~ 'negative',
                                                           T ~ 'insignificant')) %>% 
  select(met,METABOLITE,OR_1,p_1,p_adj_fdr,population,group) %>% 
  mutate(plot_name=case_when(METABOLITE=='1-(1-enyl-palmitoyl)-2-docosahexaenoyl-GPC (P-16:0/22:6)*' ~ 
                               '1-enyl-palmitoyl-2-docosahexaenoyl-GPC*',
                             METABOLITE=='1-(1-enyl-palmitoyl)-2-docosahexaenoyl-GPE (P-16:0/22:6)*' ~ 
                               '1-enyl-palmitoyl-2-docosahexaenoyl-GPE*',
                             METABOLITE=='1-(1-enyl-stearoyl)-2-docosahexaenoyl-GPE (P-18:0/22:6)*' ~ 
                               '1-enyl-stearoyl-2-docosahexaenoyl-GPE*',
                             METABOLITE=='1-(1-enyl-stearoyl)-2-oleoyl-GPC (P-18:0/18:1)' ~ 
                               '1-enyl-stearoyl-2-oleoyl-GPC',
                             T ~ METABOLITE),
         plot_name=gsub("\\(.*\\)", "", plot_name),
         plot_name=gsub("[:\\s]", "", plot_name, perl = TRUE),
         plot_name=gsub("/.*", "", plot_name)) 

results_mv_map_1<-results_mv_map %>% filter(population=='whole')
results_mv_map_2<-results_mv_map %>% filter(population=='hp')
results_mv_map_3<-results_mv_map %>% filter(population=='n1')
Pvalue=0.05
head(results_mv_map[which(results_mv_map$p_adj_fdr<0.05),])
met_show=c("X52629") # "X100009004"

#1
volcano = ggplot(results_mv_map_1, aes(x=OR_1, y=-log10(p_adj_fdr),color = group)) +
  geom_point(alpha = 0.75, size = 2) +
  labs(x = 'Adjusted Odds Ratio', 
       y = bquote(~-Log[10]~italic("FDR Adjust P-value")), 
       title = "HPFS+NHS") + # Whole population (N=7499)
  scale_colour_manual(name = "", values = alpha(c("#d8d8d8","#2DB2EB","#EB4232"), 0.7)) + 
  geom_hline(yintercept = c(-log10(Pvalue)),size = 0.7,color = "black",lty = "dashed") + 
  theme_classic()+theme(legend.text=element_text(size=14,face = "bold"),
                        axis.title =element_text(size=12),
                        title =element_text(size=12,face = 'bold'),
                        axis.text = element_text(size = 10),
                        legend.title = element_blank(),legend.position = 'top',
                        plot.margin=margin(0,1,0,1,"cm"))
# select
Up <- filter(results_mv_map_1, met %in% met_show) 
nudge_y_up = Up$p_adj_fdr

p1 <- volcano + 
  geom_point(data = Up,aes(x = OR_1, y = -log10(p_adj_fdr)),
             color = 'grey100', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Up,aes(x = OR_1, y = -log10(p_adj_fdr), label = plot_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 2, 
                  nudge_y = nudge_y_up, 
                  direction = "y", 
                  hjust = 1, 
                  force = 5,
                  force_pull = 5,
                  size = 4,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf)

#2
volcano = ggplot(results_mv_map_2, aes(x=OR_1, y=-log10(p_adj_fdr),color = group)) +
  geom_point(alpha = 0.75, size = 2) +
  labs(x = '', 
       y = '', 
       title = "HPFS") + 
  scale_colour_manual(name = "", values = alpha(c("#d8d8d8","#2DB2EB","#EB4232"), 0.7)) + 
  geom_hline(yintercept = c(-log10(Pvalue)),size = 0.7,color = "black",lty = "dashed") + 
  theme_classic()+theme(legend.text=element_text(size=14,face = "bold"),
                        axis.title =element_text(size=12),
                        title =element_text(size=12,face = 'bold'),
                        axis.text = element_text(size = 10),
                        legend.title = element_blank(),legend.position = 'top',
                        plot.margin=margin(0,1,0,1,"cm"))
# select top 3
Up <- filter(results_mv_map_2, met %in% met_show) 
nudge_y_up = Up$p_adj_fdr+0.05

p2 <- volcano + 
  geom_point(data = Up,aes(x = OR_1, y = -log10(p_adj_fdr)),
             color = 'grey100', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Up,aes(x = OR_1, y = -log10(p_adj_fdr), label = plot_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 2, 
                  nudge_y = nudge_y_up, 
                  direction = "y", 
                  hjust = 1, 
                  force = 5,
                  force_pull = 5,
                  size = 4,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf)

#3
volcano = ggplot(results_mv_map_3, aes(x=OR_1, y=-log10(p_adj_fdr),color = group)) +
  geom_point(alpha = 0.75, size = 2) +
  labs(x = '', 
       y = '', 
       title = "NHS") + 
  scale_colour_manual(name = "", values = alpha(c("#d8d8d8","#2DB2EB","#EB4232"), 0.7)) + 
  geom_hline(yintercept = c(-log10(Pvalue)),size = 0.7,color = "black",lty = "dashed") + 
  theme_classic()+theme(legend.text=element_text(size=14,face = "bold"),
                        axis.title =element_text(size=12),
                        title =element_text(size=12,face = 'bold'),
                        axis.text = element_text(size = 10),
                        legend.title = element_blank(),legend.position = 'top',
                        plot.margin=margin(0,1,0,1,"cm"))
# select top 3
Up <- filter(results_mv_map_3, met %in% met_show)
nudge_y_up = Up$p_adj_fdr+0.1

p3 <- volcano + 
  geom_point(data = Up,aes(x = OR_1, y = -log10(p_adj_fdr)),
             color = 'grey100', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Up,aes(x = OR_1, y = -log10(p_adj_fdr), label = plot_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 2, 
                  nudge_y = nudge_y_up, 
                  direction = "y", 
                  hjust = 3, 
                  force = 5,
                  force_pull = 5,
                  size = 4,
                  box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = Inf)


f_1=ggpubr::ggarrange(p1, p2, p3,labels = c('1a', '1b', '1c'), 
                      nrow=3, common.legend = TRUE, legend =  "top")


# 1d
results_mv_main_select<-results_mv_main %>% select(met, OR_1, p_1, population, METABOLITE, SUPER_PATHWAY, HMDB) %>% 
  mutate(HMDB = sub(",.*", "", HMDB)) %>% filter(p_1<0.05 & population %in% c('whole','hp','n1'))

results_mv_main_select$population<-factor(results_mv_main_select$population,levels = c('whole','hp','n1'),
                                          labels = c('Whole','HPFS','NHS'))

results_map<-results_mv_main_select
# names(results_map)
results_map<-results_map %>% group_by(met) %>% mutate(n=1,sum=sum(n), OR=OR_1) %>% select(-n) %>% ungroup() %>% 
  select(met,sum,METABOLITE,OR,SUPER_PATHWAY,population) %>% 
  mutate(plot_name=case_when(METABOLITE=='1-(1-enyl-palmitoyl)-2-docosahexaenoyl-GPC (P-16:0/22:6)*' ~ 
                               '1-enyl-palmitoyl-2-docosahexaenoyl-GPC*',
                             METABOLITE=='1-(1-enyl-palmitoyl)-2-docosahexaenoyl-GPE (P-16:0/22:6)*' ~ 
                               '1-enyl-palmitoyl-2-docosahexaenoyl-GPE*',
                             METABOLITE=='1-(1-enyl-stearoyl)-2-docosahexaenoyl-GPE (P-18:0/22:6)*' ~ 
                               '1-enyl-stearoyl-2-docosahexaenoyl-GPE*',
                             METABOLITE=='1-(1-enyl-stearoyl)-2-oleoyl-GPC (P-18:0/18:1)' ~ 
                               '1-enyl-stearoyl-2-oleoyl-GPC',
                             METABOLITE=='12,13-DiHOME' ~ 
                               '12-13-DiHOME',
                             METABOLITE=="2,2'-methylenebis" ~ 
                               '2-methylenebis',
                             METABOLITE=='2,4-di-tert-butylphenol' ~ 
                               '2-4-di-tert-butylphenol',
                             METABOLITE=="4'-hydroxypropiophenonesulfate" ~ 
                               '4-hydroxypropiophenonesulfate',
                             METABOLITE=="3-bromo-5-chloro-2,6-dihydroxybenzoicacid*" ~ 
                               '3-bromo-5-chloro-2-6-dihydroxybenzoicacid*',
                             METABOLITE=="3,5-dichloro-2,6-dihydroxybenzoicacid" ~ 
                               '3-5-dichloro-2-6-dihydroxybenzoicacid',
                             METABOLITE=="1,2-dipalmitoyl-GPE*" ~ 
                               '1-2-dipalmitoyl-GPE*',
                             METABOLITE=="1,2-dipalmitoyl-GPE*" ~ 
                               '1-2-dipalmitoyl-GPE*',
                             T ~ METABOLITE),
         plot_name = gsub("[',]", "-", plot_name),
         plot_name=gsub("\\(.*\\)", "", plot_name),
         plot_name=gsub("[:\\s]", "", plot_name, perl = TRUE),
         plot_name=gsub("/.*", "", plot_name),
         SUPER_PATHWAY=case_when(is.na(SUPER_PATHWAY) ~ 'N/A',
                               T ~ SUPER_PATHWAY)) 
#dput(table(results_map$group))
results_map$group<-factor(results_map$population)
# summary(results_map$OR)
results_map <-results_map %>% arrange(SUPER_PATHWAY, desc(sum)) 
names(results_map)[5]<-'Class'
names(results_map)[4]<-'Risk of hip fracture'
#table(results_map$group,is.na(results_map$`Risk of hip fracture`))
#display.brewer.all(type = "div")
hm.palette <- rev(brewer.pal(11, "RdYlBu"))

species <- results_map %>% filter(!duplicated(plot_name)) %>% 
  pull(plot_name)
newick_string <- paste("(", paste(species, collapse = ","), ");", sep = "")
tree <- read.tree(text = newick_string)
# summary(results_map$`Risk of hip fracture`)
f_2=ggtree(tree, branch.length = "none", layout = "circular",linetype = 0, size = 0.5) +
  layout_fan(angle = 60)+ 
  
  geom_fruit(data=results_map,geom=geom_tile,
             mapping=aes(y=plot_name,fill=`Class`),
             width=0.06,offset=0.001)+
  scale_fill_frontiers()+
  
  geom_tiplab(offset=0.5,size=4,color = "black")+
  new_scale_fill()+
  
  geom_fruit(data=results_map[which(!is.na(results_map$`Risk of hip fracture`)),],
             geom=geom_tile,
             mapping=aes(y=plot_name,x=group,fill=`Risk of hip fracture`),
             pwidth=0.3, offset=0.02,
             axis.params=list(axis="x", text.angle=-90,text.size=4,hjust=0))+
  #cale_fill_gradient2(colours = hm.palette, limits = c(0.5, 2.5))+
  scale_fill_gradient2(midpoint = 1,low = '#3C8DAD',mid = "white",high = '#FF6767',
                       limits=c(0.5,2.5)) +
  theme(plot.margin=margin(3,1,2,5,"cm"),
        legend.background = element_blank(),
        legend.margin=margin(0,0,0,3,"cm"),
        legend.spacing = unit(2,'cm'),
        legend.title = element_text(size=18))


ggpubr::ggarrange(f_1, f_2,labels = c('', '1d'), ncol=2, widths = c(1, 2.6))
ggsave(filename = "figure 1.pdf", 
       width = 16,
       height = 10,
       dpi = 600)



#-------------------------------------------------------------------------------
# (2). Figure 2
#-------------------------------------------------------------------------------
results_mv_main=results_mv_main %>% mutate(beta=log(OR), beta_1=log(OR_1))

pathway=c('Chemical', 'Acetylated Peptides', 'Benzoate Metabolism',
          'Diacylglycerol', 'Fatty Acid Metabolism (Acyl Choline)', 'Fatty Acid, Dicarboxylate',
          'Food Component/Plant', 'Pregnenolone Steroids', 
          'Phosphatidylcholine (PC)','Phosphatidylethanolamine (PE)')

# head(results_mv_main[which(results_mv_main$p_adj_fdr<0.05),])
results_mv_whole<-results_mv_main %>% filter(population=='whole')
gsea_union=results_mv_whole %>% select(METABOLITE,SUB_PATHWAY,beta,OR,p,beta_1,OR_1,p_1)

met.classes = vector("list",length = length(unique(gsea_union$SUB_PATHWAY)))
names(met.classes) = unique(gsea_union$SUB_PATHWAY)
met.classes[["NA"]] = NULL

for(i in 1:length(met.classes)){
  met.classes[[i]] = gsea_union$METABOLITE[which(gsea_union$SUB_PATHWAY == names(met.classes)[i])]
}

runFGSEA = function(pathways = met.classes, stats = ovca.mets, minSize=3, maxSize=500, file.result){
  fgsea.ovca <- fgseaMultilevel(pathways = pathways, stats = stats, minSize=minSize, maxSize=maxSize)
  fgsea.ovca = fgsea.ovca[order(fgsea.ovca$padj),]
  print(fgsea.ovca)
  to.print = fgsea.ovca[,-8]
  write.table(x = to.print, file = file.result ,col.names = T, row.names = F, sep = "\t", quote = F)
  return(to.print)
}

set.seed(1234)
m_unadj_mets = gsea_union$beta
names(m_unadj_mets) = gsea_union$METABOLITE

m_unadj_mets = m_unadj_mets[which(!is.na(m_unadj_mets))]
m_unadj_fgsea = runFGSEA(pathways = met.classes, stats = m_unadj_mets, 
                         file.result = "/udd/n2zwu/Meta_hip/gsea_m1_cont.csv")

m_adj_mets = gsea_union$beta_1
names(m_adj_mets) = gsea_union$METABOLITE
m_adj_mets = m_adj_mets[which(!is.na(m_adj_mets))]
m_adj_fgsea = runFGSEA(pathways = met.classes, stats = m_adj_mets, 
                       file.result = "/udd/n2zwu/Meta_hip/gsea_m2_cont.csv")

fgsea.res = merge(x = m_unadj_fgsea[,c(1,2,3,6)], by.x = "pathway", y = m_adj_fgsea[,c(1,2,3,6)], 
                  by.y = "pathway",all = TRUE, suffixes = c(".Unadjusted",".Adjusted"))

# plot
data.padj = fgsea.res[,c(1,3,6)]
data.nes = fgsea.res[,c(1,4,7)]
colnames(data.nes) = c("Pathway","Unadjusted","Adjusted")
colnames(data.padj) = c("Pathway","Unadjusted","Adjusted")

nes = tidyr::pivot_longer(data.nes, 
                          cols = c("Unadjusted","Adjusted"), 
                          names_to = "model", 
                          values_to = "value")
padj = tidyr::pivot_longer(data.padj, 
                           cols = c("Unadjusted","Adjusted"), 
                           names_to = "model", 
                           values_to = "pvalue")

nes_sig<-merge(nes, padj, by = c('Pathway','model'),all = FALSE, sort = TRUE)
table(nes_sig[which(nes_sig$model=='Adjusted' & nes_sig$pvalue<0.2),]$Pathway)

nes_sig_both<-nes_sig %>% mutate(stars=case_when(pvalue<=0.001  ~ '***',
                                                 pvalue<=0.05  ~ '**',
                                                 pvalue<=0.2  ~ '*', 
                                                 T ~ ''),
                                 pstar=paste0(round(pvalue, 3), " ", stars)) %>% 
  filter(model=='Adjusted' & Pathway %in% pathway) %>% 
  select(Pathway, value, pvalue, pstar,stars) %>% mutate(model='Whole')


# hp
results_mv_hp<-results_mv_main %>% filter(population=='hp')
gsea_union=results_mv_hp %>% select(METABOLITE,SUB_PATHWAY,beta,OR,p,beta_1,OR_1,p_1)
met.classes = vector("list",length = length(unique(gsea_union$SUB_PATHWAY)))
names(met.classes) = unique(gsea_union$SUB_PATHWAY)
met.classes[["NA"]] = NULL

for(i in 1:length(met.classes)){
  met.classes[[i]] = gsea_union$METABOLITE[which(gsea_union$SUB_PATHWAY == names(met.classes)[i])]
}

set.seed(1234)
m_unadj_mets = gsea_union$beta
names(m_unadj_mets) = gsea_union$METABOLITE

m_unadj_mets = m_unadj_mets[which(!is.na(m_unadj_mets))]
m_unadj_fgsea = runFGSEA(pathways = met.classes, stats = m_unadj_mets, 
                         file.result = "/udd/n2zwu/Meta_hip/gsea_m1_cont.csv")

m_adj_mets = gsea_union$beta_1
names(m_adj_mets) = gsea_union$METABOLITE
m_adj_mets = m_adj_mets[which(!is.na(m_adj_mets))]
m_adj_fgsea = runFGSEA(pathways = met.classes, stats = m_adj_mets, 
                       file.result = "/udd/n2zwu/Meta_hip/gsea_m2_cont.csv")

fgsea.res = merge(x = m_unadj_fgsea[,c(1,2,3,6)], by.x = "pathway", y = m_adj_fgsea[,c(1,2,3,6)], 
                  by.y = "pathway",all = TRUE, suffixes = c(".Unadjusted",".Adjusted"))

# plot
data.padj = fgsea.res[,c(1,3,6)]
data.nes = fgsea.res[,c(1,4,7)]
colnames(data.nes) = c("Pathway","Unadjusted","Adjusted")
colnames(data.padj) = c("Pathway","Unadjusted","Adjusted")

nes = tidyr::pivot_longer(data.nes, 
                          cols = c("Unadjusted","Adjusted"), 
                          names_to = "model", 
                          values_to = "value")
padj = tidyr::pivot_longer(data.padj, 
                           cols = c("Unadjusted","Adjusted"), 
                           names_to = "model", 
                           values_to = "pvalue")

nes_sig<-merge(nes, padj, by = c('Pathway','model'),all = FALSE, sort = TRUE)

table(nes_sig[which(nes_sig$model=='Adjusted' & nes_sig$pvalue<0.2),]$Pathway)
nes_sig_hp<-nes_sig %>% mutate(stars=case_when(pvalue<=0.001  ~ '***',
                                               pvalue<=0.05  ~ '**',
                                               pvalue<=0.2  ~ '*', 
                                               T ~ ''),
                               pstar=paste0(round(pvalue, 3), " ", stars)) %>% 
  filter(model=='Adjusted' & Pathway %in% pathway) %>% 
  select(Pathway, value, pvalue, pstar,stars) %>% mutate(model='HPFS')

# n1
results_mv_n1<-results_mv_main %>% filter(population=='n1')
gsea_union=results_mv_n1 %>% select(METABOLITE,SUB_PATHWAY,beta,OR,p,beta_1,OR_1,p_1)

met.classes = vector("list",length = length(unique(gsea_union$SUB_PATHWAY)))
names(met.classes) = unique(gsea_union$SUB_PATHWAY)
met.classes[["NA"]] = NULL

for(i in 1:length(met.classes)){
  met.classes[[i]] = gsea_union$METABOLITE[which(gsea_union$SUB_PATHWAY == names(met.classes)[i])]
}

set.seed(1234)
m_unadj_mets = gsea_union$beta
names(m_unadj_mets) = gsea_union$METABOLITE

m_unadj_mets = m_unadj_mets[which(!is.na(m_unadj_mets))]
m_unadj_fgsea = runFGSEA(pathways = met.classes, stats = m_unadj_mets, 
                         file.result = "/udd/n2zwu/Meta_hip/gsea_m1_cont.csv")

m_adj_mets = gsea_union$beta_1
names(m_adj_mets) = gsea_union$METABOLITE
m_adj_mets = m_adj_mets[which(!is.na(m_adj_mets))]
m_adj_fgsea = runFGSEA(pathways = met.classes, stats = m_adj_mets, 
                       file.result = "/udd/n2zwu/Meta_hip/gsea_m2_cont.csv")

fgsea.res = merge(x = m_unadj_fgsea[,c(1,2,3,6)], by.x = "pathway", y = m_adj_fgsea[,c(1,2,3,6)], 
                  by.y = "pathway",all = TRUE, suffixes = c(".Unadjusted",".Adjusted"))

# plot
data.padj = fgsea.res[,c(1,3,6)]
data.nes = fgsea.res[,c(1,4,7)]
colnames(data.nes) = c("Pathway","Unadjusted","Adjusted")
colnames(data.padj) = c("Pathway","Unadjusted","Adjusted")
nes = tidyr::pivot_longer(data.nes, 
                          cols = c("Unadjusted","Adjusted"), 
                          names_to = "model", 
                          values_to = "value")
padj = tidyr::pivot_longer(data.padj, 
                           cols = c("Unadjusted","Adjusted"), 
                           names_to = "model", 
                           values_to = "pvalue")
nes_sig<-merge(nes, padj, by = c('Pathway','model'),all = FALSE, sort = TRUE)

table(nes_sig[which(nes_sig$model=='Adjusted' & nes_sig$pvalue<0.2),]$Pathway)
nes_sig_n1<-nes_sig %>% mutate(stars=case_when(pvalue<=0.001  ~ '***',
                                               pvalue<=0.05  ~ '**',
                                               pvalue<=0.2  ~ '*', 
                                               T ~ ''),
                               pstar=paste0(round(pvalue, 3), " ", stars)) %>% 
  filter(model=='Adjusted' & Pathway %in% pathway) %>% 
  select(Pathway, value, pvalue, pstar,stars) %>% mutate(model='NHS')

nes_sig<-rbind(nes_sig_both,nes_sig_hp,nes_sig_n1)
#head(nes_sig[which(nes_sig$model=='Whole' & nes_sig$pvalue<0.05),])
#range(nes_sig$value)

#hm.palette <- rev(brewer.pal(11, "RdBu"))
hm.palette <- rev(brewer.pal(11, "RdYlBu"))
nes_sig$model<-factor(nes_sig$model, levels = c('Whole', 'HPFS','NHS'), labels = c('HPFS+NHS', 'HPFS','NHS'))
nes_sig_out_1=nes_sig

ggplot(data = nes_sig, aes(x = model, y = fct_reorder(Pathway, desc(pvalue)), fill = value)) + 
  geom_tile(color = "white", linewidth = 0.1) + 
  scale_fill_gradientn(colours = hm.palette, limits = c(-3,3))+
  geom_text(aes(x = model, y = Pathway, label=stars), color="white", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Metabolite classes\nassociated with hip fracture")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
                                   size = 11, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 11, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  guides(fill=guide_colorbar(title = "MSEA Enrichment Score",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = "* FDR<0.2; ** FDR<0.05")

ggsave(filename = "figure 2.pdf", 
       width = 7,
       height = 5,
       dpi = 600)

# head(nes_sig_both, 10)
# head(nes_sig_hp, 10)
# head(nes_sig_n1, 10)


#-------------------------------------------------------------------------------
# (3). Figure 2s -- stratify
#-------------------------------------------------------------------------------
results_mv_main=results_mv_main %>% mutate(beta=log(OR), beta_1=log(OR_1))

pathway=c('Chemical', 'Acetylated Peptides', 'Benzoate Metabolism',
          'Diacylglycerol', 'Fatty Acid Metabolism (Acyl Choline)', 'Fatty Acid, Dicarboxylate',
          'Food Component/Plant', 'Pregnenolone Steroids', 
          'Phosphatidylcholine (PC)','Phosphatidylethanolamine (PE)')

# head(results_mv_main[which(results_mv_main$p_adj_fdr<0.05),])
results_mv_whole<-results_mv_main %>% filter(population=='whole')
gsea_union=results_mv_whole %>% select(METABOLITE,SUB_PATHWAY,beta,OR,p,beta_1,OR_1,p_1)

met.classes = vector("list",length = length(unique(gsea_union$SUB_PATHWAY)))
names(met.classes) = unique(gsea_union$SUB_PATHWAY)
met.classes[["NA"]] = NULL

for(i in 1:length(met.classes)){
  met.classes[[i]] = gsea_union$METABOLITE[which(gsea_union$SUB_PATHWAY == names(met.classes)[i])]
}

runFGSEA = function(pathways = met.classes, stats = ovca.mets, minSize=3, maxSize=500, file.result){
  fgsea.ovca <- fgseaMultilevel(pathways = pathways, stats = stats, minSize=minSize, maxSize=maxSize)
  fgsea.ovca = fgsea.ovca[order(fgsea.ovca$padj),]
  print(fgsea.ovca)
  to.print = fgsea.ovca[,-8]
  write.table(x = to.print, file = file.result ,col.names = T, row.names = F, sep = "\t", quote = F)
  return(to.print)
}

set.seed(1234)
m_unadj_mets = gsea_union$beta
names(m_unadj_mets) = gsea_union$METABOLITE

m_unadj_mets = m_unadj_mets[which(!is.na(m_unadj_mets))]
m_unadj_fgsea = runFGSEA(pathways = met.classes, stats = m_unadj_mets, 
                         file.result = "/udd/n2zwu/Meta_hip/gsea_m1_cont.csv")

m_adj_mets = gsea_union$beta_1
names(m_adj_mets) = gsea_union$METABOLITE
m_adj_mets = m_adj_mets[which(!is.na(m_adj_mets))]
m_adj_fgsea = runFGSEA(pathways = met.classes, stats = m_adj_mets, 
                       file.result = "/udd/n2zwu/Meta_hip/gsea_m2_cont.csv")

fgsea.res = merge(x = m_unadj_fgsea[,c(1,2,3,6)], by.x = "pathway", y = m_adj_fgsea[,c(1,2,3,6)], 
                  by.y = "pathway",all = TRUE, suffixes = c(".Unadjusted",".Adjusted"))

# plot
data.padj = fgsea.res[,c(1,3,6)]
data.nes = fgsea.res[,c(1,4,7)]
colnames(data.nes) = c("Pathway","Unadjusted","Adjusted")
colnames(data.padj) = c("Pathway","Unadjusted","Adjusted")

nes = tidyr::pivot_longer(data.nes, 
                          cols = c("Unadjusted","Adjusted"), 
                          names_to = "model", 
                          values_to = "value")
padj = tidyr::pivot_longer(data.padj, 
                           cols = c("Unadjusted","Adjusted"), 
                           names_to = "model", 
                           values_to = "pvalue")

nes_sig<-merge(nes, padj, by = c('Pathway','model'),all = FALSE, sort = TRUE)
table(nes_sig[which(nes_sig$model=='Adjusted' & nes_sig$pvalue<0.2),]$Pathway)

nes_sig_both<-nes_sig %>% mutate(stars=case_when(pvalue<=0.001  ~ '***',
                                                 pvalue<=0.05  ~ '**',
                                                 pvalue<=0.2  ~ '*', 
                                                 T ~ ''),
                                 pstar=paste0(round(pvalue, 3), " ", stars)) %>% 
  filter(model=='Adjusted' & Pathway %in% pathway) %>% 
  select(Pathway, value, pvalue, pstar,stars) %>% mutate(model='Whole')


# bld_dx_hip <= 120
results_mv_hp<-results_mv_main %>% filter(population=='below')
gsea_union=results_mv_hp %>% select(METABOLITE,SUB_PATHWAY,beta,OR,p,beta_1,OR_1,p_1)
met.classes = vector("list",length = length(unique(gsea_union$SUB_PATHWAY)))
names(met.classes) = unique(gsea_union$SUB_PATHWAY)
met.classes[["NA"]] = NULL

for(i in 1:length(met.classes)){
  met.classes[[i]] = gsea_union$METABOLITE[which(gsea_union$SUB_PATHWAY == names(met.classes)[i])]
}

set.seed(1234)
m_unadj_mets = gsea_union$beta
names(m_unadj_mets) = gsea_union$METABOLITE

m_unadj_mets = m_unadj_mets[which(!is.na(m_unadj_mets))]
m_unadj_fgsea = runFGSEA(pathways = met.classes, stats = m_unadj_mets, 
                         file.result = "/udd/n2zwu/Meta_hip/gsea_m1_cont.csv")

m_adj_mets = gsea_union$beta_1
names(m_adj_mets) = gsea_union$METABOLITE
m_adj_mets = m_adj_mets[which(!is.na(m_adj_mets))]
m_adj_fgsea = runFGSEA(pathways = met.classes, stats = m_adj_mets, 
                       file.result = "/udd/n2zwu/Meta_hip/gsea_m2_cont.csv")

fgsea.res = merge(x = m_unadj_fgsea[,c(1,2,3,6)], by.x = "pathway", y = m_adj_fgsea[,c(1,2,3,6)], 
                  by.y = "pathway",all = TRUE, suffixes = c(".Unadjusted",".Adjusted"))

# plot
data.padj = fgsea.res[,c(1,3,6)]
data.nes = fgsea.res[,c(1,4,7)]
colnames(data.nes) = c("Pathway","Unadjusted","Adjusted")
colnames(data.padj) = c("Pathway","Unadjusted","Adjusted")

nes = tidyr::pivot_longer(data.nes, 
                          cols = c("Unadjusted","Adjusted"), 
                          names_to = "model", 
                          values_to = "value")
padj = tidyr::pivot_longer(data.padj, 
                           cols = c("Unadjusted","Adjusted"), 
                           names_to = "model", 
                           values_to = "pvalue")

nes_sig<-merge(nes, padj, by = c('Pathway','model'),all = FALSE, sort = TRUE)

table(nes_sig[which(nes_sig$model=='Adjusted' & nes_sig$pvalue<0.2),]$Pathway)
nes_sig_below<-nes_sig %>% mutate(stars=case_when(pvalue<=0.001  ~ '***',
                                               pvalue<=0.05  ~ '**',
                                               pvalue<=0.2  ~ '*', 
                                               T ~ ''),
                               pstar=paste0(round(pvalue, 3), " ", stars)) %>% 
  filter(model=='Adjusted' & Pathway %in% pathway) %>% 
  select(Pathway, value, pvalue, pstar,stars) %>% mutate(model='below')


# bld_dx_hip > 120
results_mv_n1<-results_mv_main %>% filter(population=='above')
gsea_union=results_mv_n1 %>% select(METABOLITE,SUB_PATHWAY,beta,OR,p,beta_1,OR_1,p_1)

met.classes = vector("list",length = length(unique(gsea_union$SUB_PATHWAY)))
names(met.classes) = unique(gsea_union$SUB_PATHWAY)
met.classes[["NA"]] = NULL

for(i in 1:length(met.classes)){
  met.classes[[i]] = gsea_union$METABOLITE[which(gsea_union$SUB_PATHWAY == names(met.classes)[i])]
}

set.seed(1234)
m_unadj_mets = gsea_union$beta
names(m_unadj_mets) = gsea_union$METABOLITE

m_unadj_mets = m_unadj_mets[which(!is.na(m_unadj_mets))]
m_unadj_fgsea = runFGSEA(pathways = met.classes, stats = m_unadj_mets, 
                         file.result = "/udd/n2zwu/Meta_hip/gsea_m1_cont.csv")

m_adj_mets = gsea_union$beta_1
names(m_adj_mets) = gsea_union$METABOLITE
m_adj_mets = m_adj_mets[which(!is.na(m_adj_mets))]
m_adj_fgsea = runFGSEA(pathways = met.classes, stats = m_adj_mets, 
                       file.result = "/udd/n2zwu/Meta_hip/gsea_m2_cont.csv")

fgsea.res = merge(x = m_unadj_fgsea[,c(1,2,3,6)], by.x = "pathway", y = m_adj_fgsea[,c(1,2,3,6)], 
                  by.y = "pathway",all = TRUE, suffixes = c(".Unadjusted",".Adjusted"))

# plot
data.padj = fgsea.res[,c(1,3,6)]
data.nes = fgsea.res[,c(1,4,7)]
colnames(data.nes) = c("Pathway","Unadjusted","Adjusted")
colnames(data.padj) = c("Pathway","Unadjusted","Adjusted")
nes = tidyr::pivot_longer(data.nes, 
                          cols = c("Unadjusted","Adjusted"), 
                          names_to = "model", 
                          values_to = "value")
padj = tidyr::pivot_longer(data.padj, 
                           cols = c("Unadjusted","Adjusted"), 
                           names_to = "model", 
                           values_to = "pvalue")
nes_sig<-merge(nes, padj, by = c('Pathway','model'),all = FALSE, sort = TRUE)

table(nes_sig[which(nes_sig$model=='Adjusted' & nes_sig$pvalue<0.2),]$Pathway)
nes_sig_above<-nes_sig %>% mutate(stars=case_when(pvalue<=0.001  ~ '***',
                                               pvalue<=0.05  ~ '**',
                                               pvalue<=0.2  ~ '*', 
                                               T ~ ''),
                               pstar=paste0(round(pvalue, 3), " ", stars)) %>% 
  filter(model=='Adjusted' & Pathway %in% pathway) %>% 
  select(Pathway, value, pvalue, pstar,stars) %>% mutate(model='above')

nes_sig<-rbind(nes_sig_both,nes_sig_below,nes_sig_above)
#head(nes_sig[which(nes_sig$model=='Whole' & nes_sig$pvalue<0.05),])
#range(nes_sig$value)

#hm.palette <- rev(brewer.pal(11, "RdBu"))
hm.palette <- rev(brewer.pal(11, "RdYlBu"))
nes_sig$model<-factor(nes_sig$model, levels = c('Whole', 'below', 'above'), 
                      labels = c('Whole Hip fracture', 'Occurring within 10 years', 'Occurring beyond 10 years'))

nes_sig_out_2=nes_sig

ggplot(data = nes_sig, aes(x = model, y = fct_reorder(Pathway, desc(pvalue)), fill = value)) + 
  geom_tile(color = "white", linewidth = 0.1) + 
  scale_fill_gradientn(colours = hm.palette, limits = c(-3,3))+
  geom_text(aes(x = model, y = Pathway, label=stars), color="white", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Metabolite classes\nassociated with hip fracture\noccurring within 10 years or beyond 10 years")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, 
                                   size = 11, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 11, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  guides(fill=guide_colorbar(title = "MSEA Enrichment Score",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = "* FDR<0.2; ** FDR<0.05")

ggsave(filename = "figure 2s.pdf", 
       width = 7,
       height = 5,
       dpi = 600)


#-------------------------------------------------------------------------------
# (4). Figure 3
#-------------------------------------------------------------------------------
pathway=c('Acetylated Peptides','Chemical','Fatty Acid Metabolism (Acyl Choline)',
          'Fatty Acid, Dicarboxylate','Food Component/Plant','Diacylglycerol','Benzoate Metabolism',
          'Phosphatidylcholine (PC)','Phosphatidylethanolamine (PE)','Pregnenolone Steroids')

met_sig<-results_mv_main %>% filter(SUB_PATHWAY %in% pathway & p_1<0.05 & population=='whole') %>% 
  filter(!duplicated(met)) %>% pull(met) # 31 metabolites

results_mv_sub<-results_mv_main %>% filter(met %in% met_sig & population %in% c('whole','hp','n1')) %>% 
  mutate(stars=case_when(p_adj_fdr<=0.05   ~ '**',
                         #p_1<=0.001  ~ '***',
                         #p_1<=0.01  ~ '**',
                         p_1<=0.05  ~ '*', 
                         T ~ ''),
         pstar=paste0(round(p_1, 2), " ", stars))

results_mv_sub$population<-factor(results_mv_sub$population,levels = c('whole','hp','n1'),labels = c('HPFS+NHS','HPFS', 'NHS'))
# table(results_mv_sub$SUB_PATHWAY)

results_mv_pos<-results_mv_sub %>% filter(SUB_PATHWAY %in% c('Acetylated Peptides','Chemical','Fatty Acid, Dicarboxylate', 
                                                             'Food Component/Plant','Benzoate Metabolism'))
results_mv_neg<-results_mv_sub %>% filter(SUB_PATHWAY %in% c('Phosphatidylcholine (PC)','Phosphatidylethanolamine (PE)',
                                                             'Fatty Acid Metabolism (Acyl Choline)'))

results_mv_pos$SUB_PATHWAY<-factor(results_mv_pos$SUB_PATHWAY,levels = c('Chemical',
                                                                         'Fatty Acid, Dicarboxylate',
                                                                         'Acetylated Peptides',
                                                                         'Food Component/Plant',
                                                                         'Benzoate Metabolism'))
results_mv_neg$SUB_PATHWAY<-factor(results_mv_neg$SUB_PATHWAY,levels = c('Phosphatidylcholine (PC)',
                                                                         'Phosphatidylethanolamine (PE)',
                                                                         'Fatty Acid Metabolism (Acyl Choline)'))

f_1<-ggplot(data = results_mv_pos,
            aes(x = population, y =fct_reorder(METABOLITE, desc(SUB_PATHWAY)), fill = OR_1)) + #data and lables
  geom_tile(color = "white", size = 0.1) + # what to plot
  #scale_fill_gradientn(colours = hm.palette, limits = c(0.5,2.5))+
  scale_fill_gradient2(midpoint = 1,low = '#3C8DAD',mid = "white",high = '#FF6767',
                       limits=c(0.5,2.5)) +
  geom_text(aes(x = population, y = METABOLITE, label=stars), color="black", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Positive metabolite classes")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
                                   size = 12, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 1, 2, "cm"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 10, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  geom_tile(data = results_mv_pos %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = 0, y = METABOLITE, col = factor(SUB_PATHWAY)),fill='white', 
            width = 0.1, height = 1, size=6, inherit.aes = FALSE) +
  scale_colour_lancet () +
  geom_text(data = results_mv_pos %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = -0.2, y = METABOLITE, label = SUB_PATHWAY), color = "black", size = 3, 
            hjust = 1, inherit.aes = FALSE)+
  guides(fill=guide_colorbar(title = "adjusted OR"),
         color = guide_legend(title = "Metabolite classes")) +
  theme(legend.box = "vertical", legend.position = "right")

f_2<-ggplot(data = results_mv_neg,
            aes(x = population, y =fct_reorder(METABOLITE, desc(SUB_PATHWAY)), fill = OR_1)) + #data and lables
  geom_tile(color = "white", size = 0.1) + # what to plot
  #scale_fill_gradientn(colours = hm.palette, limits = c(0.5,1.5))+
  scale_fill_gradient2(midpoint = 1,low = '#3C8DAD',mid = "white",high = '#FF6767',
                       limits=c(0.5,1.5)) +
  geom_text(aes(x = population, y = METABOLITE, label=stars), color="black", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Negative metabolite classes")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
                                   size = 12, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 1, 2, "cm"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 10, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  geom_tile(data = results_mv_neg %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = 0, y = METABOLITE, col = factor(SUB_PATHWAY)),fill='white', 
            width = 0.1, height = 1, size=6, inherit.aes = FALSE) +
  scale_colour_frontiers () +
  geom_text(data = results_mv_neg %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = -0.2, y = METABOLITE, label = SUB_PATHWAY), color = "black", size = 3, 
            hjust = 1, inherit.aes = FALSE)+
  labs(caption = "* P<0.05; ** FDR<0.05") +
  guides(fill=guide_colorbar(title = "adjusted OR"),
         color = guide_legend(title = "Metabolite classes")) +
  theme(legend.box = "vertical", legend.position = "right")

ggpubr::ggarrange(f_1, f_2, labels = c('3a', '3b'), ncol=1, nrow=2, 
                  heights = c(1,1.4), common.legend = F, align = 'v')

ggsave(filename = "figure 3.pdf", 
       width = 12,
       height = 14,
       dpi = 600)

#-------------------------------------------------------------------------------
# (5). Figure 3s --  stratify
#-------------------------------------------------------------------------------
pathway=c('Acetylated Peptides','Chemical','Fatty Acid Metabolism (Acyl Choline)',
          'Fatty Acid, Dicarboxylate','Food Component/Plant','Diacylglycerol','Benzoate Metabolism',
          'Phosphatidylcholine (PC)','Phosphatidylethanolamine (PE)','Pregnenolone Steroids')

met_sig<-results_mv_main %>% filter(SUB_PATHWAY %in% pathway & p_1<0.05 & population=='whole') %>% 
  filter(!duplicated(met)) %>% pull(met) # 31 metabolites

results_mv_sub<-results_mv_main %>% filter(met %in% met_sig & population %in% c('whole','below','above')) %>% 
  mutate(stars=case_when(p_adj_fdr<=0.05   ~ '**',
                         #p_1<=0.001  ~ '***',
                         #p_1<=0.01  ~ '**',
                         p_1<=0.05  ~ '*', 
                         T ~ ''),
         pstar=paste0(round(p_1, 2), " ", stars))

results_mv_sub$population<-factor(results_mv_sub$population,levels = c('whole','below','above'),
                                  labels = c('Whole Hip fracture', 'Occurring within 10 years', 'Occurring beyond 10 years'))
# table(results_mv_sub$SUB_PATHWAY)

results_mv_pos<-results_mv_sub %>% filter(SUB_PATHWAY %in% c('Acetylated Peptides','Chemical','Fatty Acid, Dicarboxylate', 
                                                             'Food Component/Plant','Benzoate Metabolism'))
results_mv_neg<-results_mv_sub %>% filter(SUB_PATHWAY %in% c('Phosphatidylcholine (PC)','Phosphatidylethanolamine (PE)',
                                                             'Fatty Acid Metabolism (Acyl Choline)'))

results_mv_pos$SUB_PATHWAY<-factor(results_mv_pos$SUB_PATHWAY,levels = c('Chemical',
                                                                         'Fatty Acid, Dicarboxylate',
                                                                         'Acetylated Peptides',
                                                                         'Food Component/Plant',
                                                                         'Benzoate Metabolism'))
results_mv_neg$SUB_PATHWAY<-factor(results_mv_neg$SUB_PATHWAY,levels = c('Phosphatidylcholine (PC)',
                                                                         'Phosphatidylethanolamine (PE)',
                                                                         'Fatty Acid Metabolism (Acyl Choline)'))
# range(results_mv_pos$OR_1)
f_1<-ggplot(data = results_mv_pos,
            aes(x = population, y =fct_reorder(METABOLITE, desc(SUB_PATHWAY)), fill = OR_1)) + #data and lables
  geom_tile(color = "white", size = 0.1) + # what to plot
  #scale_fill_gradientn(colours = hm.palette, limits = c(0.5,2.5))+
  scale_fill_gradient2(midpoint = 1,low = '#3C8DAD',mid = "white",high = '#FF6767',
                       limits=c(0.5,3.5)) +
  geom_text(aes(x = population, y = METABOLITE, label=stars), color="black", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Positive metabolite classes")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, 
                                   size = 12, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 1, 2, "cm"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 10, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  geom_tile(data = results_mv_pos %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = 0, y = METABOLITE, col = factor(SUB_PATHWAY)),fill='white', 
            width = 0.1, height = 1, size=6, inherit.aes = FALSE) +
  scale_colour_lancet () +
  geom_text(data = results_mv_pos %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = -0.2, y = METABOLITE, label = SUB_PATHWAY), color = "black", size = 3, 
            hjust = 1, inherit.aes = FALSE)+
  guides(fill=guide_colorbar(title = "adjusted OR"),
         color = guide_legend(title = "Metabolite classes")) +
  theme(legend.box = "vertical", legend.position = "right")

f_2<-ggplot(data = results_mv_neg,
            aes(x = population, y =fct_reorder(METABOLITE, desc(SUB_PATHWAY)), fill = OR_1)) + #data and lables
  geom_tile(color = "white", size = 0.1) + # what to plot
  #scale_fill_gradientn(colours = hm.palette, limits = c(0.5,1.5))+
  scale_fill_gradient2(midpoint = 1,low = '#3C8DAD',mid = "white",high = '#FF6767',
                       limits=c(0.5,1.5)) +
  geom_text(aes(x = population, y = METABOLITE, label=stars), color="black", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Negative metabolite classes")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, 
                                   size = 12, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 1, 2, "cm"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 10, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  geom_tile(data = results_mv_neg %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = 0, y = METABOLITE, col = factor(SUB_PATHWAY)),fill='white', 
            width = 0.1, height = 1, size=6, inherit.aes = FALSE) +
  scale_colour_frontiers () +
  geom_text(data = results_mv_neg %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = -0.2, y = METABOLITE, label = SUB_PATHWAY), color = "black", size = 3, 
            hjust = 1, inherit.aes = FALSE)+
  labs(caption = "* P<0.05; ** FDR<0.05") +
  guides(fill=guide_colorbar(title = "adjusted OR"),
         color = guide_legend(title = "Metabolite classes")) +
  theme(legend.box = "vertical", legend.position = "right")

ggpubr::ggarrange(f_1, f_2, labels = c('3a', '3b'), ncol=1, nrow=2, 
                  heights = c(1,1.4), common.legend = F, align = 'v')

ggsave(filename = "figure 3s.pdf", 
       width = 12,
       height = 14,
       dpi = 600)


#-------------------------------------------------------------------------------
# (6). Figure 4 -- remove
#-------------------------------------------------------------------------------
# results_mv_main_select<-results_mv_main %>% 
#  select(met,OR_1, p_1, population, METABOLITE,SUB_PATHWAY,HMDB) %>% 
#  mutate(HMDB = sub(",.*", "", HMDB))
#results_cox_select<-results_cox_metabolite %>% select(hr_adj, p_adj, population,met)
#names(results_cox_select)[4]<-'HMDB'

#results_mv_main_select$population<-factor(results_mv_main_select$population,levels = c('whole','hp','n1'),
#                                          labels = c('Whole','HPFS','NHS'))
#results_cox_select$population<-factor(results_cox_select$population,levels = c('whole','hpfs','nhs1'),
#                                          labels = c('Whole','HPFS','NHS'))
# table(results_mv_main_select$population)
# table(results_cox_select$population)
# table(results_mv_main_select$HMDB)

# results_map<-merge(results_mv_main_select, results_cox_select, by = c('population','HMDB'), all.x = TRUE, sort = TRUE)
# head(results_map[which(results_map$HMDB=='HMDB0010404'),], 10)

# pathway=c('Chemical','Fatty Acid Metabolism (Acyl Choline)',
#          'Phosphatidylcholine (PC)','Phosphatidylethanolamine (PE)')

#results_map<- results_map %>% filter(SUB_PATHWAY %in% pathway) %>% # mutate(hr_adj=case_when(p_adj<0.05 ~ hr_adj, T ~ NA)) %>% 
#  select(met,HMDB, population, OR_1, hr_adj, METABOLITE, SUB_PATHWAY)
# head(results_map[which(!is.na(results_map$hr_adj)),], 10)

#results_map = tidyr::pivot_longer(results_map, 
#                          cols = c("OR_1","hr_adj"), 
#                          names_to = "design", 
#                          values_to = "OR")

#results_map<-results_map %>% mutate(group=case_when(population=='Whole' & design=='OR_1' ~ 'Whole-case control',
#                                                    population=='Whole' & design=='hr_adj' ~ 'Whole-cohort',
#                                                    population=='HPFS' & design=='OR_1' ~ 'HPFS-case control',
#                                                    population=='HPFS' & design=='hr_adj' ~ 'HPFS-cohort',
#                                                    population=='NHS' & design=='OR_1' ~ 'NHS-case control',
#                                                    population=='NHS' & design=='hr_adj' ~ 'NHS-cohort'))
# table(is.na(results_map$SUB_PATHWAY))
# table(is.na(results_map$plot_name))
# table(results_map$plot_name)


#-------------------------------------------------------------------------------
# (6). Figure 4
#-------------------------------------------------------------------------------
# names(results_case_cox)
# results_cox<-results_cox %>% mutate(beta=log(hr),beta_1=log(hr_adj)) 
# table(results_case_cox$population)

results_case_cox<- results_case_cox %>% mutate(stars=case_when(pp<=0.001  ~ '***',
                                                          pp<=0.01  ~ '**',
                                                          pp<=0.05  ~ '*', 
                                                          T ~ ''),
                                          pstar=paste0(round(pp, 3), " ", stars)) 


results_case_cox$population<-factor(results_case_cox$population, levels = c('cohort-whole', 'hpfs', 'nhs1', 'whole'), 
                                    labels = c('Cohort: HPFS+NHS', 'Cohort: HPFS','Cohort: NHS','Case-control: HPFS+NHS'))

ggplot(data = results_case_cox,
           aes(x = population, y =fct_reorder(metabolite_name, desc(class_broad)), fill = rr)) + 
  geom_tile(color = "white", size = 0.1) + 
  scale_fill_gradientn(colours = hm.palette, limits = c(0.5,1.5))+
  #scale_fill_gradient2(midpoint = 1,low = '#3C8DAD',mid = "white",high = '#FF6767',limits=c(0.5,1.5)) +
  geom_text(aes(x = population, y =metabolite_name, label=stars), color="black", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Individual metabolites\nassociated with hip fracture")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 12, vjust = 0.5, 
                                   size = 12, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 1, 2, "cm"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 10, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  #new_scale_fill()+
  geom_tile(data = results_case_cox %>% distinct(metabolite_name, class_broad),
            aes(x = 0, y = metabolite_name, col = factor(class_broad)),fill='white', 
            width = 0.1, height = 1, size=6, inherit.aes = FALSE) +
  scale_colour_lancet() +
  #scale_fill_manual(values = pathway_colors) +
  geom_text(data = results_case_cox %>% distinct(metabolite_name, class_broad),
            aes(x = -0.2, y = metabolite_name, label = class_broad), color = "black", size = 3, 
            hjust = 1, inherit.aes = FALSE)+
  guides(fill=guide_colorbar(title = "Association: HR/OR",title.position = "top",nrow = 1, label.position = "bottom"),
         color = guide_legend(title = "Metabolite Class",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = "* P<0.05; ** P<0.01; *** FDR<0.05") 


# ggpubr::ggarrange(f_1, f_2, labels = c('3a', '3b'), ncol=1, nrow=2, heights = c(1,3), common.legend = F)
ggsave(filename = "figure 4.pdf", 
       width = 9,
       height = 9,
       dpi = 600)


#-------------------------------------------------------------------------------
# (7). Figure 6s
#-------------------------------------------------------------------------------
# table(results_cox$population)
results_cox_plot<-results_cox %>% mutate(stars=case_when(p_adj<=0.001  ~ '***',
                                                    p_adj<=0.01  ~ '**',
                                                    p_adj<=0.05  ~ '*', 
                                                    T ~ ''),
                                          pstar=paste0(round(p_adj, 3), " ", stars)) %>% 
  filter(population %in% c('cohort-whole','controls', 'hpfs','controls-hpfs','nhs1','controls-nhs'))

results_cox_plot$population<-factor(results_cox_plot$population, 
                                    levels = c('cohort-whole','controls', 'hpfs','controls-hpfs','nhs1','controls-nhs'), 
                                    labels = c('HPFS+NHS','HPFS+NHS-controls', 'HPFS', 'HPFS-controls','NHS','NHS-controls'))

ggplot(data = results_cox_plot,
       aes(x = population, y =fct_reorder(metabolite_name, desc(class_broad)), fill = hr_adj)) + 
  geom_tile(color = "white", size = 0.1) + 
  scale_fill_gradientn(colours = hm.palette, limits = c(0.5,1.5))+
  #scale_fill_gradient2(midpoint = 1,low = '#3C8DAD',mid = "white",high = '#FF6767',limits=c(0.5,1.5)) +
  geom_text(aes(x = population, y =metabolite_name, label=stars), color="black", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Individual metabolites\nassociated with hip fracture")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 12, vjust = 0.5, 
                                   size = 12, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 1, 2, "cm"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 10, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  #new_scale_fill()+
  geom_tile(data = results_cox_plot %>% distinct(metabolite_name, class_broad),
            aes(x = 0, y = metabolite_name, col = factor(class_broad)),fill='white', 
            width = 0.1, height = 1, size=6, inherit.aes = FALSE) +
  scale_colour_lancet() +
  #scale_fill_manual(values = pathway_colors) +
  geom_text(data = results_cox_plot %>% distinct(metabolite_name, class_broad),
            aes(x = -0.2, y = metabolite_name, label = class_broad), color = "black", size = 3, 
            hjust = 1, inherit.aes = FALSE)+
  guides(fill=guide_colorbar(title = "adjusted HR",title.position = "top",nrow = 1, label.position = "bottom"),
         color = guide_legend(title = "Metabolite Class",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = "* P<0.05; ** P<0.01") 

# ggpubr::ggarrange(f_1, f_2, labels = c('3a', '3b'), ncol=1, nrow=2, heights = c(1,3), common.legend = F)
ggsave(filename = "figure 6s.pdf", 
       width = 9,
       height = 9,
       dpi = 600)


#-------------------------------------------------------------------------------
# (8). Figure 5s
#-------------------------------------------------------------------------------
# table(results_cox$population)
results_cox_plot<-results_cox %>% mutate(stars=case_when(p_adj<=0.001  ~ '***',
                                                         p_adj<=0.01  ~ '**',
                                                         p_adj<=0.05  ~ '*', 
                                                         T ~ ''),
                                         pstar=paste0(round(p_adj, 3), " ", stars)) %>% 
  filter(population %in% c('within10', 'exceed10','hpfs-within10','exceed10-hpfs','nhs-within10','exceed10-nhs'))

results_cox_plot$population<-factor(results_cox_plot$population, 
                                    levels = c('within10', 'exceed10','hpfs-within10','exceed10-hpfs','nhs-within10','exceed10-nhs'), 
                                    labels = c('HPFS+NHS within 10 years','HPFS+NHS beyond 10 years', 
                                               'HPFS within 10 years','HPFS beyond 10 years',
                                               'NHS within 10 years','NHS beyond 10 years'))

ggplot(data = results_cox_plot,
       aes(x = population, y =fct_reorder(metabolite_name, desc(class_broad)), fill = hr_adj)) + 
  geom_tile(color = "white", size = 0.1) + 
  scale_fill_gradientn(colours = hm.palette, limits = c(0.5,1.5))+
  #scale_fill_gradient2(midpoint = 1,low = '#3C8DAD',mid = "white",high = '#FF6767',limits=c(0.5,1.5)) +
  geom_text(aes(x = population, y =metabolite_name, label=stars), color="black", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("Individual metabolites\nassociated with hip fracture")+
  labs(y="")+
  theme(axis.text.x = element_text(angle = 12, vjust = 0.5, 
                                   size = 12, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 1, 2, "cm"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 10, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  #new_scale_fill()+
  geom_tile(data = results_cox_plot %>% distinct(metabolite_name, class_broad),
            aes(x = 0, y = metabolite_name, col = factor(class_broad)),fill='white', 
            width = 0.1, height = 1, size=6, inherit.aes = FALSE) +
  scale_colour_lancet() +
  #scale_fill_manual(values = pathway_colors) +
  geom_text(data = results_cox_plot %>% distinct(metabolite_name, class_broad),
            aes(x = -0.2, y = metabolite_name, label = class_broad), color = "black", size = 3, 
            hjust = 1, inherit.aes = FALSE)+
  guides(fill=guide_colorbar(title = "adjusted HR",title.position = "top",nrow = 1, label.position = "bottom"),
         color = guide_legend(title = "Metabolite Class",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = "* P<0.05; ** P<0.01") 

# ggpubr::ggarrange(f_1, f_2, labels = c('3a', '3b'), ncol=1, nrow=2, heights = c(1,3), common.legend = F)
ggsave(filename = "figure 5s.pdf", 
       width = 9,
       height = 9,
       dpi = 600)


#-------------------------------------------------------------------------------
# (9). Figure S1
#-------------------------------------------------------------------------------
z_metabolite_miss=print(naniar::miss_var_summary(e_imputed))
z_metabolite_miss$x=rank(z_metabolite_miss[["pct_miss"]],ties.method='first', na.last=1)
z_metabolite_miss<-z_metabolite_miss %>% arrange(x)

f_2=ggplot(data=z_metabolite_miss,aes(x, y=pct_miss))+ geom_point(shape=1) +
  labs(x='rank % missing', y='% missing', title='652 metabolites') + theme_bw () +
  theme(axis.title = element_text(size = 10,face = 'bold'))+
  geom_hline(yintercept = 20,linewidth = 0.7,lty = "dashed")

z_metabolite_miss<-z_metabolite_miss %>% mutate(prop_miss=case_when(pct_miss>=60 ~ '60~',
                                                                    pct_miss>=50 ~ '50-59.9',
                                                                    pct_miss>=40 ~ '40-49.9',
                                                                    pct_miss>=30 ~ '30-39.9',
                                                                    pct_miss>=20 ~ '20-29.9',
                                                                    pct_miss>=10 ~ '10-19.9',
                                                                    pct_miss>0 ~ '0.1-9.9',
                                                                    T ~ '0'),
                                                overall='Metabolites')

f_1=ggplot(data = z_metabolite_miss,aes(x=overall,fill=factor(prop_miss)))+
  geom_bar(stat='count', position = "fill", width = 1) + 
  scale_fill_nejm()+ theme_classic() +
  labs(x='', y='')+
  scale_y_continuous(expand = c(0,0.04)) + coord_flip() +
  theme(axis.text.y = element_text(size = 8,face = 'bold'),
        legend.position = 'top')+
  guides(fill = guide_legend(title = 'missing proportion, %')) 

ggpubr::ggarrange(f_1, f_2, nrow = 2, heights = c(1,1.5), labels = c('a', 'b'), align='v')

ggsave(filename = "figure 1s.pdf", 
       width = 7,
       height = 5,
       dpi = 600)

#-------------------------------------------------------------------------------
# (10). Figure S2
#-------------------------------------------------------------------------------
met_cor<-results_mv_main %>% filter(SUB_PATHWAY %in% c('Phosphatidylcholine (PC)','Phosphatidylethanolamine (PE)')) %>% 
  pull(met) %>% unique()
dat_cor=d %>% select(all_of(met_cor),"bxc", "osteocalcin", "inpth", "vitd25","sclerostin", 'calcium',
                     "dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti")

# print(naniar::miss_var_summary(dat_cor))
corr<- psych::corr.test(dat_cor, method = 'spearman')
r2=data.frame(corr$r[1:34,35:47])
r2$met=rownames(r2)
rr = tidyr::pivot_longer(r2, cols = c("bxc", "osteocalcin", "inpth", "vitd25", "sclerostin", "calcium", 
                                      "dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti"), 
                         names_to = "var", values_to = "value")

p2=data.frame(corr$p[1:34,35:47])
p2$met=rownames(p2)
pp = tidyr::pivot_longer(p2, cols = c("bxc", "osteocalcin", "inpth", "vitd25", "sclerostin", "calcium", 
                                      "dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti"), 
                         names_to = "var", values_to = "pvalue")
corr=cbind(rr,pp[,3])
corr<-merge(corr, f_match, by = 'met',all.x = TRUE, sort = TRUE)

corr$var <-factor(corr$var,levels = c("bxc", "osteocalcin", "inpth", "vitd25", "sclerostin", "calcium", 
                                      "dcalci", "scalc", "dvitdi", "svitd", "caffi", "alcoi", "tproti"),
                  labels = c("BetaXlaps", "osteocalcin", "PTH", "Vitamin D-25", "sclerostin", "calcium", 
                             "dietary calcium", "supplement calcium", "dietary VD-25", 
                             "supplement VD-25", "caffine intake", "alcohol intake", "protein intake"))
corr$p_group <- cut(corr$pvalue, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), 
                    labels = c("<0.0001", "<0.001", "<0.01", "<0.05", "<0.1",">0.1"))
corr$p_group <-factor(corr$p_group, levels = c(">0.1","<0.1","<0.05","<0.01","<0.001","<0.0001"))

corr<-corr %>% mutate(value_2=case_when(value<= -0.1 | value>= 0.1 ~ value,
                                        T ~ 0))
ggplot()+
  geom_point(data = corr[which(!is.na(corr$p_group)),],
             shape = 19, #col='black', # stroke = 0,
             size=2,
             aes(x = var, y = fct_reorder(METABOLITE, desc(SUB_PATHWAY)),col = value))+
  scale_color_gradient2(midpoint = 0,low = '#3C8DAD',mid = "white",high = '#FF6767',limits=c(-0.3, 0.3)) +
  geom_tile(data = corr[which(!is.na(corr$p_group)),] %>% distinct(METABOLITE, SUB_PATHWAY),
            aes(x = 0.5, y = METABOLITE, fill = factor(SUB_PATHWAY)),
            width = 0.3, height = 1, size=6, inherit.aes = FALSE) +
  scale_fill_lancet () +
  cowplot::theme_cowplot() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size=11,angle = 45, hjust = 1),
        axis.text.y = element_text(size=10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5,size=14))+
  labs(color = "Spearman's correlation", #size = "Spearman's p",
       fill = "Sub-pathway",
       #size='Correlation p-value',
       #title = "",
       title = "Correlation of metabolites with Covariates",
       x='',
       y='')

ggsave(filename = "figure 4s.pdf", 
       width = 12,
       height = 9,
       dpi = 600)

#-------------------------------------------------------------------------------
# (8). Regression results output
#-------------------------------------------------------------------------------
# Data 1 #
# names(results_mv_main)
results_mv_main_out <- results_mv_main %>% 
  mutate(ORCI = paste0(round(OR_1, 2), " (", round(lci_1, 2), "-", round(uci_1, 2), ")"),
         Pvalue = round(p_1, 2), FDR = round(p_adj_fdr, 2)) %>% 
  select("met", "n","ORCI", "Pvalue", "FDR", 
         "population", "SUB_PATHWAY", "METABOLITE", "HMDB")
openxlsx::write.xlsx(results_mv_main_out,"Data 1.xlsx")



# Data 2 #
# names(results_cox)
nes_sig_out=bind_rows(nes_sig_out_1, nes_sig_out_2)
openxlsx::write.xlsx(nes_sig_out,"Data 2.xlsx")



# Data 3 #
# names(results_cox)
results_cox_out <- results_cox %>% 
  mutate(HRCI = paste0(round(hr_adj, 2), " (", round(ll_adj, 2), "-", round(ul_adj, 2), ")"),
         Pvalue = round(p_adj, 2)) %>% 
  select("met", "n", "HRCI", "Pvalue", 
         "population", "class_broad", "metabolite_name")
openxlsx::write.xlsx(results_cox_out,"Data 3.xlsx")



