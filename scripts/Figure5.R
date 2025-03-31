###
#Manuscript title: Landscape of mobile genetic elements and their antibiotic resistance cargo in prokaryotic genomes

#The following code can be used for creating Figure 5 in the manuscript

#19-02-2021
#supriya.khedkar@embl.de
###

####load packages####

library(tidyverse)
library(cowplot)
library(viridis)
library(scales)

#### load data ####

arg_mge <- read_tsv("raw_data/mge_arg_final.txt", col_names = T)
arg_genome <- read_tsv("raw_data/genome_arg.txt", col_names = F)

m_arg_rm <- read_tsv("raw_data/mge_arg_resistance_mechanism_final.txt", col_names = T)
g_arg_rm <- read_tsv("raw_data/genome_arg_resistance_mechanism.txt", col_names = T)
mge_solitary <- read_tsv("processed_data/solitary_mge_bins_final.txt", col_names = T)
rm_list <- read_tsv("raw_data/resistance_mechanism.list", col_names = T)
glist <- read_tsv("raw_data/genome_status_supplementary_tableS2.txt", col_names = T)

#remove low and medium quality genomes#

glist_high <- glist %>% 
  filter(genome_quality == "high")

####Figure 5A Enrichment analysis of ARGs on MGE####

arg_mge <- arg_mge %>% 
dplyr::rename(abr_count = X3, total = X4) %>% 
mutate(mge = if_else(mge == "Cellular", "Genome", mge)) %>% 
mutate(island1 = island) %>% 
separate(island1, c("g1","g2","g3")) %>% 
mutate(genome = paste(g1, g2, sep = ".")) %>% 
select(-g1,-g2,-g3) %>% 
filter(genome %in% glist_high$genome) %>% 
select(-genome)

arg_genome <- arg_genome %>% 
  dplyr::rename(mge = X1, island = X2, abr_count = X3, total = X4) %>% 
  mutate(abr_count = ifelse(abr_count > 0, abr_count, 0)) 

arg_combined <- rbind (arg_mge,arg_genome)

arg_combined_gNames <-  arg_combined %>%
#  slice(1:1000) %>%
  mutate(island1 = island) %>%
  separate(island1, c("one","two","three"), sep = "\\.", extra = "merge") %>%
  unite(.,genome, c("one","two"), sep = ".") %>%
  select(-three)

#arg_combined_gNames <- read_tsv("processed_data/arg_mge_gNames.txt", col_names = T)
arg_mge_filtered <- arg_mge %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome) %>% 
  filter(.,total > 2) %>% 
  group_by(mge) %>% 
  summarise(abr_total = sum(abr_count), all = sum(total) - sum(abr_count))

arg_combined_filtered <- arg_combined_gNames %>% 
  filter(genomeID %in% glist_high$genome) %>% 
  select(-genomeID) %>% 
  filter(.,total > 2) %>% 
  group_by(mge) %>% 
  summarise(abr_total = sum(abr_count), all = sum(total) - sum(abr_count))

#enrichment analysis#
pre_arg_combined_filtered_m <- arg_combined_filtered %>% 
  dplyr::rename(abr = abr_total, non_abr = all) %>% 
  reshape2::melt() %>% 
  spread(mge,value)

allrec <- subset(colnames(pre_arg_combined_filtered_m), colnames(pre_arg_combined_filtered_m)!="variable") 
alldomains <- as.vector(pre_arg_combined_filtered_m$variable)
rownames(pre_arg_combined_filtered_m) <- as.vector(pre_arg_combined_filtered_m$variable)
arg_combined_filtered_m <- pre_arg_combined_filtered_m %>% select(-variable)

FE <- matrix(NA, nrow = length(alldomains), ncol = length(allrec),dimnames = list(alldomains, allrec))
OD <- matrix(NA, nrow = length(alldomains), ncol = length(allrec),dimnames = list(alldomains, allrec))
for (i in alldomains) {
  for(j in allrec) {
    print(paste0(i, "_", j))
    incog <- arg_combined_filtered_m [i,j]
    outcog <- sum(arg_combined_filtered_m [,j]) - incog
    notmge <- sum(arg_combined_filtered_m [i,]) - incog
    outcogoutmge <- sum(arg_combined_filtered_m ) - incog - outcog - notmge
    m <- matrix(data = c(incog, outcog, notmge, outcogoutmge), ncol = 2)
    test <- fisher.test(m,alternative = "greater")
    FE[i,j] <- test$p.value
    OD[i,j] <- test$estimate
  }
}

FE_adj <- matrix(p.adjust(p = FE, method = "bonferroni"), ncol = ncol(FE), nrow = nrow(FE), dimnames = list(rownames(FE), colnames(FE)))

OD_mat <- matrix(OD, ncol = ncol(OD), nrow = nrow(OD), dimnames = list(rownames(OD), colnames(OD)))

od_melt <- OD_mat %>% 
  reshape2::melt() %>% 
  dplyr::rename(resistance = Var1, mge = Var2, score = value)

fe_melt <- FE_adj %>% 
  reshape2::melt() %>% 
  dplyr::rename(resistance = Var1, mge = Var2, adj_pval = value)

merge_od_pval <- left_join(od_melt,fe_melt, by = c("resistance", "mge")) %>% 
  filter(.,resistance == "abr") %>% 
  mutate(signP=ifelse(adj_pval>0.05,"p value > 0.05","p value < 0.05"))

combine <- left_join(merge_od_pval,arg_combined_filtered, by = "mge")
abr_nplot <- ggplot(combine, aes(x = reorder(mge,-abr_total,mean), y = log10(abr_total), fill = as.factor(signP))) + geom_bar(stat = "identity") + scale_fill_manual(values=c('#E69F00','#999999')) + geom_text(data = combine, aes(x = mge, y = log10(abr_total), label = abr_total, angle = 90)) 
abr_nplot + theme_cowplot(font_size = 20) + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs (x = "", y = "#ABR genes" ,fill = "") 


####Figure 5B Association of ARG resistance mechanisms####

mge_solitary <- mge_solitary %>% 
  dplyr::rename(IS_Tn = X1, Phage = X2, Phage_like = X3, CE = X4, Integron = X5, MI = X6,	Hotspot	= X7, UC = X8,	Cellular = X9, island = X10, island_size = X11, prot_count = X12,	phage_count = X13, CONJ_T4SS = X14,	mgeR = X15)

m_arg_rm <- m_arg_rm %>% 
  dplyr::rename(island1 = ptn, prot = island) %>% 
  dplyr::rename(island = island1) %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome)

m_arg_rm %>% pull(prot) %>%  n_distinct()
g_arg_rm <- g_arg_rm %>% 
  dplyr::rename(island1 = ptn, prot = island) %>% 
  dplyr::rename(island = island1) %>% 
  mutate(mge = "Genome") %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome)

mge_non_nested <- mge_solitary %>% 
  select(1:11,15) %>% 
  gather(mge, mge_pa, 1:9) %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome) %>% 
  filter(mge_pa ==1) %>% 
  filter(.,!grepl("UC", mge)) %>% 
  filter(.,!grepl("Hotspot", mge)) %>% 
  mutate(mgeRn = str_replace_all(mgeR,"_","")) %>% 
  filter(!str_detect(mgeRn, '[:alnum:] &{1,}')) %>% 
  select(-mgeRn) %>% 
  mutate(mgeR = str_replace_all(mgeR,"&","")) %>% 
  select(-mgeR,-mge_pa)

m_arg_rm_final <- left_join(m_arg_rm,mge_non_nested, by="island") %>% 
drop_na(mge)

#
rm_list <- rm_list %>% 
  mutate(abr_categories = str_replace_all(abr_categories,"\\s","_"))

m_arg_rm_final <- m_arg_rm_final %>% 
  select(-confers_resistance_to,-island_size)
g_arg_rm_final <- g_arg_rm %>% 
  select(-confers_resistance_to)  

together_rm <- rbind(m_arg_rm_final,g_arg_rm_final) %>% 
  mutate(mge = ifelse(mge == "Cellular", "Genome", mge)) %>% 
  mutate(participates_in_1 = str_replace_all(participates_in,"^ARO:(\\d+)\\s!\\s",""), is_a_1 = str_replace_all(is_a,"^ARO:(\\d+)\\s!\\s","")) %>% 
  mutate(participates_in_1 = str_replace_all(participates_in_1,"\\s","_"), is_a_1 = str_replace_all(is_a_1,"\\s","_")) %>% 
  mutate(comb =if_else(participates_in_1 %in% rm_list$abr_categories, participates_in_1, is_a_1)) %>% 
  filter(!is.na(comb)) %>%
  filter(comb != "reduced_permeability_to_antibiotic" | comb != "GO:0015291_!_secondary_active_transmembrane_transporter_activity") %>% 
  select(island,mge,comb) %>%
#  group_by(mge,island,comb) %>%
  unique() %>%
#  ungroup() %>%
  group_by(mge,comb) %>%
  summarise(count = n()) %>% 
  filter(comb != "GO:0015291_!_secondary_active_transmembrane_transporter_activity")

pre_arg_rm <- together_rm %>%
  filter(comb != "hydrolysis_of_antibiotic_conferring_resistance") %>%
  dplyr::rename(abr_categories = comb) %>% 
  spread(mge, count,fill = 0)

m_pre_arg <- pre_arg_rm %>% 
reshape2::melt()

arg_rm <- as.data.frame(pre_arg_rm)
allMges <- subset(colnames(arg_rm), colnames(arg_rm)!="abr_categories") 
allClasses <- as.vector(arg_rm$abr_categories)
rownames(arg_rm) <- as.vector(arg_rm$abr_categories)
arg_rm <- arg_rm %>% select(-abr_categories)

##fishers exact
FE_arg <- matrix(NA, nrow = length(allClasses), ncol = length(allMges),dimnames = list(allClasses, allMges))
for (i in allClasses) {
  for(j in allMges) {
    print(paste0(i, "_", j))
    inarg <- arg_rm[i,j]
    outarg <- sum(arg_rm[,j]) - inarg
    notmge <- sum(arg_rm[i,]) - inarg
    outargoutmge <- sum(arg_rm) - inarg - outarg - notmge
    m <- matrix(data = c(inarg, notmge,outarg, outargoutmge), ncol = 2)
    test <- fisher.test(m, alternative = "greater")
    FE_arg[i,j] <- test$p.value
  }
}

FE_arg_adj <- matrix(p.adjust(p = FE_arg, method = "bonferroni"), ncol = ncol(FE_arg), nrow = nrow(FE_arg), dimnames = list(rownames(FE_arg), colnames(FE_arg)))

FE_arg_adj_t <- as.tibble(FE_arg_adj) %>% 
  mutate(abr = rownames(FE_arg_adj))
m_FE_arg_adj_t <- FE_arg_adj_t %>% 
  gather(mge,p_adj,1:7) %>% 
  mutate(signP=ifelse(p_adj>0.05,"","*"))


m_pre_arg_n <- m_pre_arg %>% 
  dplyr::rename(abr = abr_categories,mge = variable)
m_FE_arg_adj_t_rs <- left_join(m_FE_arg_adj_t, m_pre_arg_n, by = c("abr", "mge")) 

arg_hmap <- ggplot(m_FE_arg_adj_t_rs,aes(x = reorder(mge,value, sum),y = reorder(abr,value,sum), fill = log10(value))) + geom_tile(colour = "black") + geom_text(aes(label = signP), colour = "white", size = 8) + scale_fill_gradientn(colours = viridis(10), values = rescale(c(0,3,4,6)), na.value = "grey80")

arg_hmap + theme_cowplot(font_size=20) + theme(axis.text.x=element_text(angle = 90, hjust = 1))  + labs(x = "", y = "", fill = "#ABR genes")
