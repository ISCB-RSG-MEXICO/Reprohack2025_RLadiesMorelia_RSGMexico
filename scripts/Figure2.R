###
#Manuscript title: Landscape of mobile genetic elements and their antibiotic resistance cargo in prokaryotic genomes

#The following code can be used for creating Figure 2 in the manuscript

#19-02-2021
#supriya.khedkar@embl.de
###

####load packages####

library(tidyverse)
library(cowplot)
library(ggpubr)
library(ggforce)


####load data####
setwd("/g/scb2/bork/khedkar/jumpAR/promge_manuscript")
mge_pg <- read_tsv("processed_data/mge_bins_per_genome_final.txt", col_names = T)
mge_solitary <- read_tsv("processed_data/solitary_mge_bins_final.txt", col_names = F)
rec_class <- read_tsv("raw_data/recombinase.list", col_names = F)
glist <- read_tsv("raw_data/genome_status_supplementary_tableS2.txt", col_names = T)

#color scheme#

colc <- c("#D55E00", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2","#CECCCC")
names(colc) <- c("IS_Tn", "Phage", "Phage_like", "CE", "Integron", "MI", "Cellular") 

#remove low and medium quality genomes#

glist_high <- glist %>% 
  filter(genome_quality == "high")

####Figure 2A - alluvial plot####

rec_class <- rec_class %>% dplyr::rename(class = X1, mgeR = X2)

mge_solitary <- mge_solitary %>% 
  dplyr::rename(IS_Tn = X1, Phage = X2, Phage_like = X3, CE = X4, Integron = X5, MI = X6,	Hotspot	= X7, UC = X8,	Cellular = X9, island = X10, island_size = X11, prot_count = X12,	phage_count = X13, CONJ_T4SS = X14,	mgeR = X15)

mge_solitary_melted <- mge_solitary %>% 
  select(1:11,15) %>% 
  gather(mge, mge_pa, 1:9) %>% 
  filter(mge_pa ==1) %>% 
  filter(.,!grepl("UC", mge)) %>% 
  filter(.,!grepl("Hotspot", mge)) %>% 
  mutate(mgeRn = str_replace_all(mgeR,"_","")) %>% 
  filter(!str_detect(mgeRn, '[:alnum:] &{1,}')) %>% 
  select(-mgeRn) %>% 
  mutate(mgeR = str_replace_all(mgeR,"&",""))

mge_solitary_melted_dw <- mge_solitary_melted %>% 
  select(-mgeR,-mge_pa) %>% 
  mutate(island1 = island) %>% 
  separate(island1, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>% 
  select(-genome)

#write.table(mge_solitary_melted_dw,file="processed_data/mge_bins_final_solitary_collapsed.txt", sep = "\t", row.names = F, col.names = T, quote = F)

mge_solitary_rclass_all <- left_join(mge_solitary_melted,rec_class, by = "mgeR") %>% 
  select(-island_size, -mge_pa, -mgeR) %>% 
  separate(island, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  group_by(class,mge) %>% 
  summarise(count = n())

#high quality genomes
mge_solitary_rclass <- left_join(mge_solitary_melted,rec_class, by = "mgeR") %>% 
  select(-island_size, -mge_pa, -mgeR) %>% 
  separate(island, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome) %>%
  group_by(class,mge) %>% 
  summarise(count = n())


#select recombinase family
colclass <- c("CDCCCC", "CDCCCC", "CDCCCC", "CDCCCC", "CDCCCC")
names(colclass) <- c( "cas", "dde","huh", "ser", "tyr")
colall <- c(colclass, colc)

order_class <- c( "cas", "dde","huh", "ser", "tyr")
order_mge <- c("IS_Tn", "Phage_like", "Phage", "CE", "MI", "Integron", "Cellular")
mge_solitary_rclass$mge <- factor(mge_solitary_rclass$mge, levels = order_mge)
mge_solitary_rclass_p <- gather_set_data(mge_solitary_rclass, 1:2)
mge_solitary_rclass_p$y <- factor(mge_solitary_rclass_p$y, levels = c(order_class, order_mge))
mge_solitary_rclass_p <- mge_solitary_rclass_p %>% add_column(col = colall[match(.$y, names(colall))])

rclass_mge_alluvial <- ggplot(mge_solitary_rclass_p, aes(x, id = id, split = y, value = count)) +
  geom_parallel_sets(aes(fill = mge), color = "black", lwd = 0.2, axis.width = 0.26) +
  geom_parallel_sets_axes(axis.width = 0.22, fill = "grey80") +
  geom_parallel_sets_labels(
    color = 'black',
    size = 12/.pt,
    angle = 0
  ) +
  scale_x_discrete(
    name = NULL,
    expand = c(0, 0.12)
  ) +
  scale_y_continuous(breaks = NULL, expand = c(0.1, 0)) +
  scale_fill_manual(
  values = colc,
    guide = "none"
  ) +
  labs(fill = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(14, 1.5, 2, 1.5)
  )
rclass_mge_alluvial

####Figure 2B - barplot####

mge_pg_melt_all <- mge_pg %>% 
  reshape2::melt() %>% 
  select(1,7,8) %>% 
  group_by(variable) %>% 
  summarise(total = sum(value)) %>% 
  filter(.,!grepl("Hotspot",variable))

#high quality genomes
mge_pg_melt <- mge_pg %>% 
  filter(Genome %in% glist_high$genome) %>% 
  reshape2::melt() %>% 
  select(1,7,8) %>% 
  group_by(variable) %>% 
  summarise(total = sum(value)) %>% 
  filter(.,!grepl("Hotspot",variable))

barplot_2b <- ggplot(mge_pg_melt,aes(x = reorder(variable,total,sum), y = total, fill = variable)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = total)) + 
  coord_flip() +
  scale_fill_manual("MGE", values = colc, guide = F) +
  theme_cowplot() + 
  labs(x ="", y = "counts")
barplot_2b

####Figure 2B - donut chart####
mge_pg_relative <- mge_pg_melt %>% 
  filter(.,!grepl("Cellular",variable)) %>% 
  mutate(rel = round((total/sum(total))*100,digits = 2), labs = paste0(variable, " (", rel, "%)")) %>%
  arrange(total)

donutchart_2b <- ggdonutchart(mge_pg_relative, "rel", label = "labs",
             lab.pos = "in",
             fill = "variable", color = "white",
             palette = c("#D55E00", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2","#CECCCC"))
donutchart_2b

####Figure 2C - boxplot####

mge_solitary_length_all <- mge_solitary %>% 
  select(1:11) %>% 
  gather(mge, mge_pa, 1:9) %>% 
  filter(mge_pa ==1) %>% 
  filter(.,!grepl("UC", mge)) %>% 
  filter(.,!grepl("Hotspot", mge))

#high quality genomes
mge_solitary_length <- mge_solitary %>% 
  select(1:11) %>% 
  gather(mge, mge_pa, 1:9) %>% 
  filter(mge_pa ==1) %>% 
  filter(.,!grepl("UC", mge)) %>% 
  filter(.,!grepl("Hotspot", mge)) %>% 
  separate(island, c("g1","g2","g3")) %>% 
  mutate(genome = paste(g1, g2, sep = ".")) %>% 
  select(-g1,-g2,-g3) %>% 
  filter(genome %in% glist_high$genome)
  
mge_length_boxplot <- ggplot(mge_solitary_length, aes(x = reorder(mge,-island_size,median), y =island_size, fill = mge)) + 
  geom_boxplot(outlier.shape = NA, notch = F, lwd = 1) + 
  scale_fill_manual("MGE", values = colc, guide = F) + 
  scale_y_continuous(limits = quantile(mge_solitary_length$island_size, c(0.1, 0.9))) + 
  coord_flip()

mge_length_boxplot + labs (y = "length (bp)", x = "") + theme_cowplot(font_size = 20)
