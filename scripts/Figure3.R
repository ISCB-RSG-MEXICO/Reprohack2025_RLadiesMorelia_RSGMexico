###
#Manuscript title: Landscape of mobile genetic elements and their antibiotic resistance cargo in prokaryotic genomes

#The following code can be used for creating Figure 3 in the manuscript

#19-02-2021
#supriya.khedkar@embl.de
###

####load packages####

library(tidyverse)
library(cowplot)
library(phytools)
library(viridisLite)
library(viridis)
library(scales)

####load data####

tax <- read_tsv("raw_data/species_with_atleast_2genomes.list", col_names=F)
db <-read_tsv("processed_data/mge_bins_per_genome_final.txt", col_names = T)
gs <-read_tsv("raw_data/genome_size.txt", col_names = T)
class_tree <- read.tree("raw_data/progenomes2_class_tree.nwk")
glist <- read_tsv("raw_data/genome_status_supplementary_tableS2.txt", col_names = T)

#remove low and medium quality genomes#

glist_high <- glist %>% 
  filter(genome_quality == "high")

#color scheme#

colc <- c("#D55E00", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2","#CECCCC")
names(colc) <- c("IS_Tn", "Phage", "Phage_like", "CE", "Integron", "MI", "Cellular") 

#genome and species stats
tax %>% 
  filter(X2 %in% glist_high$genome) %>% 
  pull(X1) %>% 
  n_distinct()
#76,228 genomes 3207species

####Figure 3B####
colnames(tax) <- c("specI","genomeID","kingdom","phylum","class","genus")

mdb <- db %>% reshape2::melt() 
mdb <- mdb %>% dplyr::rename(mge = variable, count = value, genomeID = 'Genome') %>% 
  filter(genomeID %in% glist_high$genome)

##Format data 
#create a data.frame containing for each specI for each MGE: genome counts (genomeCnt), 
#total counts (cnt_tot), 
#average counts across genomes (avg_cnt), 
#number of genomes where the MGE was present (pa), 
#the fraction of genomes with the MGE present (frac) plus taxonomic information

##generate data.frame with all combinations of genomes and MGEs with tax infos and MGE counts
mdb_cnt <- mdb %>% select(1,7,8)
all_mge <- unique(mdb$mge)

##create a big data frame with all combinations of bacteria and all mge.
db_tax <-NULL       
for(i in 1:length(all_mge)){
  db_tax_sing <- tax %>% add_column(mge=all_mge[i])
  db_tax <- rbind(db_tax, db_tax_sing)
}
head(db_tax)
db_cnt_all <- db_tax %>% 
  filter(genomeID %in% glist_high$genome) %>% 
  left_join(., mdb_cnt, by=c("genomeID", "mge"))
db_cnt_all[is.na(db_cnt_all)] <-0

#number of genomes  per species per class
db_cnt_all %>% 
  select(specI,genomeID, class) %>% 
  unique() %>% 
  group_by(class) %>% 
  summarise(count = n()) %>% 
  filter(count > 9) %>% 
  View()

#create presence absence data frame
## add presence absence information
db_pa_all <- db_cnt_all %>% mutate(presAbs=ifelse(count>0, 1,0))


##group all genomes by specI 
##summarise: avg_cnt= average count per mge per specI, frac= fraction of genomes per specI with each mge
db_specI <- db_pa_all %>% group_by(specI, mge) %>% summarise(genomeCnt=n(), avg_cnt=mean(count), pa=sum(presAbs), cnt_tot=sum(count)) %>% mutate(frac=pa/genomeCnt) %>% left_join(., tax, by="specI") %>% select(-genomeID) %>% unique(.)


## summarise: avg_cnt= average count per mge per genome, frac= fraction of genomes with each mge
db_genome <- db_pa_all %>% group_by(genomeID, mge) %>% summarise(genomeCnt=n(), avg_cnt=mean(count), pa=sum(presAbs), cnt_tot=sum(count)) %>% mutate(frac=pa/genomeCnt) %>% left_join(., tax, by="genomeID") %>% unique(.)


##create seperate data frames for mge count and mge frac data
#average Counts of each MGE per specI where MGEs are columns and specIs rows
all_cnt_specI <- db_specI %>% select(specI, mge, avg_cnt)
all_frac_specI <- db_specI %>% select(specI, mge, frac)
all_cnt_genome <- db_genome %>% select(genomeID, mge, avg_cnt)


##Analysis per specI by class
# Mann-Whitney test
#use Mann-Whitney test to analyze associations of MGEs on class level  
#Analysis is based on the fraction of genomes per SpecI that contains the MGE to avoid sampling-bias

# input tab 
all_avg_cnt_tax_class <-  db_specI %>% select(specI, mge, avg_cnt, class) %>% filter(.,!grepl("Hotspot",mge)) %>% filter(.,!grepl("Cellular",mge)) 
allMges <- unique(all_avg_cnt_tax_class$mge)

# filter for classes with at least 10 genomes
allClasses <- all_avg_cnt_tax_class %>% filter(mge==all_avg_cnt_tax_class$mge[1]) %>% group_by(class) %>% summarise(cntClass=n()) %>% filter(cntClass>9)

all_avg_cnt_tax_class_sel_pre <- all_avg_cnt_tax_class %>% filter(class %in% allClasses$class)

###genome size normalisation###
gs_int <- gs %>% group_by(SpecI_id_v3) %>% summarise(avg_gs = mean(ProteinGeneCounts)) %>% dplyr::rename(specI = SpecI_id_v3)

all_avg_cnt_tax_class_sel <- left_join(all_avg_cnt_tax_class_sel_pre,gs_int,by = "specI") %>% mutate(norm_count = avg_cnt/avg_gs) %>% dplyr::rename(count = avg_cnt, avg_cnt = norm_count)

# run mann whitney test
MW_all <- NULL
for(i in 1:length(allMges)){
  mgeX <- allMges[i]
  mgeDat <- all_avg_cnt_tax_class_sel %>% filter(mge == mgeX)
  MW <- sapply(seq_along(allClasses$class), function(j){
    class_sel <-allClasses$class[j]
    ingroup <- mgeDat %>% filter(class == class_sel)
    outgroup <- mgeDat %>% filter(class != class_sel)
    MWout <- wilcox.test(ingroup$avg_cnt,outgroup$avg_cnt, alternative = "greater")
    out <-c(class_sel ,MWout$p.value)
    return(out)
  })
  MWOut <- as.data.frame(t(MW[2,]))
  colnames(MWOut) <- MW[1,]
  rownames(MWOut) <- mgeX
  MW_all <- rbind(MW_all, MWOut) 
}

MW_all_mod <- as.data.frame(t(MW_all))

##adjust p_values
MW_all_adj <- NULL
for(k in 1:ncol(MW_all_mod)){
  MW_all_mod2 <- as.matrix(MW_all_mod)
  pVec <- c(MW_all_mod2[,k])
  adjP <- p.adjust(pVec, "BH")
  MW_all_adj <- cbind(MW_all_adj, adjP)
}

colnames(MW_all_adj) <- colnames(MW_all_mod)

MW_all_fin <- as.data.frame(MW_all_adj)

## binary tab with pVal < 0.1 --> 1 else 0
MW_all_stat <- NULL
for(l in 1:ncol(MW_all_fin)){
  MW_test <- MW_all_fin %>% mutate(signP=ifelse(MW_all_fin[,l]>0.1,0,1))
  MW_all_stat <- cbind(MW_all_stat, MW_test$signP)
}

MW_all_stat <- as.data.frame(MW_all_stat)
colnames(MW_all_stat) <- colnames(MW_all_fin)
rownames(MW_all_stat) <- rownames(MW_all_fin)
MW_all_stat_sum <- MW_all_stat %>% mutate(sumSign=rowSums(.))

MW_all_stat2 <- cbind(rownames(MW_all_stat), MW_all_stat)
colnames(MW_all_stat2)[1] <- "class"
MW_all_fin2 <- cbind(rownames(MW_all_fin), MW_all_fin)
colnames(MW_all_fin2)[1] <- "class"

melted_specI_class <- reshape2::melt(MW_all_fin2)
melted_specI_class1 <- melted_specI_class %>% mutate(signP=ifelse(value>0.05,"","*")) 
all_avg_cnt_tax_class_sel1 <- all_avg_cnt_tax_class_sel %>% group_by(class,mge) %>% summarise(avg_class = mean(count)) %>% dplyr::rename(variable = mge)

plot_table_w_g9 <- left_join(melted_specI_class1, all_avg_cnt_tax_class_sel1, by = c("class", "variable"))
plot_table <- left_join(melted_specI_class1, all_avg_cnt_tax_class_sel1, by = c("class", "variable"))

##with g > 9 (presence of at least 10 genomes)
specI_class_hmap_w_g9 <- ggplot(plot_table_w_g9,aes(y = reorder(variable,avg_class,sum), x = reorder(class,avg_class,sum), fill = avg_class)) + geom_tile(colour = "black") + geom_text(aes(label = signP), colour = "white") + scale_fill_gradientn(colours = alpha(viridis(10),0.75), values = rescale(c(0,10,100,500)), na.value = "grey80")
specI_class_hmap_w_g9 +  theme_cowplot() + theme(axis.text.x=element_text(angle=90,hjust = 1, vjust = 0.5)) + labs(x = "", y = "", fill = "Average counts per species")

##make phyloheatmap##
normalize <- function(x) {
  x/max(x)
}

phylo_heatmap_mat <- all_avg_cnt_tax_class_sel1 %>% 
  filter(!grepl("NA ", class)) %>% 
  spread(variable, avg_class) %>% 
  as_tibble(column_to_rownames(var = "class")) %>%
  select(class, IS_Tn, Phage, Phage_like, CE, Integron, MI) %>% 
  column_to_rownames(var = "class")


class_tree_w_g9 <- keep.tip(class_tree, tip = rownames(phylo_heatmap_mat))

##tweak colors in phyloheatmap for skewed data##
small_value <- unique(sort(unlist(phylo_heatmap_mat)))[2]
small_value <- small_value - small_value %% 0.001
hist(unlist(phylo_heatmap_mat), nclass = 50)
heatmap_breaks <- c(0,small_value,1,2,3,4,seq(5,55,10))
hist(unlist(phylo_heatmap_mat), breaks = heatmap_breaks, freq = T)
plot_colors <- c("white",colorRampPalette((viridis(10)), bias = 5)(length(heatmap_breaks)-2))

sig_mat <- plot_table_w_g9 %>% 
  select(class, variable, signP)

sig_mat <- sig_mat %>% 
  filter(!grepl("NA ", class)) %>%
  spread(variable, signP) %>% 
  column_to_rownames(var = "class")

sig_mat <- sig_mat[class_tree_w_g9$tip.label,colnames(phylo_heatmap_mat)]

##custom functions to extract x and y co-ordinates from the phylo heatmap and open up interactive legend plotting
phylo.heatmap.coords <- function (tree, X, fsize = 1, colors = NULL, standardize = FALSE, ...) 
{
  if (length(fsize) != 3) 
    fsize <- rep(fsize, 3)
  if (hasArg(legend)) 
    legend <- list(...)$legend
  else legend <- TRUE
  if (hasArg(labels)) 
    labels <- list(...)$labels
  else labels <- TRUE
  if (hasArg(split)) 
    split <- list(...)$split
  else split <- c(0.5, 0.5)
  split <- split/sum(split)
  if (is.null(colnames(X))) 
    colnames(X) <- paste("var", 1:ncol(X), sep = "")
  if (standardize) {
    sd <- apply(X, 2, function(x) sqrt(var(x, na.rm = TRUE)))
    X <- (X - matrix(rep(1, Ntip(tree)), Ntip(tree), 1) %*% 
            colMeans(X, na.rm = TRUE))/(matrix(rep(1, Ntip(tree)), 
                                               Ntip(tree), 1) %*% sd)
  }
  if (hasArg(xlim)) 
    xlim <- list(...)$xlim
  else xlim <- c(-0.5, (2 - 0.5) * split[2]/split[1] + 0.5)
  if (hasArg(ylim)) 
    ylim <- list(...)$ylim
  else ylim <- if (legend) 
    c(if (standardize) -0.15 else -0.1, if (labels) 1.1 else 1)
  else c(0, if (labels) 1.1 else 1)
  if (hasArg(mar)) 
    mar <- list(...)$mar
  else mar <- rep(1.1, 4)
  if (is.null(colors)) 
    colors <- heat.colors(n = 20)[20:1]
  if (hasArg(grid)) 
    add.grid <- list(...)$grid
  else add.grid <- FALSE
  cw <- untangle(tree, "read.tree")
  plot.new()
  par(mar = mar)
  plot.window(xlim = xlim, ylim = ylim)
  h <- phylogram(cw, fsize = fsize[1], ...)
  START <- h + 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - 
                        h)/(ncol(X) - 1) + 0.5 * strwidth("W") * fsize[1]
  END <- (2 - 0.5) * split[2]/split[1] + 0.5 - 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - START)/(ncol(X) - 1)
  nTips <- length(tree$tip.label)
  y <- c(-1/(2 * (nTips - 1)), seq(0, 1, length = nTips) + 2/(2 * (nTips - 1)))
  x <- seq(START, END, by = (END - START)/(ncol(X) - 1))
  list(xx = x, yy = y)
}
environment(phylo.heatmap.coords) <- environment(phylo.heatmap)

phylo.heatmap.legendmod <- function (tree, X, fsize = 1, colors = NULL, standardize = FALSE, ...) 
{
  if (length(fsize) != 3) 
    fsize <- rep(fsize, 3)
  if (hasArg(legend)) 
    legend <- list(...)$legend
  else legend <- TRUE
  if (hasArg(labels)) 
    labels <- list(...)$labels
  else labels <- TRUE
  if (hasArg(split)) 
    split <- list(...)$split
  else split <- c(0.5, 0.5)
  split <- split/sum(split)
  if (is.null(colnames(X))) 
    colnames(X) <- paste("var", 1:ncol(X), sep = "")
  if (standardize) {
    sd <- apply(X, 2, function(x) sqrt(var(x, na.rm = TRUE)))
    X <- (X - matrix(rep(1, Ntip(tree)), Ntip(tree), 1) %*% 
            colMeans(X, na.rm = TRUE))/(matrix(rep(1, Ntip(tree)), Ntip(tree), 1) %*% sd)
  }
  if (hasArg(xlim)) 
    xlim <- list(...)$xlim
  else xlim <- c(-0.5, (2 - 0.5) * split[2]/split[1] + 0.5)
  if (hasArg(ylim)) 
    ylim <- list(...)$ylim
  else ylim <- if (legend) 
    c(if (standardize) -0.15 else -0.1, if (labels) 1.1 else 1)
  else c(0, if (labels) 1.1 else 1)
  if (hasArg(mar)) 
    mar <- list(...)$mar
  else mar <- rep(1.1, 4)
  if (is.null(colors)) 
    colors <- heat.colors(n = 20)[20:1]
  if (hasArg(grid)) 
    add.grid <- list(...)$grid
  else add.grid <- FALSE
  cw <- untangle(tree, "read.tree")
  plot.new()
  par(mar = mar)
  plot.window(xlim = xlim, ylim = ylim)
  h <- phylogram(cw, fsize = fsize[1], ...)
  START <- h + 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - h)/(ncol(X) - 1) + 0.5 * strwidth("W") * fsize[1]
  END <- (2 - 0.5) * split[2]/split[1] + 0.5 - 1/2 * ((2 - 0.5) * split[2]/split[1] + 0.5 - START)/(ncol(X) - 1)
  X <- X[cw$tip.label, ]
  image(x = seq(START, END, by = (END - START)/(ncol(X) - 1)), 
        z = t(X[cw$tip.label, ]), add = TRUE, col = colors, right = T, ...)
  if (add.grid) {
    dx <- (END - START)/(ncol(X) - 1)
    x <- seq(START - dx/2, END + dx/2, by = dx)
    nTips <- length(tree$tip.label)
    y <- c(-1/(2 * (nTips - 1)), seq(0, 1, length = nTips) + 
             1/(2 * (nTips - 1)))
    segments(x, y[1], x, y[length(y)])
    segments(x[1], y, x[length(x)], y)
  }
  if (legend) 
    add.color.bar(leg = END - START, cols = colors, lims = range(X, na.rm = TRUE), title = if (standardize) "standardized value" else "value", subtitle = if (standardize) "SD units" else "", prompt = F, digits = if (max(abs(X),na.rm = TRUE) < 1) round(log10(1/max(abs(X), na.rm = TRUE))) + 1 else 2, fsize = fsize[3], lwd = 20, outline = T, x = START, y = -1/(2 * (Ntip(cw) - 1)) - 3 * fsize[3] * strheight("W"))
  if (labels) 
    text(x = seq(START, END, by = (END - START)/(ncol(X) - 1)), y = rep(1 + 1/(2 * (Ntip(cw) - 1)) + 0.4 * fsize[2] * strwidth("I"), ncol(X)), colnames(X), srt = 70, adj = c(0,0.5), cex = fsize[2])
  if (any(is.na(X))) {
    ii <- which(is.na(X), arr.ind = TRUE)
    x.na <- seq(START, END, by = (END - START)/(ncol(X) - 1))[ii[, 2]]
    y.na <- seq(0, 1, by = 1/(nrow(X) - 1))[ii[, 1]]
    for (i in 1:length(x.na)) {
      xx <- x.na[i] + c(1/2, -1/2) * (END - START)/(ncol(X) - 1)
      yy <- y.na[i] + c(-1/2, 1/2) * 1/(nrow(X) - 1)
      lines(xx, yy)
    }
  }
}
environment(phylo.heatmap.legendmod) <- environment(phylo.heatmap)

##extract xx and yy from phylo.heatmap custom function
xx_yy <- phylo.heatmap.coords(class_tree_w_g9, phylo_heatmap_mat, fsize = c(0.8, 0.9, 0.7), colors = plot_colors, grid = T, split = c(0.7, 0.3), lwd = 1, breaks = heatmap_breaks, mar = c(1.2,1.2,1.2,1.2))

##plot##
phylo.heatmap.legendmod(class_tree_w_g9, phylo_heatmap_mat, fsize = c(0.8, 0.9, 0.7), colors = plot_colors, grid = T, split = c(0.7, 0.3), lwd = 1, breaks = heatmap_breaks, mar = c(1.2,1.2,1.2,1.2))
for(i in 1:nrow(sig_mat)) {
  for(j in 1:ncol(sig_mat)) {
    if(sig_mat[i,j]=="*") {
      rnm <- rownames(sig_mat)[i]
      cnm <- colnames(sig_mat)[j]
      points(xx_yy$xx[j], xx_yy$yy[i], cex = 1.5, col = if(phylo_heatmap_mat[rnm,cnm] < heatmap_breaks[3]) "white" else "black", pch = "*")
    }
  }
}


####Figure 3A - barplot####
allclass_all <- all_avg_cnt_tax_class %>% filter(mge==all_avg_cnt_tax_class$mge[1]) %>% group_by(class) %>% summarise(cntclass=n()) %>% filter(.,cntclass > 9)

all_avg_cnt_tax_class_sel_all <- all_avg_cnt_tax_class %>% filter(.,!grepl("Cellular",mge)) %>% filter(.,!grepl("Hotspot",mge)) %>% filter(class %in% allclass_all$class) %>% group_by(mge,class) %>% summarise(avg_class = mean(avg_cnt))

all_avg_cnt_tax_class_sel_all_mge_all <- all_avg_cnt_tax_class_sel_all %>% group_by(class) %>% summarise(total_mge = sum(avg_class))

all_avg_cnt_tax_class_sel_4bar_all <- left_join(all_avg_cnt_tax_class_sel_all, all_avg_cnt_tax_class_sel_all_mge_all) %>% group_by(class, mge, total_mge) %>% summarise(frac = avg_class/total_mge)

barplot_mat <- all_avg_cnt_tax_class_sel_4bar_all %>% filter(.,!grepl("NA ",class))
barplot_mat$class <- factor(barplot_mat$class, levels = class_tree_w_g9$tip.label)

barplot_mat_for_n <- all_avg_cnt_tax_class_sel_all_mge_all %>% filter(.,!grepl("NA ",class))
barplot_mat_for_n$class <- factor(barplot_mat_for_n$class, levels = class_tree_w_g9$tip.label)

##plot##
barplot_specI_class_count_all <- ggplot(barplot_mat, aes(y = class, x = frac, fill = mge)) +  
  geom_barh(stat="identity",color = "grey60") + 
  scale_fill_manual("MGE", values = colc) 
barplot_specI_class_count_all +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0), limits = c(0, 1.1)) +
  geom_text(data = barplot_mat_for_n, mapping = aes(y = class, x = 1.05, label = paste("",round(total_mge))), inherit.aes = F, col = "black", size = 6) + 
  theme_cowplot(font_size = 20) + 
  theme(axis.text.x=element_text(angle = 0), legend.position = "bottom") +
  theme(axis.text.y=element_text(hjust=0,vjust=0.5)) + 
  labs(y = "", x="relative proportion")
