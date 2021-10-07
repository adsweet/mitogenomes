library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(nlme)
library(ape)
library(phytools)
library(rstatix)
library(geiger)
library(caper)

#Read in data
lice <- read.csv('at_comp_lice.csv')

mean.insect.at <- 76.0 #Average AT% of other insects

#Plot box plots
p2 <- ggplot(lice, aes(Structure, AT)) #all sites
p3 <- ggplot(lice, aes(Structure, X4Fold)) #four-fold degenerate sites
#sample(cbp2,3)
p3 + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_jitter(width=0.25, size=4, aes(color=Groupings))+
  scale_color_manual(values = cbp2[c(2,5,6)])+
  ylab("AT%")+
  theme_classic()+
  geom_hline(yintercept = mean.insect.at, linetype="dashed",color=cbp2[7])
  scale_shape_manual(values = c(18,19,19))+
  theme(aspect.ratio = 1)

#Melt matrix
at_m <- melt(lice, id=c("Genus","Groupings","Code","Higher_taxonomy","GenBank","Study","Structure","Species","Chromosomes"), variable.name = "Source")

#Plot AT% from to different codon positions
ggplot(at_m, aes(Source, value, fill=Structure))+
  geom_boxplot()+
  theme_classic()+
  ylab("AT%")+
  scale_fill_manual(values = cbp2)+
  geom_hline(yintercept = mean.insect.at, linetype="dashed",color=cbp2[7])
ggsave("Lice_mt_AT_seq_codon_structure.pdf")

#Pair-wise t-tests (with correction) for different subsets of mitogenome data
t.tests <- at_m %>% t_test(value ~ Source) %>%
  adjust_pvalue(method = "BH") %>% add_significance()

#Plot results from t-tests
myplot <- ggboxplot(
  at_m, x = "Source", y = "value",
  fill = "Source", palette = cbp2[c(2:8)], legend = "none",
  ggtheme = theme_pubr(border = TRUE)
)
t.tests <- t.tests %>% add_xy_position(x = "Source")
myplot + stat_pvalue_manual(t.tests, label = "p.adj.signif")

#Phylogenetically-informed tests

#read in tree
louse_tree <- read.tree("lice_mt_concat_pcgs_partitions.tre")
louse_tree <- root(louse_tree, "Likel")
louse_tree <- drop.tip(louse_tree, "Likel")

#Recode data
code <- lice$Code
group <- as.vector(lice$Groupings)
struct <- as.vector(lice$Structure)
at <- as.vector(lice$AT)
at_4 <- as.vector(lice$X4Fold)
at_s <- lice$AT_skew
gc_s <- lice$GC_skew
names(singv) <- scode
names(at) <- code
names(at_4) <- code
names(at_s) <- code
names(gc_s) <- code
names(group) <- code
names(struct) <- code

#AT%
shapiro.test(lice$AT)
ggqqplot(lice$AT)
t.test(AT~Structure, data=lice)

phylANOVA(louse_tree, struct, at)

#AT% 4-fold degenerate sites
shapiro.test(lice$X4Fold)
ggqqplot(lice$X4Fold)
t.test(X4Fold~Structure, data=lice)

phylANOVA(louse_tree, struct, at_4)

#Correlation between fragmented length and AT%

len.tab <- read.csv('at_comp_lice.csv', row.names = 2)
single <- filter(len.tab, Architecture == "single") #Only taxa with single mitogenomes
len.tab.frag <- filter(len.tab, Structure == "Fragmented") #Only taxa with fragmented mitogenomes
#trim tree; only taxa with fragmented mitogenomes
mt_tree_frag <- drop.tip(louse_tree, c("Atict","Brant","Btmac","Cacom","Dobre","Famar","Fulon","Htspi","Ibbid","Mgtat","Oscro","Risp","Tfbab","Amyr","Colp","Faqua","Het","Ibbid","TfbabPrver"))

#Check names
name.check(louse_tree_frag,len.tab.frag)
name.check(mt_tree_frag,len.tab.frag)
name.check(louse_tree, len.tab)

#Run models for all (len.tab) and fragmented-only (len.tab.frag) mitogenomes for total and four-fold degenerate sites
#Pagel's lambda
model1 <- gls(avg_length~AT, data=len.tab, correlation=corPagel(1,louse_tree, form = ~Species))
summary(model1)

model2 <- gls(avg_length~X4Fold, data=len.tab, correlation=corPagel(1,louse_tree, form = ~Species))
summary(model1)

model3 <- gls(avg_length~AT, data=len.tab.frag, correlation=corPagel(1,mt_tree_frag, form = ~Species))
summary(model1)

model4 <- gls(avg_length~X4Fold, data=len.tab.frag, correlation=corPagel(1,mt_tree_frag, form = ~Species))
summary(model1)

#Brownian motion
model5 <- gls(avg_length~AT, data=len.tab, correlation=corBrownian(1,louse_tree, form = ~Species))
summary(model2)

model6 <- gls(avg_length~X4Fold, data=len.tab, correlation=corBrownian(1,louse_tree, form = ~Species))
summary(model2)

model7 <- gls(avg_length~AT, data=len.tab.frag, correlation=corBrownian(1,mt_tree_frag, form = ~Species))
summary(model2)

model8 <- gls(avg_length~X4Fold, data=len.tab.frag, correlation=corBrownian(1,mt_tree_frag, form = ~Species))
summary(model2)

#Plot
ggscatter(len.tab, x="avg_length", y = "AT", add="reg.line", conf.int=T, cor.coef=T, cor.method='spearman',
          color="black",shape = 16, add.params = list(color = "black", fill = "lightgray", linetype=1),
          size=3, repel=T, xlab = "Sequence length", ylab="AT%")


