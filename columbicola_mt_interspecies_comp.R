library(ape)
library(ggplot2)
library(RColorBrewer)
library(adegenet)
library(dplyr)
library(reshape2)
library(phangorn)

#Read in gene tables
bin <- read.table('columbicola_mt_comp_gene_boundaries.txt', sep='\t',header=T, row.names = 1) #With outgroup taxa
bin <- read.table('copas2_mt_comp_gene_boundaries.txt', sep='\t',header=T, row.names = 1, na.strings = "M", check.names = F) #Only C. passerinae 2
copas2.gene.dist <- dist(bin, method = "binary") #Binary distance matrix from gene boundary data

#UPGMA distance tree from gene boundary matrix
copas2.gene.tree <- upgma(copas2.gene.dist)
copas2.gene.tree.root <- root(copas2.gene.tree, "Cocol")
plot(copas2.gene.tree)

#Neighbor-joining distance tree from genetic data
copas2.mt.seq <- read.FASTA("copas2_mt_alignment.fasta")
copas2.mt.dist <- dist.dna(copas2.mt.seq, pairwise.deletion = T, model = "JC69")
copas2.mt.tree <- nj(copas2.mt.dist)
copas2.mt.tree.root <- root(copas2.mt.tree, "Cocol")
plot(copas2.mt.tree.root)

gene.d <- as.matrix(copas2.gene.dist)
seq.d <- as.matrix(copas2.mt.dist)

#Order genetic boundary matrix according to the gene genetic distance matrix
gene.d <- gene.d[rownames(seq.d), rownames(seq.d)]

#Compare gene boundary and genetic distance matrices
mantel.test(gene.d, seq.d, nperm = 999, graph = T)
mantel.randtest(as.dist(gene.d), as.dist(seq.d), nrepet = 999)

#Plot gene boundaries
bound <- sapply(bin, function(x) sum(x, na.rm = T)/length(x[!is.na(x)]))
bound <-data.frame(lapply(bound, type.convert), stringsAsFactors=F, check.names = F)
bound.m <- data.frame(Boundary=colnames(bound), Proportion=as.vector(t(bound)))
write.csv(bound.m, "copas2_mt_gene_boundaries_prop.csv")

prop <- read.csv("copas2_mt_gene_boundaries_prop.csv", header=T)
prop$Boundary <- factor(prop$Boundary, levels = unique(prop$Boundary))

cCount <- length(unique(prop$Chrom))
getP <- colorRampPalette(brewer.pal(8, "Set2"), bias=0.9)

ggplot(prop, aes(Boundary, Proportion, fill=Chrom))+
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
  scale_fill_manual(values = getP(cCount))

ggsave("copas2_boundaries_classic.pdf")

chrom.count <- read.csv("copas2_mt_chrome_arrangment_counts.csv")

ggplot(chrom.count, aes(Chrom, Count))+
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
  scale_fill_manual(values = getP(cCount))
ggsave("copas2_chromosome_boundary_count_classic.pdf")
