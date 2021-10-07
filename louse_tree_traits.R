library(ape)
library(phyloch)
library(strap)
library(phylotate)
library(MCMCtreeR)
library(phytools)
library(geiger)
library(ggtree)
library(caper)

#Read in Johnson et al. (2018) tree
tr <- read.nexus('Dated200LabelledTree.tre')
tr <- root(tr, "Stjap") #root
#remove outgroup and extra taxa
tr_trim <- drop.tip(tr, c("Shlar","Pnbad","Ptgot","Ptpub","Nepac","Hbarb","Ltspi","Hieur","Aamic","Prflu"))
tr_trim$tip.label[tr_trim$tip.label=="Geaur"] <- "Bobov" #change name of a tip
tr_trim$root.time <- 119.5644 #root age
tr_trim$edge.length <- (tr_trim$edge.length)*100 #convert ages to years
#Plot tree with geologic scale
geoscalePhylo(tree=ladderize(tr_trim,right=FALSE), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-15,120), lwd=3, width=2,quat.rm = T)

write.tree(tr_trim,"biol_lett_louse_tree_trim.tre")

#read in table with mitogenome structure
trait <- read.csv("louse_mt_traits.csv",header=T)
arch <- setNames(trait[,2],trait[,1])
arch <- setNames(as.factor(trait[,5]),trait[,1])
#arch <- as.matrix(trait)[,1]

cols<-setNames(palette()[1:length(unique(arch))],levels(arch)) #set colors
#plot tip labels
tiplabels(pie=to.matrix(arch[tr_trim$tip.label],levels(arch)),piecol=cols,cex=0.4)
add.simmap.legend(colors=cols,prompt=F,x=8,
                  y=8,fsize=0.5,shape='circle',vertical = F)

#Trait reconstruction model comparisons
fitER <- ace(arch,tr_trim,model='ER',type='discrete')
fitARD <- ace(arch,tr_trim,model='ARD',type='discrete')
fitSYM <- ace(arch,tr_trim,model='SYM',type='discrete')
fitUSR <- ace(arch, tr_trim, model=matrix(c(0,1,0,0),2),type='discrete') #custom model, not allowing transition from fragmented to single
model_geiger_ER <- fitDiscrete(tr_trim,arch,type='discrete',model="ER")
model_geiger_ARD <- fitDiscrete(tr_trim,arch,type='discrete',model="ARD")
model_geiger_SYM <- fitDiscrete(tr_trim,arch,type='discrete',model="SYM")
model_geiger_USR <- fitDiscrete(tr_trim,arch,type='discrete',model=matrix(c(0,1,0,0)))

ERvSYM <- abs(model_geiger_ER$opt$aicc - model_geiger_SYM$opt$aicc)
ERvARD <- abs(model_geiger_ER$opt$aicc - model_geiger_ARD$opt$aicc)
ERvUSR <- abs(model_geiger_ER$opt$aicc - model_geiger_USR$opt$aicc)
SYMvARD <- abs(model_geiger_SYM$opt$aicc - model_geiger_ARD$opt$aicc)

compare AIC and AICc for different models
comp <- data.frame(model = c("ER","SYM","ARD","USR"),AIC = c(model_geiger_ER$opt$aic,model_geiger_SYM$opt$aic,model_geiger_ARD$opt$aic,model_geiger_USR$opt$aic),
                   AICc = c(model_geiger_ER$opt$aicc,model_geiger_SYM$opt$aicc,model_geiger_ARD$opt$aicc,model_geiger_USR$opt$aicc))

#Other model comparison
1-pchisq(2*abs(fitER$loglik - fitARD$loglik), 1)
1-pchisq(2*abs(fitER$loglik - fitUSR$loglik), 1)
1-pchisq(2*abs(fitARD$loglik - fitUSR$loglik), 1)

#Plot pie charts on nodes
round(fitER$lik.anc,3)
round(fitARD$lik.anc,3)
round(fitUSR$lik.anc,3)

plotTree(tr_trim)

nodelabels(node=1:tr_trim$Nnode+Ntip(tr_trim),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(arch[tr_trim$tip.label],levels(as.factor(arch))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=F,x=8,
                  y=8,fsize=0.5,shape='circle',vertical = F)

tr_trim <- ladderize(tr_trim,right = F)
is_tip <- tr_trim$edge[,2] <= length(tr_trim$tip.label)
ordered_tips <- tr_trim$edge[is_tip, 2]
o <- tr_trim$tip.label[ordered_tips]

#Plot barplot of # of chromosomes
barplot(co,col=sapply(co,function(co) if(co==1) "red" else "blue"),xlim=c(0,max(co)),space=0,xlab="Chromosomes",horiz = T)

#Tests for phylogenetic signal
datasig_full <- read.csv("louse_mt_traits.csv")
sig <- phylo.d(datasig,tr,names.col = X,binvar=Architecture)
summary(sig)
full_data <- comparative.data(tr2,datasig_full, X, na.omit = F)
sig_full <- phylo.d(full_data, binvar = Architecture)
summary(sig_full)


#Stochastic character mapping
mtrees<-make.simmap(tr_trim,arch,model=matrix(c(0,1,0,0),2),nsim=1000)
pd<-summary(mtrees,plot=F)
pd
plot(pd,fsize=0.6,ftype="i")
add.simmap.legend(colors=cols,prompt=F,x=8,
                  y=8,fsize=0.5,shape='circle',vertical = F)

plot(mtrees[[10]],cols,fsize=0.8,ftype="i")
nodelabels(pie=pd$ace,piecol=cols,cex=0.5)

