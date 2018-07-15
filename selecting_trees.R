#load necessary package
require(ape)
require(phytools)

setwd("~/Desktop/Analyses/Phylogenies")
#create tip for Gorilla beringei beringei and identify sister
tip <- "Gorilla_beringei_beringei"
sister <- "Gorilla_beringei_graueri"
#create a function so that we can quickly bind the G.b.beringei tip to each tree at half way along the bl of its sister species G.b.graueri
gorilla.bind <- function(x) {bind.tip(x,tip,where=which(x$tip.label==sister), 
                                      position=0.5*x$edge.length[which(x$edge[,2]==which(x$tip.label==sister))])
}

#read in the trees
alltrees <- read.nexus("trimmed_samples_tree.trees")

#apply gorilla.bind to alltree
alltrees <- lapply(alltrees, FUN = function(x) gorilla.bind(x))

#drop the first 20% of trees as burnin
alltrees <- alltrees[201:1001] #ensure odd number so that medians are actually represented in posterior
write.nexus(alltrees, file = "alltrees_burnt.trees")

#create a list containing the sum edge lengths of each sampled tree
tree.lengths <- lapply(alltrees, FUN = function(x) sum(x$edge.length))
#extract the longest tree
long.tree <- which.max(tree.lengths)
#long.tree <- names(long.tree)
long.tree <- alltrees[[long.tree]]
#long.tree <- gorilla.bind(long.tree)


#now do the same but for the shortest tree
short.tree <- which.min(tree.lengths)
#short.tree <- names(short.tree)
short.tree <- alltrees[[short.tree]]
#short.tree <- gorilla.bind(short.tree)

#create a list containing the ages of all the trees - age of the root
tree.ages <- lapply(alltrees, FUN = function(x) branching.times(x))
tree.ages <- lapply(tree.ages, FUN = function(x) x[1])
#extract the oldest tree - ie the tree with the deepest root
old.tree <- which.max(tree.ages)
#old.tree <- names(old.tree)
old.tree <- alltrees[[old.tree]]
#old.tree <- gorilla.bind(old.tree)

#extract the youngest tree
young.tree <- which.min(tree.ages)
#young.tree <- names(young.tree)
young.tree <- alltrees[[young.tree]]
#young.tree <- gorilla.bind(young.tree)

#create a list containing the gamma stat values for each of the trees
tree.gamma <- lapply(alltrees, FUN = function(x) gammaStat(x))

#get tree with highest gamma value
max.gamma.tree <- which.max(tree.gamma) #the shallowest tree
max.gamma.tree <- names(max.gamma.tree)
max.gamma.tree <- alltrees[[max.gamma.tree]]
#max.gamma.tree <- gorilla.bind(max.gamma.tree)

#get tree with minimum gamma value
min.gamma.tree <- which.min(tree.gamma)#the deepest tree
min.gamma.tree <- names(min.gamma.tree)
min.gamma.tree <- alltrees[[min.gamma.tree]]
#min.gamma.tree <- gorilla.bind(min.gamma.tree)

#get tree with median gamma value
median.gamma.tree <- as.numeric(tree.gamma)
median.gamma.tree <- median(median.gamma.tree)
median.gamma.tree <- which(tree.gamma == median.gamma.tree)#this will only work if there is an odd number of trees in the distribution, otherwise the median gamma value will be the avreage of the two on either side of the centre and thus will not exist in the original distribution
median.gamma.tree <- alltrees[[median.gamma.tree]]
#median.gamma.tree <- gorilla.bind(median.gamma.tree)

#get the Maximum Clade Credibility Tree
#mcc.tree <- maxCladeCred(alltrees) #why is it different to TreeAnnotator?
#mcc.tree <- gorilla.bind(mcc.tree)
#no just extract with tree annotator
#add gorilla to mcc tree#
#mcc.tree <- read.nexus("mcc.tree")
#mcc.tree <- gorilla.bind(mcc.tree)

write.nexus(long.tree, file = "long.tree")
write.nexus(short.tree, file = "short.tree")
write.nexus(old.tree, file = "old.tree")
write.nexus(young.tree, file = "young.tree")
#write.nexus(mcc.tree, file = "mcc.tree")
write.nexus(max.gamma.tree, file = "max.gamma.tree")
write.nexus(min.gamma.tree, file = "min.gamma.tree")
write.nexus(median.gamma.tree, file = "median.gamma.tree")

#creating some summary figures
#first we convert the ages and alpha vectors from lists to numeric
tree.ages.raw <- as.numeric(tree.ages)
shapiro.test(tree.ages.raw)
median(tree.ages.raw)
IQR(tree.ages.raw)
mean(tree.ages.raw)
sd(tree.ages.raw)

tree.gamma.raw <- as.numeric(tree.gamma)
shapiro.test(tree.gamma.raw)
median(tree.gamma.raw)
IQR(tree.gamma.raw)
mean(tree.gamma.raw)
sd(tree.gamma.raw)

png(filename = "tree.ages.png", width = 1000, height = 800, pointsize = 16)
age.density <- density(tree.ages.raw)
plot(age.density, xlab = "Tree ages", main = "Distribution of tree ages")
polygon(age.density, col = rgb(0,0,1, .5))
abline(v = median(tree.ages.raw), lty = 2, lwd = 2)
dev.off()

png(filename = "tree.ages.hist.png", width = 500, height = 500, pointsize = 20)
hist(tree.ages.raw, xlim = c(0,50), ylim = c(0,250), 
     col = "lightblue", xlab = "Tree Ages", main = "")
dev.off()

png(filename = "gamma.density.png", width = 500, height = 500, pointsize = 16)
gamma.density <- density(tree.gamma.raw)
plot(gamma.density, main = "Distribution of tree gamma values", xlab = "Gamma value",
     xlim = c(2, 3.8), ylim = c(0, 3.5))
polygon(gamma.density, col = rgb(0,0,1, .5))
abline(v = median(tree.gamma.raw), lty = 2, lwd = 2)
dev.off()

png(filename = "gamma.scores.hist.png", width = 500, height = 500, pointsize = 20)
hist(tree.gamma.raw, breaks = 12, ylim = c(0,250), main = "", xlim = c(2,4), 
     col = "lightblue", xlab = "Gamma Scores")
dev.off()

par(xpd = F)
png(filename = "age.v.gamma.png", height = 800, width = 800, pointsize = 20)
plot(tree.ages, tree.gamma, xlab = "Tree Age", ylab = "Gamma Value", pch = 20)
abline(lm(tree.gamma.raw~tree.ages.raw), col = "red")
dev.off()

#get branching times
makeit <- function(x) {colbind(branching.times(x))}
split.times <- data.frame()
split.times <- lapply(trees, FUN = function(x) makeit(x))
rownames(split.times) <- c("Old_Tree", "Young_Tree", "Deep_Tree", "Shallow_Tree", "MCC_Tree")
