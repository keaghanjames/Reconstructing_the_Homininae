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


#create a list containing the ages of all the trees - age of the root
tree.ages <- lapply(alltrees, FUN = function(x) branching.times(x))
tree.ages <- lapply(tree.ages, FUN = function(x) x[1])
#extract the oldest tree - ie the tree with the deepest root
old.tree <- which.max(tree.ages)
old.tree <- alltrees[[old.tree]]
#extract the youngest tree
young.tree <- which.min(tree.ages)
young.tree <- alltrees[[young.tree]]

#create a list containing the gamma stat values for each of the trees
tree.gamma <- lapply(alltrees, FUN = function(x) gammaStat(x))

#get tree with highest gamma value
max.gamma.tree <- which.max(tree.gamma) #the shallowest tree
max.gamma.tree <- names(max.gamma.tree)
max.gamma.tree <- alltrees[[max.gamma.tree]]

#get tree with minimum gamma value
min.gamma.tree <- which.min(tree.gamma)#the deepest tree
min.gamma.tree <- names(min.gamma.tree)
min.gamma.tree <- alltrees[[min.gamma.tree]]



#save each of the selected trees
write.nexus(long.tree, file = "long.tree")
write.nexus(short.tree, file = "short.tree")
write.nexus(old.tree, file = "old.tree")
write.nexus(young.tree, file = "young.tree")
write.nexus(max.gamma.tree, file = "max.gamma.tree")
write.nexus(min.gamma.tree, file = "min.gamma.tree")
write.nexus(median.gamma.tree, file = "median.gamma.tree")

#the MCC tree was extracted from the posterior sample using TreeAnnotator v2.4.8 part of the BEAST software package                     
                     
#creating some summary figures
#first we convert the ages and alpha vectors from lists to numeric values
tree.ages.raw <- as.numeric(tree.ages)
tree.gamma.raw <- as.numeric(tree.gamma)
                     
#then check for normality and extract some summary statistics
shapiro.test(tree.ages.raw)
median(tree.ages.raw)
IQR(tree.ages.raw)
mean(tree.ages.raw)
sd(tree.ages.raw)

shapiro.test(tree.gamma.raw)
median(tree.gamma.raw)
IQR(tree.gamma.raw)
mean(tree.gamma.raw)
sd(tree.gamma.raw)

#create a histogram of these values
png(filename = "tree.ages.hist.png", width = 500, height = 500, pointsize = 20)
hist(tree.ages.raw, xlim = c(0,50), ylim = c(0,250), 
     col = "lightblue", xlab = "Tree Ages", main = "")
dev.off()

png(filename = "gamma.scores.hist.png", width = 500, height = 500, pointsize = 20)
hist(tree.gamma.raw, breaks = 12, ylim = c(0,250), main = "", xlim = c(2,4), 
     col = "lightblue", xlab = "Gamma Scores")
dev.off()

#Create a plot comparing gamma and lambda values                     
par(xpd = F)
png(filename = "age.v.gamma.png", height = 800, width = 800, pointsize = 20)
plot(tree.ages, tree.gamma, xlab = "Tree Age", ylab = "Gamma Value", pch = 20)
abline(lm(tree.gamma.raw~tree.ages.raw), col = "red")
dev.off()

#extract the branching times to a data fram
makeit <- function(x) {colbind(branching.times(x))}
split.times <- data.frame()
split.times <- lapply(trees, FUN = function(x) makeit(x))
rownames(split.times) <- c("Old_Tree", "Young_Tree", "Deep_Tree", "Shallow_Tree", "MCC_Tree")
