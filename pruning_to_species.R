
#load in the mcc tree
setwd("~/Desktop/Analyses/Phylogenies")
tree <- read.nexus("mcc.tree")

#create a vector containing the names of the tips we want to conserve
species <- c("Gorilla_beringei_graueri","Gorilla_gorilla_gorilla", 
"Homo_sapiens", "Pan_paniscus", "Pan_troglodytes_ellioti", "Pongo_abelii")

#prune the tree such that only single tip per species remains
tree <- drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])

#rename those tips to the species names
tips <- c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes", "Pongo_abelii")
tree$tip.label <- tips

#save the new tree
setwd("species_level")
write.nexus(tree, file = "species.level.tree")
