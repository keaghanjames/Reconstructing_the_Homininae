install.packages("vegan")
require(vegan)

setwd("~/Desktop/Analyses/Cultural_Traits/Culture_shannon/wVegan")
Pan <- read.csv("Pan.csv")
rownames(Pan) <- Pan$X
Pan[1] <- NULL
Pan <- t(Pan)
Pan <- diversity(Pan, index = "shannon")

Gorilla <- read.csv("Gorilla.csv")
rownames(Gorilla) <- Gorilla$X
Gorilla[1] <- NULL
Gorilla <- t(Gorilla)
Gorilla <- diversity(Gorilla, index = "shannon")

Pongo <- read.csv("Pongo.csv")
rownames(Pongo) <- Pongo$X
Pongo[1] <- NULL
Pongo <- t(Pongo)
Pongo <- diversity(Pongo, index = "shannon")

Shannon <- c(Pongo, Pan, Gorilla)
Shannon <- as.data.frame(Shannon)
write.csv(Shannon, file = "Shannon.csv")
