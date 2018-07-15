install.packages("vegan")
require(vegan)

#set workding directory to location of cultural frequency data
setwd("~/Desktop/Analyses/Cultural_Traits/Culture_shannon/wVegan")

#because we wanted the the number of possible traits to be at the level genera but not across all the hominids I set up each 
#genera's data in it's own csv.
Pan <- read.csv("Pan.csv")# read in the files
rownames(Pan) <- Pan$X #give the data frame rownames
Pan[1] <- NULL #remove the taxa names, they are now the rownames
Pan <- t(Pan) #transpose the dataframe
Pan <- diversity(Pan, index = "shannon") #calculate shannon's index

#repeate for the other two genera
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

#now we want a workable output
Shannon <- c(Pongo, Pan, Gorilla) #combine the genera
Shannon <- as.data.frame(Shannon)#turn into data frame
write.csv(Shannon, file = "Shannon.csv")#save as csv
