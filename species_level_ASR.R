#load in required packages
require(RColorBrewer)
require(ape)
require(phytools)
require(Rphylopars)

#create an empty vector
numbers <- vector()

#because this procdure draws values at random we set the seed to 888 for reproducibility
set.seed(888)

#set working director the location of the data - below is an example
setwd("~/Desktop/Analyses/Cultural_Traits/Culture_count/no_homo")
#read in the traits  
the_traits <- read.csv("Culture_count_binary.csv")

#set working direction to location of the phlogenies
setwd("~/Desktop/Analyses/Phylogenies/species_level/")
#read in the species level tree
homtree <- read.nexus("species.level.tree")

#now we enter a new working direction
setwd("culture_count_no_homo/averages")

#create some emptyr vectors to staw the values of interest from each run
lambda_values <- vector()
p_values <- vector()
X2 <- vector()
df <- vector()

#now we create species level subsets for our trait data
#if data is missing for an entire species (typically H. sapiens or one of the gorilla species) it may be necessary to assign
#a value of NA to these manually - for example sapiens <- NA - otherwise the code will not run

trog <- subset(the_traits, species == "Pan_troglodytes_verus" | 
                 species == "Pan_troglodytes_troglodytes" | 
                 species == "Pan_troglodytes_schweinfurthii" |
                 species == "Pan_troglodytes_elloti") #subset by subspecies names, uniting their values in one vectore
trog <- na.omit(trog) #delete any NAs
trog <- trog[[2]]# ensure that the vector is a single string of values

#process is repeated for each species
paniscus <- subset(the_traits, species == "Pan_paniscus")
paniscus <- na.omit(paniscus)
paniscus <- paniscus[[2]]

sapiens <- subset(the_traits, species == "Homo_sapiens")
sapiens <- na.omit(sapiens)
sapiens <- sapiens[[2]]

gorilla <- subset(the_traits, species == "Gorilla_gorilla_gorilla" | 
                 species == "Gorilla_gorilla_diehli")
gorilla <- na.omit(gorilla)
gorilla <- gorilla[[2]]

beringei <- subset(the_traits, species == "Gorilla_beringei_beringei" | 
                    species == "Gorilla_beringei_graueri")
beringei <- na.omit(beringei)
beringei <- beringei[[2]]

abelii <- subset(the_traits, species == "Pongo_abelii" | species == "Pongo_pygmaeus")
abelii <- na.omit(abelii)
abelii <- abelii[[2]]

#creat a vector containing the species names
species <- c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapiens", 
             "Pan_paniscus", "Pan_troglodytes", "Pongo_abelii")

#now the loop begins
for (i in 1:100) {
traits <- c(sample(beringei, 1), sample(gorilla, 1), sample(sapiens, 1), 
           sample(paniscus, 1), sample(trog, 1), sample(abelii, 1)) #for each species a value is sampled at randomly and placed in traits

#traits <- c(median(beringei), median(gorilla), median(sapiens), 
            #median(paniscus), median(trog), median(abelii)) #for taking cultural averages rather than sampling at random

traits <- data.frame(species, traits) #we then combine the list of species and the randomly sampled traits into a data frame
rownames(traits) <- traits$species #this is necessary step for phylopars otherwise it will crash to desktop
traits[2] <- log(traits[2]) #log the trait values as per section 3.4.2 - turn off for Shannon index as the value is logged

#now we model lambda
lambda <- phylopars(traits, homtree, model = "lambda") #use the phylopars command to fit lambda to the tree and trait data
lambda_values <- c(lambda_values, lambda$model$lambda) # add the lambda score to the vector 
star <- phylopars(traits, homtree, model = "star") #fit the star phlogeny
chi_square <- as.double(2*(logLik(lambda) - logLik(star))) # perform the log-likelihood test i.e. 2*(logLik_alt - logLik_null)
                              
X2 <- c(X2, chi_square) #extract the chi square 
degrees_freedom <- lambda$npars - star$npars # df = difference in # model parameters
df <- c(df, degrees_freedom) #extract these too
p_value <- pchisq(q = chi_square,df = degrees_freedom,lower.tail = FALSE) # p-value
p_values <- c(p_values, p_value) #extract the p-value
print(i) #allows us to track the progress of the loop
}

#create a histogram containing the lambda values
png(filename = "lambda_values.png", width = 1000, height = 1000, pointsize = 20)
hist(lambda_values, ylim = c(0,100), xlim = c(0,1), 
     main = "", xlab = "Lambda values", col = "lightblue")
dev.off()

#create a historgram containing the p-values
png(filename = "p_values.png", width = 1000, height = 1000, pointsize = 20)
hist(p_values, ylim = c(0,60), xlim = c(0,1), main = "", xlab = "p-values",
     col = "lightblue")
dev.off()

#get the summary statistics included in Table 4.2
the_goods <- c(median(lambda_values), quantile(lambda_values)[[2]], 
               quantile(lambda_values)[[4]], length(which(p_values < 0.05)))
numbers <- rbind(numbers, the_goods)
results <- rbind(lambda_values, quantile(lambda)[[2]], quantile(lambda)[[2]], X2, df, p_values)
results <- t(results)
write.csv(results, file = "results.csv")

