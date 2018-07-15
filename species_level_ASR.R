#load in requisite packages
require(RColorBrewer)
require(ape)
require(phytools)
require(Rphylopars)

numbers <- vector()

set.seed(888)

setwd("~/Desktop/Analyses/Cultural_Traits/Culture_count/no_homo")
  the_traits <- read.csv("Culture_count_binary.csv")
setwd("~/Desktop/Analyses/Phylogenies/species_level/")
homtree <- read.nexus("species.level.tree")

setwd("culture_count_no_homo/averages")

lambda_values <- vector()
p_values <- vector()
X2 <- vector()
df <- vector()

trog <- subset(the_traits, species == "Pan_troglodytes_verus" | 
                 species == "Pan_troglodytes_troglodytes" | 
                 species == "Pan_troglodytes_schweinfurthii" |
                 species == "Pan_troglodytes_elloti")
#trog$species <- "Pan_troglodytes"
trog <- na.omit(trog)
trog <- trog[[2]]

paniscus <- subset(the_traits, species == "Pan_paniscus")
paniscus <- na.omit(paniscus)
paniscus <- paniscus[[2]]

sapiens <- subset(the_traits, species == "Homo_sapiens")
sapiens <- na.omit(sapiens)
sapiens <- sapiens[[2]]
sapiens <- NA #just for the cultural traits

gorilla <- subset(the_traits, species == "Gorilla_gorilla_gorilla" | 
                 species == "Gorilla_gorilla_diehli")
gorilla <- na.omit(gorilla)
gorilla <- gorilla[[2]]
#gorilla <- NA #just for age at first repo

beringei <- subset(the_traits, species == "Gorilla_beringei_beringei" | 
                    species == "Gorilla_beringei_graueri")
beringei <- na.omit(beringei)
beringei <- beringei[[2]]
#beringei <- NA #just for com. size

abelii <- subset(the_traits, species == "Pongo_abelii" | species == "Pongo_pygmaeus")
abelii <- na.omit(abelii)
abelii <- abelii[[2]]

species <- c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapiens", 
             "Pan_paniscus", "Pan_troglodytes", "Pongo_abelii")

for (i in 1:100) {
traits <- c(sample(beringei, 1), sample(gorilla, 1), sample(sapiens, 1), 
           sample(paniscus, 1), sample(trog, 1), sample(abelii, 1))

#traits <- c(median(beringei), median(gorilla), median(sapiens), 
            #median(paniscus), median(trog), median(abelii)) #for taking cultural averages

traits <- data.frame(species, traits)
rownames(traits) <- traits$species
traits[2] <- log(traits[2]) #turn off for Shannon index

lambda <- phylopars(traits, homtree, model = "lambda")
lambda_values <- c(lambda_values, lambda$model$lambda)
star <- phylopars(traits, homtree, model = "star")
chi_square <- as.double(2*(logLik(lambda) - # 2*(logLik_alt - logLik_null)
                             logLik(star))) 
X2 <- c(X2, chi_square)
degrees_freedom <- lambda$npars - star$npars # df = difference in # model parameters
df <- c(df, degrees_freedom)
p_value <- pchisq(q = chi_square,df = degrees_freedom,lower.tail = FALSE) # p-value
p_values <- c(p_values, p_value)
print(i)
}


png(filename = "lambda_values.png", width = 1000, height = 1000, pointsize = 20)
hist(lambda_values, ylim = c(0,100), xlim = c(0,1), 
     main = "", xlab = "Lambda values", col = "lightblue")
dev.off()

png(filename = "p_values.png", width = 1000, height = 1000, pointsize = 20)
hist(p_values, ylim = c(0,60), xlim = c(0,1), main = "", xlab = "p-values",
     col = "lightblue")
dev.off()
the_goods <- c(median(lambda_values), quantile(lambda_values)[[2]], 
               quantile(lambda_values)[[4]], length(which(p_values < 0.05)))
numbers <- rbind(numbers, the_goods)
results <- rbind(lambda_values, quantile(lambda)[[2]], quantile(lambda)[[2]], X2, df, p_values)
results <- t(results)
write.csv(results, file = "results.csv")

