#load in requisite packages
require(RColorBrewer)
require(ape)
require(phytools)
require(Rphylopars)
#establish some global colours for future figures
CP <- brewer.pal(n = 5, name = "RdBu")

#read in the trees
setwd("~/Desktop/Analyses/Phylogenies/")
mcc.tree <- read.nexus("mcc.tree")
old.tree <- read.nexus("old.tree")
young.tree <- read.nexus("young.tree")
deep.tree <- read.nexus("min.gamma.tree")
shallow.tree <- read.nexus("max.gamma.tree")
trees <- c(mcc.tree, old.tree, young.tree, deep.tree, shallow.tree)
treenames <- c("mcc.tree.Rdata", "old.tree.Rdata", "young.tree.Rdata", "deep.tree.Rdata", "shallow.tree.Rdata")
#homtree <- mcc.tree
#set teh director to were you want the experimen to run
setwd("~/Desktop/Analyses/Non_Cultural_Traits/body_size_integrated/")
output <- vector()#a vector into which we will place the reulst of each analysis

for (i in 1:length(trees)) {
  tryCatch({
#read in the trait data - need to change according to the csv
  traits <- read.csv("body_size_integrated.csv")
  traits$X <- NULL #just for the integrated run
homtree <- trees[[i]]
numbers <- vector()#a vector to which the results of interest will be placed
#rearrange the rows so they are in the same order and have exactly the same names as the tips
#rownames(traits) <- traits$species
#traits <- traits[homtree$tip.label,]
#rownames(traits) <- homtree$tip.label #rename them all because a few get lost for not apparent reason
#traits$species <- homtree$tip.label

#having set up the data frame we now log-transform the trait values
traits[2:3] <- log(traits[2:3])

#next we check whether the trait actually shows any phylogenetic signal
#first calculate Pagel's lambda for the trait
lambda <- phylopars(traits, homtree, model = "lambda")
numbers <- c(numbers, lambda$model$lambda)
numbers <- c(numbers, lambda$logLik)
#SAVE LAMBDA
#then we model the traits but assume a star phylogeny - this acts as a null
star <- phylopars(traits, homtree, model = "star")
numbers <- c(numbers, star$logLik)
#SAVE STAR
#now we can do a log-likelihood to test the two models
chi_square <- as.double(2*(logLik(lambda) - # 2*(logLik_alt - logLik_null)
                             logLik(star))) 
numbers <- c(numbers, chi_square)
degrees_freedom <- lambda$npars - star$npars # df = difference in # model parameters
numbers <- c(numbers, degrees_freedom)
p_value <- pchisq(q = chi_square,df = degrees_freedom,lower.tail = FALSE) # p-value
numbers <- c(numbers, p_value)
print(p_value)

#reconstrcut under various phylogenetic models
BM <- phylopars(trait_data = traits, homtree, model = "BM")
#save(BM, file = "BM.Rdata")
mvOU <- phylopars(traits, homtree, model = "mvOU")
#save(mvOU, file = "mvOU.Rdata")
OU <- phylopars(traits, homtree, model = "OU")
#save(OU, file = "OU.Rdata")
EB <- phylopars(traits, homtree, model = "EB")
#save(EB, file = "EB.Rdata")

#combine outputs into single list
models <- list(BM, mvOU, OU, EB)
#get AIC for each of the models
best_model <- c(AIC(BM), AIC(mvOU), AIC(OU), AIC(EB))

#now identify the model with the lowest AIC
x <- which.min(best_model)
numbers <- c(numbers, best_model[[x]])
best_model <- models[[x]]
numbers <- c(numbers, best_model$model$model, best_model$mu, best_model$npars)
#print(best_model)
save(best_model, file = paste(treenames[i]))
#now we view the imputed missing data and the imputation variance 
node_states <- best_model$anc_recon # Data with imputed species means
node_variance <- best_model$anc_var # Variances for each estimate

#and compute 95%CI by using the square root
lowerCI <- best_model$anc_recon - sqrt(best_model$anc_var)*1.96 # Lower 95% CI
upperCI <- best_model$anc_recon + sqrt(best_model$anc_var)*1.96 # Upper 95% CI

#if you log transformed the traits you'll want to reverse that before you make the figures
node_states <- exp(node_states) #convert back to meaningful values
node_variance <- exp(node_variance) #convert back to meaningful value
lowerCI <- exp(lowerCI)
upperCI <- exp(upperCI)

#to make some figues we first need to organise the tip and node values into matrices
node_states <- as.matrix(node_states)[,1]
tip_states <- node_states[1:length(homtree$tip.label)]

tips.ancestors <- node_states
names(tips.ancestors) <- NULL
lower <- as.vector(lowerCI)
upper<- as.vector(upperCI)
tips.ancestors <- as.vector(rbind(tips.ancestors, lower, upper))
numbers <- c(numbers, tips.ancestors)
#now we make a phenogram
fancyTree(homtree, type = "phenogram95", x = tip_states)
phenogram(homtree, node_states, add = T, ftype = "off")

#and we make that phylogeny with the coloured nodes
tips_states_cut <- cut(node_states[1:length(homtree$tip.label)], breaks = 5, labels = F)
node_states_cut <- cut(node_states[length(homtree$tip.label)+1:length(node_states)], breaks = 5, labels = F)
plot(homtree, type="p", use.edge.length = TRUE, label.offset=1,cex=0.6) 
tiplabels(pch=21,bg=CP[tips_states_cut],col="black",cex=2,adj=0.505) 
nodelabels(pch=21,bg=CP[node_states_cut],col="black",cex=2,adj=0.505) 
op<-par(xpd=TRUE)
legend(0,0, legend=levels(cut(node_states, breaks=5)),
       col=CP,pch=20,bty="n",cex=0.7,pt.cex=1.5,ncol=5)
  }, error=function(e){})
  output <- rbind(output, numbers)
}
rownames(output) <- c("mcc.tree", "old.tree", "young.tree", "deep.tree", "shallow.tree")
thestats <- c("lambda", "lambda_loglik", "star_loglik", "chi_square", "df", "p_value", "AIC", "Model", "mu", "npars",
                      "Gorilla_beringei_graueri", "lower95", "upper95", "Gorilla_beringei_beringei", "lower95", "upper95", 
                      "Gorilla_gorilla_diehli", "lower95", "upper95", "Gorilla_gorilla_gorilla",  "lower95", "upper95",
                      "Homo_sapiens", "lower95", "upper95", "Pan_paniscus", "lower95", "upper95", "Pan_troglodytes_ellioti", 
                      "lower95", "upper95", "Pan_troglodytes_verus", "lower95", "upper95","Pan_troglodytes_schweinfurthii", 
                      "lower95", "upper95","Pan_troglodytes_troglodytes", "lower95", "upper95",	"Pongo_abelii", "lower95", 
                      "upper95", "root", "lower95", "upper95", "LCAGHP", "lower95", "upper95", "LCAG", "lower95", "upper95", 
                      "LCAGb", "lower95", "upper95", "LCAGg", "lower95", "upper95", "LCAPH", "lower95", "upper95", "LCAP", 
                      "lower95", "upper95", "LCAPt", "lower95", "upper95", "LCAPtve", "lower95", "upper95", "LCAPtts", "lower95", 
                      "upper95")
colnames(output) <- thestats
write.csv(output, file = "output.csv")

