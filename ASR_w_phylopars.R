#the following script executes the procedure described in section 3.4.2.

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
#set the directory to wherever you want the experiment to run from - wherever your trait data is stored
setwd("~/Desktop/Analyses/Non_Cultural_Traits/body_size_male/")
output <- vector()#a vector into which we will place the reulst of each analysis

#we begin the loop - it repeates the analysis across of all the trees included in trees - theoretically this could be
#used to analyse many more trees however here we limit the analysis to the five selected
for (i in 1:length(trees)) {
  tryCatch({
#read in the trait data - this will need to change dependening on the reconstructed trait
  traits <- read.csv("body_size_male.csv")
homtree <- trees[[i]]
numbers <- vector()#a vector to which the results of interest will be placed

#having set up the data frame we now log-transform the trait values
traits[2:3] <- log(traits[2:3]) #do not do this when using the cultural diversity index, these are already log values

#next we check whether the trait actually shows any phylogenetic signal
#first calculate Pagel's lambda for the trait
lambda <- phylopars(traits, homtree, model = "lambda")
numbers <- c(numbers, lambda$model$lambda)#save the value of lambda
numbers <- c(numbers, lambda$logLik)#save the log-likelihood

#then we model the traits but assume a star phylogeny - this will act as our null
star <- phylopars(traits, homtree, model = "star")
numbers <- c(numbers, star$logLik) #again we save the log-likelihood

#now we can do a log-likelihood ratio test using the two models
chi_square <- as.double(2*(logLik(lambda) - logLik(star))) # 2*(logLik_alt - logLik_null)
                             
numbers <- c(numbers, chi_square)# save the chi_square values
degrees_freedom <- lambda$npars - star$npars # df = difference in # model parameters
numbers <- c(numbers, degrees_freedom) #save the degress of freedom
p_value <- pchisq(q = chi_square,df = degrees_freedom,lower.tail = FALSE) # save the p-value
numbers <- c(numbers, p_value)

#reconstrcut under various phylogenetic models
BM <- phylopars(trait_data = traits, homtree, model = "BM") #the Brownian motion model
#save(BM, file = "BM.Rdata")
OU <- phylopars(traits, homtree, model = "OU") #the Ornstein-Uhlenbeck model
#save(OU, file = "OU.Rdata")
EB <- phylopars(traits, homtree, model = "EB")# the early burst model
#save(EB, file = "EB.Rdata")

#combine outputs into single list
models <- list(BM, OU, EB)
#get AIC for each of the models
best_model <- c(AIC(BM), AIC(OU), AIC(EB))

#now identify the model with the lowest AIC
x <- which.min(best_model)
numbers <- c(numbers, best_model[[x]]) #record which model it is in numbers
best_model <- models[[x]] #select that model
numbers <- c(numbers, best_model$model$model, best_model$mu, best_model$npars) #save its most important parameters 
save(best_model, file = paste(treenames[i]))#save the data from this model incase we need it later

#now we extract the ancestral state estimates and the variance
node_states <- best_model$anc_recon # Data with ancestral state estimates
node_variance <- best_model$anc_var # Variances for each estimate

#and compute 95%CI by using the square root
lowerCI <- best_model$anc_recon - sqrt(best_model$anc_var)*1.96 # Lower 95% CI
upperCI <- best_model$anc_recon + sqrt(best_model$anc_var)*1.96 # Upper 95% CI

#if you log transformed the traits you'll want to reverse that before you make the figures
node_states <- exp(node_states) #convert back to meaningful values
node_variance <- exp(node_variance) #convert back to meaningful value
lowerCI <- exp(lowerCI)#convert back to meaningful value
upperCI <- exp(upperCI)#convert back to meaningful value

#to make some figures we first need to organise the tip and node values into matrices
node_states <- as.matrix(node_states)[,1]
tip_states <- node_states[1:length(homtree$tip.label)]

tips.ancestors <- node_states
names(tips.ancestors) <- NULL
lower <- as.vector(lowerCI)
upper<- as.vector(upperCI)
tips.ancestors <- as.vector(rbind(tips.ancestors, lower, upper))

#and finally we place our ancestral state estimates into numbers
numbers <- c(numbers, tips.ancestors)
#now we make a phenogram
fancyTree(homtree, type = "phenogram95", x = tip_states)
phenogram(homtree, node_states, add = T, ftype = "off")

#and we make  phylogeny with the coloured nodes
tips_states_cut <- cut(node_states[1:length(homtree$tip.label)], breaks = 5, labels = F)
node_states_cut <- cut(node_states[length(homtree$tip.label)+1:length(node_states)], breaks = 5, labels = F)
plot(homtree, type="p", use.edge.length = TRUE, label.offset=1,cex=0.6) 
tiplabels(pch=21,bg=CP[tips_states_cut],col="black",cex=2,adj=0.505) 
nodelabels(pch=21,bg=CP[node_states_cut],col="black",cex=2,adj=0.505) 
op<-par(xpd=TRUE)
legend(0,0, legend=levels(cut(node_states, breaks=5)),
       col=CP,pch=20,bty="n",cex=0.7,pt.cex=1.5,ncol=5)
  }, error=function(e){}) #this ensures that the loop continues even if it encounters an error
  output <- rbind(output, numbers) #place numbers in output
}
rownames(output) <- c("mcc.tree", "old.tree", "young.tree", "deep.tree", "shallow.tree") #name the rows by tree
thestats <- c("lambda", "lambda_loglik", "star_loglik", "chi_square", "df", "p_value", "AIC", "Model", "mu", "npars",
                      "Gorilla_beringei_graueri", "lower95", "upper95", "Gorilla_beringei_beringei", "lower95", "upper95", 
                      "Gorilla_gorilla_diehli", "lower95", "upper95", "Gorilla_gorilla_gorilla",  "lower95", "upper95",
                      "Homo_sapiens", "lower95", "upper95", "Pan_paniscus", "lower95", "upper95", "Pan_troglodytes_ellioti", 
                      "lower95", "upper95", "Pan_troglodytes_verus", "lower95", "upper95","Pan_troglodytes_schweinfurthii", 
                      "lower95", "upper95","Pan_troglodytes_troglodytes", "lower95", "upper95",	"Pongo_abelii", "lower95", 
                      "upper95", "root", "lower95", "upper95", "LCAGHP", "lower95", "upper95", "LCAG", "lower95", "upper95", 
                      "LCAGb", "lower95", "upper95", "LCAGg", "lower95", "upper95", "LCAPH", "lower95", "upper95", "LCAP", 
                      "lower95", "upper95", "LCAPt", "lower95", "upper95", "LCAPtve", "lower95", "upper95", "LCAPtts", "lower95", 
                      "upper95") #name the columns - there's a lot of them. 
colnames(output) <- thestats #give the columns their names
write.csv(output, file = "output.csv") #save the output of the ancestral state reconstruction
