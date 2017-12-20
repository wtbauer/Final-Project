
setwd("~/Desktop/SFSU files/Bio 738/1Final Project")

install.packages("ape")   #installing packages

install.packages("phangorn")

install.packages("seqinr")

library(ape)

library(phangorn)

library(seqinr)


### UPGMA tree for concatenated alignments

TwoGene <- read.dna("bothgenes.fasta", format="fasta") #reads fasta alignment

TwoGene_phyDat <- phyDat(TwoGene, type = "DNA", levels = NULL) #calls functions of ape and phangorn packages

print(TwoGene_phyDat) #test to count character numbers
 
dm <- dist.ml(TwoGene_phyDat)  #creates distance matrix for maximum likelihood

treeUPGMA <- upgma(dm) #creates UPGMA tree


layout(matrix(c(.01,.02), 2, 1), height=c(.01,.02)) #plots UPGMA tree
par(mar = c(0,0,2,0)+ .00001)
plot(treeUPGMA, main="UPGMA")

title("T.sirtalis UPGMA Tree")


############# Maximum Likelihood tree for concatenated alignments

TwoGene <- read.dna("bothgenes.fasta", format="fasta") #reads fasta alignment

TwoGene_phyDat <- phyDat(TwoGene, type = "DNA", levels = NULL) #calls functions of ape and phangorn packages

print(TwoGene_phyDat) #test to count character numbers

dm <- dist.ml(TwoGene_phyDat) #creates distance matrix for maximum likelihood
treeNJ <- NJ(dm) #creates Neighbor Joining tree

treeNJ$edge.length[which(treeNJ$edge.length <= 0)] <- 0.0000001  #sets parameters to >0 
#(this prevents errors while rounding all zeros up to a very small #)

fit = pml(treeNJ, data=TwoGene_phyDat) #computes likelihood for this tree

AIC(fitJC)  #using tests to see which model fits the data the best
3888.323
AIC(fitGTR)
3716.371  #note that GTR returns lowest value, so we will use this
AICc(fitGTR)
3720.064
BIC(fitGTR)
3944.17


fitGTR <- update(fit, k=4, inv=0.2) #fitting GTR model to our tree
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR

bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
                   control = pml.control(trace = 0)) #calculating bootstrap values for branches
par(mfrow=c(1.5,1)) 
par(mar=c(1,1,3,1))

plotBS(midpoint(fitGTR$tree), bs, p = 50, type="p") #plotting bootstrap values to tree

title("T.sirtalis Maximum Likelihood Tree")


