# Figure 4: Phylogenetic analysis

#### Phylogeny - log model
library(seqinr)
library(phangorn)
library(phytools)

## Import a tree
my.tree = read.tree(file = "pruned tree")

## Import alignment data
my.alignment = read.alignment(file="gtsB alignment.fasta", format = "fasta")
## Convert it to PhyDat format
my.phyDat = as.phyDat(my.alignment, type = "DNA", return.index = T)


N = 100 # number of times to randomly resolve the multichotomy

for ( n in 1:N ){
  print(n)
  # Randomly resolve the tree
  my.tree <- multi2di(my.tree, random = TRUE)
  
  ## Calculate ancestral states at each node (here I'm using the parsimony method, but there are other ways to do this)
  anc.acctran <- ancestral.pars(my.tree, my.phyDat, type = "ACCTRAN")
  
  ## Number of divergent sites in the sequence alignment
  nsites = nrow(anc.acctran[1][[1]]) 
  
  ## Go through each divergent site in the alignment and count the number of evolutionary events on the tree
  for ( i in 1:nsites){
    # print(i)
    all.changes = matrix(data = 0, nrow = 4, ncol = 4)
    rownames(all.changes) = c("a", "c", "g", "t") ## ancestral state (state 0)
    colnames(all.changes) = rownames(all.changes) ## derived state (state 1)
    
    ## Check each branch for a substitution
    for (e in 1 : nrow(my.tree$edge) ) {
      
      state0 = anc.acctran[my.tree$edge[e,1]][[1]][i,] # state at the start of the edge
      state1 = anc.acctran[my.tree$edge[e,2]][[1]][i,] # state at the end of the edge
      
      ## Compare the states and count the transitions (if the state are probabilities, the counts can be a non-integer)
      changes = matrix(data = 0, nrow = 4, ncol = 4)
      changes[, which(state1!=0)] = state0
      for( p in 1:4 ){
        changes[,p] = changes[,p]*state1[p]
      }
      
      all.changes = all.changes + changes
      
    }
    
    all.changes[is.nan(all.changes)] = 0
    if(i == 1) { 
      all.edge.changes = list(all.changes) 
    }
    
    if(i > 1) { 
      all.edge.changes = append(all.edge.changes,list(all.changes))
    }
    
  }
  if(n == 1) {
    mean.edge.changes = all.edge.changes
  }
  if(n > 1) {
    for( j in 1:nsites ){
      mean.edge.changes[[j]] = mean.edge.changes[[j]] + all.edge.changes[[j]]
    }
  }
}




mut.data = read.csv(file = "mutant fitness for comparison.csv")

mut.data$wt.nt = as.character(mut.data$wt.nt)
mut.data$wt.nt[mut.data$wt.nt=="A"] = 1
mut.data$wt.nt[mut.data$wt.nt=="C"] = 2
mut.data$wt.nt[mut.data$wt.nt=="G"] = 3
mut.data$wt.nt[mut.data$wt.nt=="T"] = 4
mut.data$wt.nt = as.numeric(mut.data$wt.nt)

mut.data$mut.nt = as.character(mut.data$mut.nt)
mut.data$mut.nt[mut.data$mut.nt=="A"] = 1
mut.data$mut.nt[mut.data$mut.nt=="C"] = 2
mut.data$mut.nt[mut.data$mut.nt=="G"] = 3
mut.data$mut.nt[mut.data$mut.nt=="T"] = 4
mut.data$mut.nt = as.numeric(mut.data$mut.nt)


no.changes = rep(0, nrow(mut.data))
for( m in 1:nrow(mut.data) ){
  no.changes[m] = all.edge.changes[[attributes(anc.acctran)$index[mut.data$nt.pos[m]]]][mut.data$wt.nt[m], mut.data$mut.nt[m]]/sum(all.edge.changes[[attributes(anc.acctran)$index[mut.data$nt.pos[m]]]][mut.data$wt.nt[m], ])
}
mut.data$no.changes = no.changes


all.changes = rep(0, nrow(mut.data))
for( m in 1:nrow(mut.data) ){
  all.changes[m] = sum(all.edge.changes[[attributes(anc.acctran)$index[mut.data$nt.pos[m]]]][, ]*(upper.tri(all.edge.changes[[attributes(anc.acctran)$index[mut.data$nt.pos[m]]]][, ])+lower.tri(all.edge.changes[[attributes(anc.acctran)$index[mut.data$nt.pos[m]]]][, ])))
}


y = no.changes[order(mut.data$mean.w)]
x = sort(mut.data$mean.w)
z = as.factor(mut.data$effect[order(mut.data$mean.w)])
numeric.z = (z=='S')+0

colors = c( (col=rgb(0,0,1,1/2)), (col=rgb(1,0,0,1/2)), "grey50" )
plot(x, y, pch = 19, col = colors[(z!='S')+1], xlab = "Relative fitness in SBW25", 
     ylab = "Probability of observing mutation")
abline(v = 1, lty = 3)


mod <- nls(y ~ exp(a + b * x), start = list(a = 0, b = 0), control = list(maxiter = 500))
summary(mod)

lines(x, predict(mod, list(x)), lwd = 2, col = "black")


# binary: mutation observed or not
y.bin = y>0
plot(x, y.bin, pch = 19, col = colors[as.numeric(z)], xlab = "Relative fitness", ylab = "Mutation observed (T/F)")
mod.bin = glm(y.bin ~ x, family = "binomial")
summary(mod.bin)
lines(x, predict(mod.bin, list(x), type = "response"), lwd = 2, col = "black", lty = 2)
abline(v = 1, lty = 3)
