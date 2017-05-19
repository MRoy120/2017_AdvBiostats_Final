## This example of our Bayes test for model correctness uses the
## set of 5 five-variable models (fivevariable1.dia - fivevariable5.dia),
## originally developed by Dambacher et al. (2003, Ecologial
## Modelling).
library("XML")

source("dia.r")
source("community.r")

files <- paste("fivevariable", 1:5, ".dia", sep="")

## Read model specification
edges <- lapply(files, QPress::model.dia)
names(edges) <- files

## Examine unweighted adjacency matrices
A <- lapply(edges, adjacency.matrix)
A

## Function to generate community matrices
s <- lapply(edges,community.sampler)

## Function to estimate posterior probabilities of model correctness, 
## assuming uniform priors
Bayes.test <- function(edges,samplers,model,perturb,monitor,n.samples=1000) {

    ## Simulate a stable weight matrix from the selected model
    W <- s[[model]]$community()
    while(!stable.community(W)) {
        W <- s[[model]]$community()
    }

    ## Determine the outcome of the press perturbation
    impact <- press.impact(edges[[model]],perturb=perturb,monitor=monitor)
    observed <- signum(impact(W))
    names(observed) <- names(monitor)

    prop <- double(length(edges))
    
    ## Loop over models
    for(k in 1:length(edges)) {
        n.stable <- 0
        n.valid <- 0

        ## Function to validate press condition
        press <- press.validate(edges[[k]],perturb=perturb, monitor=observed)

        while(n.stable < n.samples) {

            ## Sample community matrix
            W <- s[[k]]$community()

            ## Check stability
            if(!stable.community(W)) next
            n.stable <- n.stable+1

            ## Check press condition
            if(press(W)) n.valid <- n.valid+1
        }
        prop[k] <- n.valid/n.stable
    }
    list(p=prop/sum(prop),
         perturb=perturb,
         monitor=monitor)
}


n.sims <- 10
## Pick one node to perturb and 2 others to monitor
n.perturb <- 1
n.monitor <- 2
P <- array(0,c(length(edges),n.sims,length(edges)))
for(model in 1:5) {
  for(k in 1:n.sims) {
    ## Random perturbation and monitoring
    nodes <- levels(edges[[model]]$From)
    perturb <- rep(1,length=n.perturb)
    names(perturb) <- sample(nodes,length(perturb))
    monitor <- rep(1,length=n.monitor)
    names(monitor) <- sample(setdiff(nodes,names(perturb)),length(monitor))

    ## Compute posterior probabilities
    fit <- Bayes.test(edges,s,model,perturb,monitor,n.samples=1000)
    P[,k,model] <- fit$p
  }
}

## Generate a confusion matrix
CM <- matrix(0,dim(P)[3],dim(P)[1])
for(model in 1:(dim(P)[3]))
  CM[model,] <- tabulate(unlist(apply(P[,,model],2,function(x) which(x==max(x)))),dim(P)[1])
## Columns: correct model, Rows: highest ranking model
t(CM)

## Determine how often the correct model is most likely, second most
## likely etc.
Q <- matrix(0,dim(P)[3],dim(P)[1])
for(model in 1:(dim(P)[3]))
  Q[model,] <- tabulate(apply(-P[,,model],2,
	function(x) which(order(x,c(model,rep(0,dim(P)[1]-1)))==model)),dim(P)[1])

## Plot how often the correct model is most likely, second most likely
## etc, separately for each correct model and overall.
opar <- par(mfrow=c(3,2))
for(model in 1:nrow(Q)) {
  p <- Q[model,]
  names(p) <- 1:ncol(Q)
  barplot(p/sum(p),main=model,ylim=c(0,1))
}
p <- colSums(Q)
names(p) <- 1:ncol(Q)
barplot(p/sum(p),main="all",ylim=c(0,1))

## Use clustering to show model similarity
opar <- par(mfrow=c(3,2))
method <- "complete"
for(model in 1:(dim(P)[3])) {
  plot(hclust(dist(P[,,model]),method=method),ann=F)
  title(model)
}
plot(hclust(dist(matrix(P,dim(P)[1],prod(dim(P)[-1]))),method=method),ann=F)
title("All")
par(opar)

#Test out Bayes Test
Bayes.test(edges, samplers, model, 1, 2, n.samples=1000)







