############################################################################################
## cMonkey - version 4, Copyright (C) Christopher L Plaisier, Institute for Systems Biology
##                      cplaisier@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###########################################################################################

## get.enrichment.sets
## Reads in the sets for enrichment analysis during a cMonkey run.
## Parameters:
##  set.types <- types of sets to use for calculating enrichment as a named vector where the key is the name of the set (what it will be called by cMonkey) and the value is the location of a csv file of the sets to be read in by cMonkey. Csv files can be 2 columns: first column is the element name, and the second column a gene id matching to the ratios matrix row ids where equal weighting is assumed, or 3 columns:  first column is the element name, and the second column a gene id matching to the ratios matrix row ids, and a weight for the connection to the element.
## Returns a list indexed by set.type of a list indexed by set names of a named vector of gene keys with weights as the values.
get.enrichment.sets <- function(set.types=NA) {
  # Make sure that there is something to intialize
  if ( ! set.types=='NA') {
    stop('Error in get.enrichment.sets -> No set.types.')
  }
  # Read in file
  enrichment.sets <- list()
  for (set.type.i in names(set.types)) {
    d1 = read.csv(file=set.types[[set.type.i]]['file'], header=F)
    d2 = list()
    for(i in unique(d1[,1])) {
      tmp.vector.names <- d1[which(d1[,1]==i),2]
      # If no weights are in file then default them to 1
      if(ncol(d1)==2) {
        tmp.named.vector <- rep(1,length(tmp.vector.names))
      # Otherwise set them as the value in the file
      } else if(ncol(d1)==3) {
        tmp.named.vector <- as.numeric(d1[which(d1[,1]==i),3])
      }
      names(tmp.named.vector) <- tmp.vector.names
      d2[[i]] <- tmp.named.vector
    }
    enrichment.sets[[factor.i]] <- d2
  }
  return(enrichment.sets)
}

## get.set.enrichment.scores
## Calculates enrichment scores for all sets and returns the best scoring element for each bicluster.
## Parameters:
##  k <- can either be the cluster index or an array or probe ids.
##  set <- the name of the set to be used to calculate enrichment.
get.set.enrichment.scores <- function( k, set=NA) {
  # Make sure that things are right for anlaysis (took some of this from Daves Network scoring function)
  if ( length( k ) <= 0 ) {
    stop('get.set.enrichment.scores <- k is not set.')
  }
  if ( set==NA && set %in% enrichment.sets) {
    stop('get.set.enrichment.scores <- set==NA or set not present in enrichment.sets: need to give a valid set.')
  }
  if ( is.numeric( k[ 1 ] ) ) {
    k.rows <- get.rows( k )
  } else {
    k.rows <- k
  }
  all.rows <- attr( ratios, "rnames" )
  if ( length( rows ) < 1 ) {
    return( rep( NA, length( all.rows ) ) )
  }

  ###   Hypergeometric - Choose best set   ###
  # Now calculate the hypergeometric enrichment p-values for all elements in the specified set
  enrichment.p.values <- c(1:length(names(enrichment.sets[[set]])))
  names(enrichment.p.values) <- names(enrichment.sets[[set]])
  # Get a list of all the genes for the enrichment set, the background set
  all.genes.enrichment.set <- unique(as.vector(unlist(enrichment.sets[[set]])))
  # This screens out genes in the cluster that are not in the background set
  cluster.genes <- k.rows[k.rows %in% all.genes.enrichment.set]
  for(i in names(enrichment.p.values)) {
    # Calculate hypergeometric for an element of a set
    element.genes <- unique(enrichment.sets[[set]][[i]])
    overlap.genes <- intersect(cluster.genes,element.genes)
    # q <- overlap of the genes for the element with the genes from the bicluster, m <- genes for the element, n <- background set of genes minus genes for element, k <- genes in the bicluster (fitting, no?) 
    enrichment.p.values[i] <- phyper(length(overlap.genes),length(element.genes),length(all.genes.enrichment.set)-length(element.genes),length(cluster.genes),lower.tail=F)
  }
  # Get the minimum p-value element for the enrichment analysis (assuming for the time being that no such thing as a tie). Could put a dealy here that collects all tieing top enrichments and then randomly selects the top. Adding to the stocasticity of cMonkey. Likely issues like these would sort themselves out over a run.
  min.element <- names(enrichment.p.values)[which(enrichment.p.value==min(enrichment.p.value))]
  
  ###   Scoring function   ###
  # Now set up the scores for all the rows
  if (mean(  #  1. Genes from the current bicluster and from the minimum element (overlap.genes) are given full log10 of the minimum element p-value
  #  2. Genes from the minimum element and not in the current bicluster (element.genes) are given half the log10 of the minimum element p-value
  #  3. Genes not from either the current bicluster nor the minimum element are given a zero score (i.e. no information)
  element.genes <- unique(enrichment.sets[[set]][[min.element]]
  overlap.genes <- intersect(cluster.genes,element.genes)
  # If done in this order fist all default to zero, then minimum element genes are set to half score, and finally overlapping genes are set at full score
  scores <- rep(length(all.rows),0) # Default all to zero
  names(scores) <- all.rows
  scores[element.genes] <- log10(enrichment.p.values[min.element])/2 # In minimum scoring element, but not
  scores[overlap.genes] <- log10(enrichment.p.values[min.element]) # In bicluster and from minimum scoring element
  # Attach the minimum element name so this can be stored
  attr(scores,'min.element') <- min.element
  return( scores )
}

## seed.cluster.enrichment
## Calculates enrichment scores for all sets and returns the best scoring element for each bicluster.
## Parameters:
##  k <- can either be the cluster index or an array or probe ids.
##  set <- the name of the set to be used to calculate enrichment.
#seed.cluster <- function( k.clust, seed.method='rnd', col.method='rnd') {
#  
#}

