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
##  set.types <- types of sets to use for calculating enrichment as a list of set types with a named vector where the keys are file and cutoff for the set (what it will be called by cMonkey) and the value is the location of a csv file of the sets to be read in by cMonkey. Csv files can be 2 columns: first column is the element name, and the second column a gene id matching to the ratios matrix row ids where equal weighting is assumed, or 3 columns:  first column is the element name, and the second column a gene id matching to the ratios matrix row ids, and a weight for the connection to the element.
##  e.g. set.types <- list(miRNA = c(file='miRNAset.csv',cutoff=0.05), TF = c(file='tfSet.csv',cutoff=0.01))
## Returns a list indexed by set.type of a list indexed by set names of a named vector of gene keys with weights as the values.
get.enrichment.sets <- function(set.types=NA) {
  cat('Reading in enrichment sets...')
  # Make sure that there is something to intialize
  if ( is.na(set.types) ) {
    stop('Error in get.enrichment.sets -> No set.types.')
  }
  # Load data structure
  enrichment.sets <- list()
  for (set.type.i in names(set.types)) {
    # Read in set.type file of set identities, can be accompanied by weights (3rd column)
    # 2 column e.g. (no header):
    # hsa-miR-i1,gene1
    # hsa-miR-i2,gene2
    #      ...
    # hsa-miR-in,genem
    #
    # or
    #
    # 3 column e.g. (no header):
    # hsa-miR-i1,gene1,weight1
    # hsa-miR-i2,gene2,weight2
    #      ...
    # hsa-miR-im,genen,weighto
    d1 = read.csv(file=set.types[[set.type.i]]['file'], header=F)
    d2 = list()
    # Build a cutoff array to attach as an attribute for discretizing the data for hypergeometric p-value calculation
    if ( ! is.na(set.types[[set.type.i]]['cutoff']) ) {
      # A two column file with cutoffs for each set:  Column 1 = set name, Column 2 = cutoff
      # e.g. hsa-miR-1,0.05
      #      hsa-miR-2,0.05
      #            ...
      #      hsa-miR-i,j
      tmp1 <- read.csv(file=set.types[[set.type.i]]['cutoff'], header=F)
      cutoffs <- tmp1[,2]
      names(cutoffs) <- tmp1[,1]
    } else {
      # Set them to a default cutoff of 0.5, assuming that already discretized
      cutoffs <- rep('discrete',length(unique(d1[,1])))
      names(cutoffs) <- unique(d1[,1])
    }
    for(set.i in unique(d1[,1])) {
      tmp.vector.names <- toupper(d1[which(d1[,1]==set.i),2])
      # If no weights are in file then default them to 1
      if(ncol(d1)==2) {
        tmp.named.vector <- rep(1,length(tmp.vector.names))
      # Otherwise set them as the value in the file
      } else if(ncol(d1)==3) {
        tmp.named.vector <- as.numeric(d1[which(d1[,1]==i),3])
      }
      # Add gene names to weights
      names(tmp.named.vector) <- tmp.vector.names
      # Add cutoff to vector so we can discretize it if not already discretized
      attr(tmp.named.vector,'cutoff') <- cutoffs[set.i]
      d2[[set.i]] <- tmp.named.vector
    }
    enrichment.sets[[set.type.i]] <- d2
  }
  cat('Done.')
  return(enrichment.sets)
}

## get.set.enrichment.scores
## Calculates enrichment scores for all sets and returns the best scoring set for each bicluster.
## Parameters:
##  k <- can either be the cluster index or an array or probe ids.
##  set.type <- the name of the set type to be used to calculate enrichment.
get.set.enrichment.scores <- function( k, set.type=NA, ...) {
  # Make sure that things are right for anlaysis (took some of this from Daves Network scoring function)
  if ( length( k ) <= 0 ) {
    stop('get.set.enrichment.scores <- k is not set.')
  }
  if ( is.na(set.type) && set.type %in% enrichment.sets) {
    stop('get.set.enrichment.scores <- set.type==NA or set not present in enrichment.sets: need to give a valid set.')
  }
  if ( is.numeric( k[ 1 ] ) ) {
    k.rows <- get.rows( k )
  } else {
    k.rows <- k
  }
  all.rows <- attr( ratios, "rnames" )
  if ( length( all.rows ) < 1 ) {
    return( rep( NA, length( all.rows ) ) )
  }

  ###   Hypergeometric - Choose best set   ###
  # Now calculate the hypergeometric enrichment p-values for all sets in the specified set type
  # Get a list of all the genes for the set type, the background set
  # This screens out genes in the cluster that are not in the background set
  cluster.genes <- k.rows[k.rows %in% all.genes.enrichment.set[[set.type]]]
  
  # Calculate hypergeometric for a set
  if( attr(enrichment.sets[[set.type]][[1]],'cutoff')=='discrete' ) {
    # Already a discrete set
    set.genes <- sapply(names(enrichment.sets[[set.type]]), function(x) { names(enrichment.sets[[set.type]][[x]]) } )
  } else {
    # Need to discretize into a set for enrichment
    # Assuming that values are from a similarity metric i.e. higher values indicate more confidence for inclusion in set
    set.genes = sapply(names(enrichment.sets[[set.type]]), function(x) { set.i <- enrichment.sets[[set.type]][[x]]; return(names(set.i)[which(as.numeric(set.i)>=attr(set.i,'cutoff'))]) } )
  }
  overlap.genes <- sapply(1:length(set.genes), function(x) {length(intersect(cluster.genes,set.genes[[x]]))})
  # q <- overlap of the genes for the set with the genes from the bicluster, m <- genes for the set, n <- background set of genes minus genes for set, k <- genes in the bicluster (fitting, no?) 
  enrichment.p.values <- phyper(overlap.genes,sapply(set.genes,length),rep(length(all.genes.enrichment.set[[set.type]]),length(set.genes))-sapply(set.genes,length),rep(length(cluster.genes),length(set.genes)),lower.tail=F)
  names(enrichment.p.values) <- names(enrichment.sets[[set.type]])
  # Get the minimum p-value set for the enrichment analysis (assuming for the time being that no such thing as a tie). Could put a dealy here that collects all tieing top enrichments and then randomly selects the top. Adding to the stocasticity of cMonkey. Likely issues like these would sort themselves out over a run.
  min.set <- (names(enrichment.p.values)[which(enrichment.p.values==min(enrichment.p.values))])[1]
  
  ###   Scoring function   ###
  # Now set up the scores for all the rows
  scores <- rep(NA,length(all.rows)) # Default all to zero
  names(scores) <- all.rows
  if (attr(enrichment.sets[[set.type]][[min.set]],'cutoff')=='discrete') {
    #  1. Genes from the current bicluster and from the minimum set (overlap.genes) are given full log10 of the minimum set p-value
    #  2. Genes from the minimum set and not in the current bicluster (set.genes) are given half the log10 of the minimum set p-value
    #  3. Genes not from either the current bicluster nor the minimum set are given a zero score (i.e. no information)
    min.genes <- names(enrichment.sets[[set.type]][[min.set]])
    min.genes <- min.genes[min.genes %in% all.rows]
    overlap.genes <- intersect(cluster.genes,min.genes)
    # If done in this order fist all default to zero, then minimum set genes are set to half score, and finally overlapping genes are set at full score
    scores[min.genes] <- 0.5 # In minimum scoring set, but not in the bicluster
    scores[overlap.genes] <- 1 # In bicluster and from minimum scoring set
  } else {
    # Combine the p-value with the quantitative weights for the set with the minimum p-value
    set.min.weights <- enrichment.sets[[set.type]][[set]]
    set.min.weights <- set.min.weights[names(set.min.weights) %in% all.rows]
    # Insert the weights and scale the them to be between zero and one
    scores[names(set.min.weights)] <- as.numeric(set.min.weights)-min(set.min.weights)
    scores[names(set.min.weights)] <- scores[names(set.min.weights)]/max(scores[names(set.min.weights)]) 
  }
  scores <- enrichment.p.values[min.set] / scores
  scores[is.na(scores)] <- 1
  scores <- log10(scores)
  # Attach the minimum set name so this can be stored
  attr(scores,'min.set') <- min.set
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

