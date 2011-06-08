###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

## This func. is the meat of cmonkey - run one iteration and update the scores/parameters.
## This is called from cmonkey() - see there for its usage.
## Warning - this relies on lots of global variables, that are updated inside the function and
##   the updated variables are returned. Again, see cmonkey() for its usage.

## BIG TODO: make it so all big.memory are updated in-place
cmonkey.one.iter <- function( env, dont.update=F, ... ) {  
  ##if ( ! exists( "clust.changed", envir=env ) ) env$clust.changed <- rep( FALSE, k.clust )
  if ( ! exists( "row.memb", envir=env ) ) {
    env$row.memb <- t( apply( row.membership[], 1, function( i ) 1:k.clust %in% i ) )
    env$col.memb <- t( apply( col.membership[], 1, function( i ) 1:k.clust %in% i ) )
    env$row.memb <- matrix.reference( env$row.memb, backingfile="row.memb", backingpath=env$cmonkey.filename )
    env$col.memb <- matrix.reference( env$col.memb, backingfile="col.memb", backingpath=env$cmonkey.filename )
  } else {
    env$row.memb[,] <- t( apply( row.membership[], 1, function( i ) 1:k.clust %in% i ) )
    env$col.memb[,] <- t( apply( col.membership[], 1, function( i ) 1:k.clust %in% i ) )
  }

  tmp1 <- get.all.scores( ... )
  env$row.scores <- tmp1$r ##[,];
  env$col.scores <- tmp1$c
  env$mot.scores <- tmp1$m
  env$meme.scores <- tmp1$ms
  env$net.scores <- tmp1$n
  env$set.scores <- tmp1$se
  if ( ! is.null( tmp1$cns ) ) env$cluster.net.scores <- tmp1$cns
  if ( ! is.null( tmp1$cses ) ) env$cluster.set.scores <- tmp1$cses
  rm( tmp1 )
  
  # Quantile normalize the scores
  tmp2 <- get.combined.scores( quant=T )
  env$r.scores <- tmp2$r
  env$c.scores <- tmp2$c
  env$m.scores <- tmp2$m
  env$n.scores <- tmp2$n
  env$s.scores <- tmp2$s
  if ( ! is.null( env$row.scores ) ) attr( env$row.scores, "changed" ) <- FALSE
  if ( ! is.null( env$col.scores ) ) attr( env$col.scores, "changed" ) <- FALSE
  if ( ! is.null( env$mot.scores ) ) attr( env$mot.scores, "changed" ) <- FALSE
  if ( ! is.null( env$net.scores ) ) attr( env$net.scores, "changed" ) <- FALSE
  if ( ! is.null( env$set.scores ) ) attr( env$set.scores, "changed" ) <- FALSE
  
  if ( length( tmp2$scalings ) > 0 ) {
    env$row.scaling[ iter ] <- tmp2$scalings[ "row" ]
    env$mot.scaling[ iter ] <- tmp2$scalings[ "mot" ]
    env$net.scaling[ iter ] <- tmp2$scalings[ "net" ]
    env$set.scaling[ iter ] <- tmp2$scalings[ "set" ]
  }
  rm( tmp2 )

  ## Fuzzify scores a bit for stochasticity! (fuzz should be between 0.2 and 0 (decreasing with iter)
  if ( fuzzy.index[ iter ] > 1e-5 ) {
    env$r.scores[,] <- env$r.scores[,] +
      rnorm( length( env$r.scores[,] ), sd=sd( env$r.scores[,][ row.memb[,] == 1 ], na.rm=T ) * fuzzy.index[ iter ] )
    if ( ! is.null( env$c.scores ) ) env$c.scores[,] <- env$c.scores[,] +
      rnorm( length( env$c.scores[,] ), sd=sd( env$c.scores[,][ col.memb[,] == 1 ], na.rm=T ) * fuzzy.index[ iter ] )
  }

  tmp3 <- get.density.scores( ks=1:k.clust ) ## r.scores, col.scores, 
  env$rr.scores <- tmp3$r
  env$cc.scores <- tmp3$c
  rm( tmp3 )

  ## OUTPUT
  if ( iter %in% stats.iters ) { 
    env$clusterStack <- get.clusterStack( ks=1:k.clust ) 
    env$stats <- rbind( stats, get.stats() )
    cat( organism, as.matrix( stats[ nrow( stats ), ] ), '\n' )
  } else {
    cat( sprintf( "==> %04d %.3f %.3f %.3f\n", iter, ##%5d sum( row.memb != old.row.memb, na.rm=T ),
                 mean( row.scores[,][ row.memb[,] == 1 ], na.rm=T ), ##mean( col.scores[ col.memb ], na.rm=T, trim=0.05 ),
                 if ( ! is.null( mot.scores ) ) mean( mot.scores[,][ row.memb[,] == 1 & mot.scores[,] < 0 ], na.rm=T, trim=0.05 )
                 else NA,
                 if ( ! is.null( net.scores ) ) mean( net.scores[,][ row.memb[,] == 1 ##& net.scores < 0
                                                                 ], na.rm=T, trim=0.05 ) else NA
                 if ( ! is.null( set.scores ) ) mean( set.scores[,][ row.memb[,] == 1], na.rm=T, trim=0.05 ) else NA ) )
  }
    
  ## NEW - will it work? -- help shrink big clusters, grow small clusters, both in rows and cols
  size.compensation.func.rows <- function( n ) exp( -n / ( attr( ratios, "nrow" ) * n.clust.per.row / k.clust ) )
  size.compensation.func.cols <- function( n ) exp( -n / ( attr( ratios, "ncol" ) * n.clust.per.col / k.clust ) )
  for ( k in 1:k.clust ) {
    tmp4 <- sum( row.memb[ ,k ] )
    if ( tmp4 > 0 ) env$rr.scores[ ,k ] <- env$rr.scores[ ,k ] * size.compensation.func.rows( tmp4 ) 
    else env$rr.scores[ ,k ] <- env$rr.scores[ ,k ] * size.compensation.func.rows( cluster.rows.allowed[ 1 ] )
    if ( ! is.null( env$cc.scores ) ) {
      tmp5 <- sum( col.memb[ ,k ] )
      if ( tmp5 > 0 ) env$cc.scores[ ,k ] <- env$cc.scores[ ,k ] * size.compensation.func.cols( tmp5 ) 
      else env$cc.scores[ ,k ] <- env$cc.scores[ ,k ] * size.compensation.func.cols( attr( ratios, "ncol" ) / 10 )
    }
  }
  
  ## Fuzzify it along the lines of fuzzy c-means clustering
  ##   -- see http://en.wikipedia.org/wiki/Data_clustering#Fuzzy_c-means_clustering
  ## No -- it doesnt affect things - same ordering (and updated memberships are based on ordering)
  ##   but should use these scores to weight the centroids that are selected in the next iteration.
  ##   for this, fuzz should vary between e.g. 10 and 2
  ##rr.scores <- ( rr.scores / sum( rr.scores, na.rm=T ) )^( 2 / ( fuzz - 1 ) )
  ##cc.scores <- ( cc.scores / sum( cc.scores, na.rm=T ) )^( 2 / ( fuzz - 1 ) )

  if ( ! dont.update ) {
  
    ## Make a matrix of m[i,k] = whether row/col i is in cluster k
    if ( exists( "row.membership" ) ) { 
      env$old.row.membership <- row.membership 
      env$old.col.membership <- col.membership 
    }

    tmp6 <- get.updated.memberships() ## rr.scores, cc.scores )
    env$row.membership <- tmp6$r
    env$col.membership <- tmp6$c
#ifndef PACKAGE
    tmp6 <- filter.updated.memberships()
#endif
    if ( ! is.null( tmp6 ) ) { env$row.membership <- tmp$r; env$col.membership <- tmp$c }
    
#ifndef PACKAGE
    if ( TRUE && iter %% 10 == sample( 0:9, 1 ) && iter < n.iter * 0.8 ) { ## Dont do this near the end of the run
      ## Merge clusters that are dupes and then re-seed the newly-emptied ones
      if ( merge.cutoffs[ "n" ] > 0 && merge.cutoffs[ "cor" ] < 1 ) {
        tmp.m <- merge.cutoffs[ "n" ]; if ( tmp.m < 1 && runif( 1 ) <= tmp.m ) tmp.m <- 1
        if ( tmp.m >= 1 ) {
          tmp7 <- consolidate.duplicate.clusters( scores=r.scores, cor.cutoff=merge.cutoffs[ "cor" ],
                                                n.cutoff=tmp.m, motif=F )
          env$row.membership <- tmp7$r; env$meme.scores <- tmp7$ms
        }
      }
      
      ## Re-seed that are too small (including those zero-ed out by the dupe-consolidation above)
      tmp8 <- re.seed.empty.clusters( toosmall.r=cluster.rows.allowed[ 1 ],
                                    toosmall.c=min(attr(ratios,"ncol")/10,50),
                                    n.r=cluster.rows.allowed[ 1 ] * 2, n.c=attr(ratios,"ncol")/5 )
      env$row.membership <- tmp8$r; env$col.membership <- tmp8$c; env$meme.scores <- tmp8$ms 
    }
  ## }
#endif
    
    ## PLOTTING
    if ( ! is.na( plot.iters ) && iter %in% plot.iters ) {
      env$clusterStack <- get.clusterStack( ks=1:k.clust ) 
      try( plotStats( iter, plot.clust=env$favorite.cluster(), new.dev=T ), silent=T ) ## Can be set for your given organism
    }
  
    if ( exists( "cm.func.each.iter" ) ) try( cm.func.each.iter(), silent=T ) ## User-defined func. to run each iteration
    
    ## Allow temp source file to be sourced (e.g. to change a param in the middle of a run, or print or plot or save
    ## some intermediate results). If first line of file is '## QUIET' then this is done quietly.
    if ( any( cm.script.each.iter != "" ) ) {
      for ( f in cm.script.each.iter ) {
        if ( file.exists( f ) && file.info( f )$size > 1 ) {
          tmp9 <- readLines( f )
          if ( all( substr( tmp9, 1, 1 ) == "#" ) ) next ## All commented-out code
          if ( tmp9[ 1 ] != "## QUIET" ) cat( "Sourcing the script '", f, "' ...\n", sep="" )
          try( source( f, echo=tmp9[ 1 ] != "## QUIET", local=T ), silent=T )
        }
      }
    }
  }

#ifndef PACKAGE
  ##if ( big.memory == TRUE || big.memory > 0 )
  ##ffify.env( env ) ## Big matrices and lists go to filebacked version
#endif
  
  ## Note: with the above code, when the env is saved via save.image(), all ff obj's are "closed" but their
  ##   filestores still exist, so you can "open.ff(x)" each of them after the env is re-loaded,
  ##   and that will reconnect them with their files.
  
  if ( get.parallel()$mc ) {
    if ( getDoParName() == "doMC" ) { ##require( multicore, quietly=T ) ) { ## Clean up any multicore spawned processes (as doc'ed in mclapply help)
      chld <- multicore::children()
      if ( length( chld ) > 0 ) { try( { multicore::kill( chld ); tmp10 <- multicore::collect( chld ) }, silent=T ) }
    } else if ( getDoParName() == "doSNOW" && "data" %in% ls( pos=foreach:::.foreachGlobals ) ) {
      cl <- get( "data", pos=foreach:::.foreachGlobals ) ## Tricky, eh?
      if ( ! is.null( data ) ) stopCluster( cl )
    }
  }
  
  if ( ! dont.update ) env$iter <- env$iter + 1
  invisible( env )
}

## TODO: store/plot indiv. resids and network scores for each separate component
get.stats <- function( mean.func=median ) {
  ##if ( ! exists( "row.memb" ) ) row.memb <- t( apply( row.membership, 1, function( i ) 1:k.clust %in% i ) )
  ##if ( ! exists( "col.memb" ) ) col.memb <- t( apply( col.membership, 1, function( i ) 1:k.clust %in% i ) )
  changed <- NA
  if ( ! is.null( old.row.membership ) )
    changed <- sum( row.memb[,] != t( apply( old.row.membership, 1, function( i ) 1:k.clust %in% i ) ), na.rm=T )

  cs <- as.list( clusterStack )
  resids <- sapply( cs, "[[", "resid" )
  if ( is.matrix( resids ) ) resids <- apply( resids, 1, function( r ) mean.func( r[ r != 1.0 ], na.rm=T ) )
  else resids <- mean.func( resids[ resids != 1.0 ], na.rm=T )
  p.clusts <- sapply( cs, "[[", "p.clust" )
  if ( is.matrix( p.clusts ) ) p.clusts <- apply( p.clusts, 1, mean.func, na.rm=T ) 
  else p.clusts <- mean.func( p.clusts, na.rm=T )
  out <- data.frame( iter=iter, changed=changed,
             row.scores=mean( row.scores[,][ row.memb[,] == 1 ], na.rm=T, trim=0.05 ),
             col.scores=mean( col.scores[,][ col.memb[,] == 1 ], na.rm=T, trim=0.05 ),
             mot.scores=if ( ! is.null( mot.scores ) ) mean.func( mot.scores[,][ row.memb[,] == 1 ], na.rm=T ) else NA, 
             net.scores=if ( ! is.null( net.scores ) ) mean( net.scores[,][ row.memb[,] == 1 ], na.rm=T, trim=0.05 ) else NA,
             set.scores=if ( ! is.null( set.scores ) ) mean( set.scores[,][ row.memb[,] == 1 ], na.rm=T, trim=0.05 ) else NA,
             resid=weighted.mean( resids, row.weights, na.rm=T ),
             nrow=mean.func( sapply( cs, "[[", "nrows" ), na.rm=T ),
             ncol=mean.func( sapply( cs, "[[", "ncols" ), na.rm=T ),
             p.clust=if ( ! all( is.na( p.clusts ) ) ) weighted.mean( p.clusts, mot.weights, na.rm=T ) else NA
             )
  if ( length( resids ) > 1 ) for ( i in names( resids ) ) {
    out <- cbind( out, resids[ i ] )
    names( out )[ ncol( out ) ] <- paste( "resid", i, sep="." )
  }
  if ( length( p.clusts ) > 1 ) for ( i in names( p.clusts ) ) {
    out <- cbind( out, p.clusts[ i ] )
    names( out )[ ncol( out ) ] <- paste( "p.clust", i, sep="." )
  }
  ## net.scores weighted by n.rows of bicluster
  if ( length( networks ) > 1 ) {
    for ( i in names( net.weights ) ) {
      if ( exists( "cluster.net.scores" ) && i %in% colnames( cluster.net.scores ) )
        out <- cbind( out, weighted.mean( cluster.net.scores[ ,i ], sapply( cs, "[[", "nrows" ), na.rm=T ) )
      else out <- cbind( out, rep( NA, nrow( out ) ) )
      names( out )[ ncol( out ) ] <- paste( "net", i, sep="." )
    }
    if ( exists( "cluster.net.scores" ) && "net.scores" %in% colnames( cluster.net.scores ) )
      out[ ,"net.scores" ] <- weighted.mean( cluster.net.scores[ ,"net.scores" ],
                                            sapply( cs, "[[", "nrows" ), na.rm=T )
  }
  ## set.scores weighted by n.rows of bicluster
  if ( length( enrichment.sets ) > 1 ) {
    for ( i in names( set.weights ) ) {
      if ( exists( "cluster.set.scores" ) && i %in% colnames( cluster.set.scores ) ) {
        out <- cbind( out, weighted.mean( cluster.set.scores[ ,i ], sapply( cs, "[[", "nrows" ), na.rm=T ) )
      } else {
        out <- cbind( out, rep( NA, nrow( out ) ) )
      }
      names( out )[ ncol( out ) ] <- paste( "set", i, sep="." )
    }
    if ( exists( "cluster.set.scores" ) && "set.scores" %in% colnames( cluster.set.scores ) ) {
      out[ ,"set.scores" ] <- weighted.mean( cluster.set.scores[ ,"set.scores" ], sapply( cs, "[[", "nrows" ), na.rm=T )
    }
  }
  out
}

## Only override global parameters if they were not already set previously ----
##   but even if they were, copy them into the 'cmonkey.params' environment
set.param <- function( name, val, env=cmonkey.params, override=F, quiet=F ) {
  if ( ! exists( name, envir=env ) || override ) {
    if ( ! quiet ) try( { cat( name, "-> " ); str( val, digits.d=9, no.list=T ) } )
    assign( name, val, envir=env )
  } else {
    val <- get( name, envir=env )
    if ( ! quiet ) try( { cat( name, "= " ); str( val, digits.d=9, no.list=T ) } ) 
    assign( name, val, envir=env )
  }
  assign( name, val, envir=parent.frame() )
}

## Make sure all processes are killed via kill(children(),SIGKILL) ??
## Use "parallel.cores" param to determine if parallelization is desired - if it is FALSE or 0 or 1 or NA, then no.
## If it is >1 or TRUE then yes. If TRUE, then determine # of cores via multicore:::detectCores()
## This should allow running on Windows boxes when there is no multicore package.
## Note that with the use of "foreach", this can now run on Windows if we use doSMP instead of doMC !!!
##    Set the "foreach.register.backend" function to a different backend if desired.

##is.windows.system <- function() FALSE

foreach.register.backend <- function( par ) { ##, force.snow=is.windows.system() ) {
  if ( ! require( foreach ) ) return( NULL )
  if ( par > 1 && require( doMC, quietly=T ) ) registerDoMC( cores=par )
  ##if ( par == 1 ) { registerDoSEQ(); return( NULL ) }
  ##if ( ! is.windows.system() && ! force.snow ) {
  ## else if require( doSNOW ) {
  ##  if ( ! "data" %in% ls( pos=foreach:::.foreachGlobals ) ) {
  ##    cl <- makeCluster( rep( 'localhost', par ), "SOCK" )
  ##    registerDoSNOW( cl )
  ##  }
  ##}
  else registerDoSEQ()
}

get.parallel <- function( X=k.clust, verbose=F, para.cores=get( "parallel.cores" ) ) {
  if ( is.na( para.cores ) || ( is.logical( para.cores ) && para.cores == FALSE ) ||
      ( is.numeric( para.cores ) && para.cores <= 1 ) ) {
    out <- list( mc=FALSE, par=para.cores, apply=lapply )
    if ( verbose ) cat( "NOT PARALLELIZING\n" )
  } else {
    try( has.multi <- require( multicore, quietly=T ), silent=T )
    if ( ! has.multi || ( has.multi && multicore:::isChild() ) ) {    
      out <- list( mc=FALSE, par=para.cores, apply=lapply )
      if ( verbose ) cat( "NOT PARALLELIZING\n" )
    } else {
      mc <- has.multi && ! multicore:::isChild() && X > 1 && ! is.na( para.cores ) &&
      ( is.numeric( para.cores ) && para.cores > 1 ) ||
      ( is.logical( para.cores ) && para.cores == TRUE )
      par <- para.cores
      out.apply <- lapply 
      if ( mc ) {
        if ( is.logical( par ) && par == TRUE ) par <- multicore:::detectCores() ## all.tests=TRUE )
        par <- min( c( X, par, multicore:::detectCores() ) ) ## all.tests=TRUE ) ) )
        if ( verbose ) cat( "PARALLELIZING:", par, ": " )
        ## if ( ! exists( "foreach.register.backend" ) || is.null( foreach.register.backend ) ||
        ##     is.null( foreach.register.backend( par ) ) ) { ##use.foreach ) {
        ##   out.apply <- mclapply
        ## } else {
        foreach.register.backend( par )
        if ( verbose ) cat( getDoParName(), getDoParWorkers(), "\n" )
        out.apply <- function( list, FUN, ... ) foreach( l=list ) %dopar% { FUN( l, ... ) }
        ## }
        ##if ( attr( ratios, "nrow" ) > big.run ) print( gc() ) ## gc() before we spawn new copied processes
      } else {
        par <- 1
        if ( verbose ) cat( "NOT PARALLELIZING:", par, "\n" )
      }
      out <- list( mc=mc, par=par, apply=out.apply )
    }
  }
  ##if ( attr( ratios, "nrow" ) > big.run ) print( gc() ) ## gc() before we spawn new copied processes
  if ( is.numeric( out$par ) && ! is.na( out$par ) ) options( cores=out$par ) ## getOption("cores") is the default of how mclapply gets its mc.cores number
  else if ( is.na( out$par ) || ( is.logical( out$par ) && out$par == TRUE ) ) options( cores=NULL )
  else options( cores=1 )
  out
}

## Quick and dirty quantile normalization of scores matrices against each other!!!
## Optionally weight the scores matrices!!!
quantile.normalize.scores <- function( scores, weights=NULL, keep.nas=F ) {
  ## scores is a list of matrices of same dimension
  if ( ! is.list( scores ) || sum( ! sapply( scores, is.null ) ) <= 1 ) return( scores )
  scores <- scores[ ! sapply( scores, is.null ) ]
  ##ns <- names( scores )
  d <- dim( scores[[ 1 ]] )
  dn <- dimnames( scores[[ 1 ]] )
  ##scores <- sapply( scores, function( s ) as.vector( s[,] ) )
  ##tmp <- apply( scores, 2, rank, ties="min", na="keep" )
  ##tmp2 <- apply( scores, 2, sort, na.last=T ) ##sapply( 1:ncol( scores ), function( i ) scores[ tmp[ ,i ], i ] )
  tmp2 <- sapply( scores, function( i ) sort( i[,], na.last=T ) ) ##sapply( 1:ncol( scores ), function( i ) scores[ tmp[ ,i ], i ] )
  if ( is.null( weights ) ) tmp2.mn <- rowMeans( tmp2, na.rm=T )
  else { ##tmp2.mn <- apply( tmp2, 1, weighted.mean, weights, na.rm=T )
    for ( i in 1:ncol( tmp2 ) ) tmp2[ ,i ] <- tmp2[ ,i ] * weights[ i ]
    tmp2.mn <- rowMeans( tmp2, na.rm=T ) / sum( weights, na.rm=T )
  }
  rm( tmp2 )
  ##out <- lapply( 1:ncol( scores ), function( i ) tmp2.mn[ tmp[ ,i ] ] )
  ##out <- lapply( out, function( i ) { z <- matrix( i, nrow=d[ 1 ], ncol=d[ 2 ] ); dimnames( z ) <- dn; z } ) 
  ##names( out ) <- ns
  ##out
  ##tmp <- sapply( scores, function( i ) rank( i[,], ties="min", na="keep" ) )
  out <- list()
  for ( n in names( scores ) ) {
    tmp <- rank( scores[[ n ]][,], ties="min", na="keep" )
    ##scores[[ n ]][,] <-
    z <- matrix( tmp2.mn[ tmp ], nrow=d[ 1 ], ncol=d[ 2 ] ) ##tmp[ ,n ] ]
    dimnames( z ) <- dn
    if ( keep.nas ) z[ is.na( scores[[ n ]][,] ) ] <- NA
    out[[ n ]] <- z
  }
  out ##scores
}

get.all.scores <- function( ks=1:k.clust, force.row=F, force.col=F, force.motif=F, force.net=F, force.set=F,
                           quantile.normalize=T ) {
  mc <- get.parallel( length( ks ) )

  ## Compute row.scores (microarray data)
  if ( force.row || ( row.scaling[ iter ] > 0 && ! is.na( row.iters[ 1 ] ) && iter %in% row.iters ) ) {
    ##if ( ! exists( "row.scores" ) || is.null( row.scores ) || nrow( row.scores ) != attr( ratios, "nrow" ) ||
    ##    ncol( row.scores ) != max( ks ) ) {
    if ( is.null( row.scores ) ) {
      row.scores <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=max( ks ) )
      rownames( row.scores ) <- attr( ratios, "rnames" )
      row.scores <- matrix.reference( row.scores )
    } else row.scores[ ,ks ] <- 0
    ## TODO: Try to quantile normalize the ratios scores too (as is done w/ mot.scores, net.scores below) but
    ##    this is harder because each different "tmp.row" matrix may have different nrows.
    ##rs.func <- function() { ## for profiling
    for ( i in names( ratios ) ) { 
      if ( row.weights[ i ] == 0 || is.na( row.weights[ i ] ) ) next
      tmp.row <- do.call( cbind, mc$apply( ks, get.row.scores, ratios=ratios[[ i ]]
#ifndef PACKAGE
                                          , method=row.score.func
#endif
                                          ) )
      tmp <- is.infinite( tmp.row ) | is.na( tmp.row )
      if ( any( tmp ) ) tmp.row[ tmp ] <-
        quantile( tmp.row[ row.memb[,][ rownames( tmp.row ), ] == 1 & ! tmp ], 0.95 ) ##0 ## No measurement is NOT an NA - it is an observation
      tmp <- rownames( row.scores )[ rownames( row.scores ) %in% rownames( tmp.row ) ]
      row.scores[ tmp, ks ] <- row.scores[ tmp, ks ] + tmp.row[ tmp, ] * row.weights[ i ]
      rm( tmp.row, tmp )
    }
    attr( row.scores, "changed" ) <- TRUE
    ##row.scores }; row.scores <- rs.func()    
  }
  
  ## Compute col.scores (microarray data)
  ##if ( ! exists( "col.scores" ) || is.null( col.scores ) || nrow( col.scores ) != attr( ratios, "ncol" ) ||
  ##    ncol( col.scores ) != max( ks ) ) {
  if ( n.clust.per.col < k.clust && ## If using all columns, can skip this step...  
      ( force.col || ( row.scaling[ iter ] > 0 && ! is.na( col.iters[ 1 ] ) && iter %in% col.iters ) ) ) {
    if ( is.null( col.scores ) ) {
      col.scores <- matrix( 0, nrow=attr( ratios, "ncol" ), ncol=max( ks ) )
      rownames( col.scores ) <- attr( ratios, "cnames" )
      col.scores <- matrix.reference( col.scores )
    } else col.scores[ ,ks ] <- 0
    ##col.scores <- matrix( 0, nrow=attr( ratios, "ncol" ), ncol=max( ks ) )
    ## cs.func <- function() { ## for profiling
    for ( i in names( row.weights ) ) { 
      if ( row.weights[ i ] == 0 || is.na( row.weights[ i ] ) ) next
      tmp.col <- do.call( cbind, mc$apply( ks, get.col.scores, ratios=ratios[[ i ]]
#ifndef PACKAGE
                                          , method=col.score.func
#endif
                                          ) )
      tmp <- is.infinite( tmp.col ) | is.na( tmp.col )
      if ( any( tmp ) ) tmp.col[ tmp ] <-
        quantile( tmp.col[ col.memb[,][ rownames( tmp.col ), ] == 1 & ! tmp ], 0.95 ) ##0 ## No measurement is NOT an NA - it is an observation
      tmp <- rownames( col.scores )[ rownames( col.scores ) %in% rownames( tmp.col ) ]
      col.scores[ tmp, ks ] <- col.scores[ tmp, ks ] + tmp.col[ tmp, ] * row.weights[ i ]
      rm( tmp.col, tmp )
    }
    attr( col.scores, "changed" ) <- TRUE
    ## col.scores }; col.scores <- cs.func()
  }
  
  ## Run meme on each cluster (every meme.iters iterations)
  for ( i in names( mot.weights ) ) {
    if ( force.motif == "run.meme" || ( mot.scaling[ iter ] > 0 && ! is.na( meme.iters[[ i ]][ 1 ] ) &&
           iter %in% meme.iters[[ i ]] && exists( "genome.info" ) && ! no.genome.info ) ) {
      if ( mot.weights[ i ] == 0 || is.na( mot.weights[ i ] ) ) next
      tmp <- motif.all.clusters( ks, seq.type=i, verbose=T ) ##strsplit( i, " " )[[ 1 ]][ 1 ],
                                               ##algo=strsplit( i, " " )[[ 1 ]][ 2 ] )
#ifndef PACKAGE
      tmp <- list.reference( tmp, file=sprintf( "%s/meme.scores.%s", cmonkey.filename, i ), type="RDS" )
#endif
      meme.scores[[ i ]] <- tmp
    }
  }
  
  ## Compute mot.scores from meme output
  if ( force.motif == TRUE || force.motif == "run.meme" || ( mot.scaling[ iter ] > 0 && ! is.na( mot.iters[ 1 ] ) &&
                         ##iter %in% c( meme.iters, mot.iters ) && exists( "genome.info" ) && ! no.genome.info ) ) {
                                iter %in% mot.iters && exists( "genome.info" ) && ! no.genome.info ) ) {
    if ( is.null( mot.scores ) ) {
      mot.scores <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=max( ks ) ) 
      rownames( mot.scores ) <- attr( ratios, "rnames" )
      mot.scores <- matrix.reference( mot.scores )
    } else mot.scores[ ,ks ] <- 0
    tmp.mots <- list()
    for ( i in names( mot.weights ) ) {
      if ( mot.weights[ i ] == 0 || is.na( mot.weights[ i ] ) ) next
      tmp.mot <- do.call( cbind, mc$apply( ks, get.motif.scores, seq.type=i ) )
      tmp.mot[ is.infinite( tmp.mot ) | is.na( tmp.mot ) ] <- 0 ## No measurement is NOT an NA - it is an observation
      if ( quantile.normalize && sum( mot.weights > 0 & ! is.na( mot.weights ) ) > 1 ) tmp.mots[[ i ]] <- tmp.mot
      else mot.scores[ ,ks ] <- mot.scores[ ,ks ] + tmp.mot[,] * mot.weights[ i ]
      rm( tmp.mot )
    }
    if ( quantile.normalize && length( tmp.mots ) > 1 ) {
      tmp.mots <- quantile.normalize.scores( tmp.mots, weights=mot.weights[ mot.weights != 0 ] )
      for ( i in names( tmp.mots ) ) mot.scores[ ,ks ] <- mot.scores[ ,ks ] + tmp.mots[[ i ]][,] * mot.weights[ i ]
      rm( tmp.mots )
    }
    attr( mot.scores, "changed" ) <- TRUE
  }

  cluster.ns <- NULL
  ## Compute net.scores from STRING ... add weighted scores for other networks (if exist)
  if ( force.net || ( net.scaling[ iter ] > 0 && ! is.na( net.iters[ 1 ] ) && exists( "genome.info" ) &&
                     iter %in% net.iters ) ) {
    ##net.scores <- NULL
    if ( is.null( net.scores ) ) {
      net.scores <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=max( ks ) ) 
      rownames( net.scores ) <- attr( ratios, "rnames" )
      net.scores <- matrix.reference( net.scores )
    } else net.scores[ ,ks ] <- 0
    tmp.nets <- list()
    ## ns.func <- function() { ## for profiling
    for ( i in names( networks ) ) { 
      if ( net.weights[ i ] == 0 || is.na( net.weights[ i ] ) ) next
      if ( nrow( subset( networks[[ i ]], protein1 %in% attr( ratios, "rnames" ) & protein2 %in%
                        attr( ratios, "rnames" ) ) ) <= 0 ) next
      tmp.net <- do.call( cbind, mc$apply( ks, get.network.scores, net=networks[[ i ]] ) )
      if ( all( is.na( tmp.net ) ) || all( is.character( tmp.net ) ) ) next
      ##tmp.net[,] <- tmp.net[,] - max( tmp.net[ ! is.infinite( tmp.net ) ], na.rm=T ) -
      ##  abs( diff( range( tmp.net[ ! is.infinite( tmp.net ) ], na.rm=T ) ) ) / 10
      tmp.net[ is.infinite( tmp.net ) | is.na( tmp.net ) ] <- 0 ## No edge is NOT an NA - it is an observation
      if ( quantile.normalize && sum( net.weights > 0 & ! is.na( net.weights ) ) > 1 ) tmp.nets[[ i ]] <- tmp.net
      else net.scores[ ,ks ] <- net.scores[ ,ks ] + tmp.net[,] * net.weights[ i ]
      cluster.ns <- cbind( cluster.ns, do.call( c, mc$apply( ks, function( k ) mean( tmp.net[ get.rows( k ), k ],
                                                                                    na.rm=T, trim=0.05 ) ) ) )
      colnames( cluster.ns )[ ncol( cluster.ns ) ] <- i
      rm( tmp.net )
    }
    
    if ( quantile.normalize && length( tmp.nets ) > 1 ) {
      tmp.nets <- quantile.normalize.scores( tmp.nets, weights=net.weights[ net.weights != 0 ] )
      for ( i in names( tmp.nets ) ) net.scores[ ,ks ] <- net.scores[ ,ks ] + tmp.nets[[ i ]][,] * net.weights[ i ]
      rm( tmp.nets )
    }
    attr( net.scores, "changed" ) <- TRUE

    cluster.ns <- cbind( cluster.ns, do.call( c, mc$apply( ks, function( k ) mean( net.scores[ get.rows( k ), k ],
                                                                                  na.rm=T, trim=0.05 ) ) ) )
    colnames( cluster.ns )[ ncol( cluster.ns ) ] <- "net.scores"
    ##list( net.scores, cluster.ns ) }; tmp <- ns.func(); net.scores <- tmp[[ 1 ]]; cluster.ns <- tmp[[ 2 ]]
  }
  
  ## Compute set.scores for an iteration
  if ( force.set || ( set.scaling[ iter ] > 0 && ! is.na( set.iters[ 1 ] ) && exists( "genome.info" ) && iter %in% set.iters ) ) {
    if ( is.null( set.scores ) ) {
      set.scores <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=max( ks ) ) 
      rownames( set.scores ) <- attr( ratios, "rnames" )
      set.scores <- matrix.reference( set.scores )
    } else {
      set.scores[ ,ks ] <- 0
    }
    tmp.sets <- list()
    for ( set.type in names( enrichment.sets ) ) { 
      ## If the weights for the set.type are 0 then skip this set type
      if ( set.weights[ set.type ] == 0 || is.na( set.weights[ set.type ] ) ) {
        next
      }
      ## Run the analysis
      tmp.set <- do.call( cbind, mc$apply( ks, get.set.enrichment.scores, set=set.type ) )
      tmp.set[ is.infinite( tmp.set ) | is.na( tmp.set ) ] <- 0 ## Ensures that no edge is an NA
      # If there is an issue with the weights then don't use them
      if ( quantile.normalize && sum( set.weights > 0 & ! is.na( set.weights ) ) > 1 ) {
        tmp.sets[[ set.type ]] <- tmp.set
      # Or do if they are working
      } else {
        set.scores[ ,ks ] <- set.scores[ ,ks ] + tmp.set[,] * set.weights[ set.type ]
      }
      cluster.ses <- cbind( cluster.ses, do.call( c, mc$apply( ks, function( k ) mean( tmp.set[ get.rows( k ), k ], na.rm=T, trim=0.05 ) ) ) )
      colnames( cluster.ses )[ ncol( cluster.ses ) ] <- set.type
      rm( tmp.set )
    }
    
    if ( quantile.normalize && length( tmp.sets ) > 1 ) {
      tmp.sets <- quantile.normalize.scores( tmp.sets, weights=set.weights[ set.weights != 0 ] )
      for ( set.type in names( tmp.sets ) ) {
        set.scores[ ,ks ] <- set.scores[ ,ks ] + tmp.sets[[ set.type ]][,] * set.weights[ set.type ]
      }
      rm( tmp.sets )
    }
    attr( set.scores, "changed" ) <- TRUE
    cluster.ses <- cbind( cluster.ses, do.call( c, mc$apply( ks, function( k ) mean( set.scores[ get.rows( k ), k ], na.rm=T, trim=0.05 ) ) ) )
    colnames( cluster.ses )[ ncol( cluster.ses ) ] <- "set.scores"
  }
  list( r=row.scores, c=col.scores, m=mot.scores, ms=meme.scores, n=net.scores, se=set.scores, cns=cluster.ns, cses=cluster.ses )
}

get.rows <- function( k, rm=get("row.membership") ) { out <- unique( rownames( which( rm[] == k, arr=T ) ) );
                                                      if ( is.null( out ) ) out <- character(); out } 

get.cols <- function( k, cm=get("col.membership") ) { out <- unique( rownames( which( cm[] == k, arr=T ) ) );
                                                      if ( is.null( out ) ) out <- character(); out } 

#' Create old clusterStack-style object from row.membership and col.membership
#'  SD Notes: force is a little mysterious to me.
#'
#' @params ks The clusters to include (DEFAULT: 1:k.clust)
#' @params force If FALSE, will return the existing clusterStack if possible (DEFAULT: F)
#' @params ...  Must include atleast "iter"
#' @export
#' @usage clusterStack <- get.clusterStack( ks=1:k.clust, force=F, ... ) 
get.clusterStack <- function( ks=1:k.clust, force=F, ... ) {
  if ( ! force && ! is.null( attr( clusterStack, "iter" ) ) && attr( clusterStack, "iter" ) == iter )
    return( clusterStack )
  mc <- get.parallel( length( ks ) )
  clusterStack <- mc$apply( ks, get.clust, ... ) ###, mc.cores=mc$par )
#ifndef PACKAGE
  tmp <- list.reference( clusterStack, file=sprintf( "%s/clusterStack", cmonkey.filename ), type="RDS" )
  clusterStack <- tmp
#endif
  attr( clusterStack, "iter" ) <- iter
  clusterStack
}

#' Get the unpreprocessed ratios matrix (i.e. ratios.raw)
#'  
#' @export
#' @usage ratios <- get.unpreprocessed.ratios()
get.unpreprocessed.ratios <- function( ... ) {
  return( ratios.raw )
}

cluster.resid <- function( k, rats.inds="COMBINED", varNorm=F, in.cols=T, ... ) {
  ## FLOC residual number is a good statistic
  residual.submatrix <- function( rats, rows, cols, varNorm=F, ... ) {
    rows <- rows[ rows %in% rownames( rats ) ]
    cols <- cols[ cols %in% colnames( rats ) ]
    if ( length( rows ) <= 1 || length( cols ) <= 1 ) return( 1 )
    maxRowVar <- attr( rats, "maxRowVar" )
    rats <- rats[ rows, cols ]
    if ( is.vector( rats ) || any( dim( rats ) <= 1 ) || mean( is.na( rats ) ) > 0.95 ) return( 1 )

    d.rows <- rowMeans( rats, na.rm=T )
    d.cols <- colMeans( rats, na.rm=T )
    d.all <- mean( d.rows, na.rm=T )

    rij <- rats + d.all
    rij[,] <- rij[,] - matrix( d.cols, nrow=nrow( rij ), ncol=ncol( rij ), byrow=T )
    rij[,] <- rij[,] - matrix( d.rows, nrow=nrow( rij ), ncol=ncol( rij ), byrow=F )
    ##rij[,] <- rij[,] - outer( d.rows, d.cols, '+' ) ## more elegant but slower!

    average.r <- mean( abs( rij ), na.rm = TRUE )
    if ( varNorm && ! is.null( maxRowVar ) ) {
      ##maxRowVar <- attr( rats, "maxRowVar" )
      row.var <- mean( apply( rats, 1, var, use = "pairwise.complete.obs" ), na.rm=T )
      if ( is.na( row.var ) || row.var > maxRowVar ) row.var <- maxRowVar
      average.r <- average.r / row.var
    }
    average.r
  }
  
  inds <- rats.inds
  if ( rats.inds[ 1 ] == "COMBINED" ) inds <- names( get( "row.weights" ) )
  resids <- sapply( ratios[ inds ], function( rn ) {
#ifndef PACKAGE
    if ( row.score.func == "orig" ) { ## Original FLOC resid score works for co-expressed genes
#endif
      if ( in.cols ) residual.submatrix( rn, get.rows( k ), get.cols( k ), varNorm=varNorm )
      else residual.submatrix( rn, get.rows( k ), colnames( rn )[ ! colnames( rn ) %in% get.cols( k ) ],
                              varNorm=varNorm )
#ifndef PACKAGE
    }
    else { ## using cor2 or mi will results in pos- and neg- correlated genes, FLOC score wont work!
      if ( in.cols ) mean( get.row.scores( k, for.rows=get.rows( k ), ratios=rn, method=row.score.func ) )
      else mean( get.row.scores( k, cols=get.cols( k ), for.rows=get.rows( k ), ratios=rn, method=row.score.func ) )
    }
#endif
  } )
  if ( rats.inds[ 1 ] == "COMBINED" ) resids <- weighted.mean( resids, row.weights[ inds ], na.rm=T )
  if ( rats.inds[ 1 ] != "COMBINED" && length( resids ) < length( inds ) && all( is.na( resids ) ) ) {
    resids <- rep( NA, length( inds ) ); names( resids ) <- inds }
  resids
}

cluster.pclust <- function( k, mot.inds="COMBINED" ) { ## actually returns a list with e-vals and p-vals
  inds <- mot.inds
  if ( is.null( inds ) ) return( list( p.clusts=NA, e.vals=NA ) )
  if ( mot.inds[ 1 ] == "COMBINED" ) inds <- names( get( "mot.weights" ) )
  if ( is.null( inds ) ) return( list( p.clusts=NA, e.vals=NA ) )
  rows <- get.rows( k )
  p.clusts <- sapply( inds, function( n ) {
    ##ms <- meme.scores[[ n ]][[ k ]]
    ##if ( length( rows ) > 0 && ! is.null( ms$pv.ev ) && ! is.null( ms$pv.ev[[ 1 ]] ) )
    ##  mean( log10( ms$pv.ev[[ 1 ]][ rownames( ms$pv.ev[[ 1 ]] ) %in% rows, "p.value" ] ), na.rm=T ) else NA
    pvs <- meme.scores[[ n ]]$all.pv
    if ( is.null( pvs ) ) return( NA )
    rows <- rows[ rows %in% rownames( pvs ) ]
    if ( length( rows ) > 0 ) mean( log10( pvs[ rows, k ] ), na.rm=T ) else NA
  } )
  
  e.vals <- sapply( inds, function( n ) {
    ms <- meme.scores[[ n ]][[ k ]]
    sapply( 1:length( ms$meme.out ),
           function( i ) if ( length( rows ) > 0 && ! is.null( ms$meme.out ) && ! is.null( ms$meme.out[[ i ]] ) )
           ms$meme.out[[ i ]]$e.value else NA )
  } )
  if ( ! is.matrix( e.vals ) ) e.vals <- t( t( e.vals ) )

  if ( mot.inds[ 1 ] == "COMBINED" ) {
    p.clusts <- weighted.mean( p.clusts, mot.weights[ inds ], na.rm=T )
    e.vals <- apply( e.vals, 2, weighted.mean, mot.weights[ inds ], na.rm=T )
  } else if ( mot.inds[ 1 ] != "COMBINED" ) {
    if ( length( p.clusts ) < length( inds ) && all( is.na( p.clusts ) ) ) p.clusts <- rep( NA, length( inds ) )
    e.vals <- apply( e.vals, 1, function( i ) if ( length( i ) < length( inds ) && all( is.na( i ) ) )
                    rep( NA, length( inds ) ) else i )
  }
  if ( ! is.matrix( e.vals ) ) e.vals <- t( t( e.vals ) )
  if ( is.matrix( e.vals ) && ncol( e.vals ) != length( inds ) ) e.vals <- t( e.vals )
  if ( mot.inds[ 1 ] != "COMBINED" ) names( p.clusts ) <- colnames( e.vals ) <- inds
  else e.vals <- as.vector( e.vals )
  list( p.clusts=p.clusts, e.vals=e.vals )
}

## Create old clusterStack-style cluster object from row.membership and col.membership
get.clust <- function( k, fill=T, fill.motif=T, seq.type=names( mot.weights ), varNorm=F, ... ) { ##"upstream" ) {
  gen.clust <- function( rowNames, colNames=NA, fill=F, motif=F, n.motifs=3, ... ) {
    rowNames <- rowNames[ rowNames %in% attr( ratios, "rnames" ) ]
    if ( ! is.null( colNames ) && length( colNames ) > 1 && ! is.na( colNames ) )
      colNames <- colNames[ colNames %in% attr( ratios, "cnames" ) ]
    c.tmp <- list( nrows=length( rowNames ), ncols=length( colNames ), rows=rowNames, cols=colNames,
                  k=999, p.clust=1.0, e.val=rep( 999, n.motifs ), resid={out=rep(NA,length(row.weights));names(out)<-names(row.weights);resid=out} )
    if ( fill && c.tmp$nrows > 0 && c.tmp$ncols > 0 && ! all( is.na( colNames ) ) )
      c.tmp$resid <- cluster.resid( k, names( row.weights ), varNorm=varNorm, ... ) 
    return( c.tmp )
  }
  
  cols <- get.cols( k ); rows <- get.rows( k )
  if ( length( cols ) <= 0 ) cols <- NA
  clust <- gen.clust( rows, cols, fill=fill, motif=F, n.motifs=max( unlist( n.motifs ) ) ) 
  clust$k <- k

  if ( fill.motif ) {
    tmp <- cluster.pclust( k, seq.type )
    clust$e.val <- tmp$e.vals
    clust$p.clust <- tmp$p.clusts
  }  
  clust
}

## Allow weighting of columns (used in get.row.scores) (this can be "overridden")
get.col.weights <- function( rows, cols, ratios ) {
  NA ## Default: unweighted
##   if ( length( rows ) == 1 ) return( NA ) ## Below is example code weighting columns by their "index of dispersion" over the rows in the cluster
##   rats <- ratios[ rows, cols, drop=F ]
##   vars <- apply( rats, 2, var, na.rm=T )
##   vars <- vars / ( abs( apply( rats, 2, mean, na.rm=T ) ) + 0.05 )
##   vars[ is.na( vars ) | vars == 0 ] <- 1
##   weights <- 1 / vars
##   weights <- weights / sum( weights, na.rm=T ) * length( weights )
##   names( weights ) <- cols
##   weights
}

## Allow weighting of rows (used in get.col.scores) (this can be "overridden")
get.row.weights <- function( rows, cols, ratios ) NA

get.row.scores <- function( k, cols=get.cols( k ), for.rows="all", ratios=ratios[[ 1 ]],
#ifndef PACKAGE
                           method=c("cor2","dist","orig")[1],
#endif
                           ... ) {
  ## Compute scores for ALL rows (over just the cols IN the cluster)
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k ) 
  else rows <- k
  if ( for.rows[ 1 ] == "all" ) for.rows <- rownames( ratios ) ##attr( ratios, "rnames" )
  rows <- rows[ rows %in% rownames( ratios ) ]
  cols <- cols[ cols %in% colnames( ratios ) ]
  if ( length( rows ) < 1 || length( cols ) <= 1 ) return( rep( NA, length( for.rows ) ) )

#ifndef PACKAGE
  if ( method == "orig" ) { ## ORIGINAL, GOOD, FAST SCORE FUNCTION
#endif
    rats <- ratios[ for.rows, cols, drop=F ]
    rats.mn <- colMeans( rats[ rows, , drop=F ], na.rm=T )
    rats.mn <- matrix( rats.mn, nrow=nrow( rats ), ncol=ncol( rats ), byrow=T )
    rats[,] <- ( rats[,] - rats.mn )^2 ## abs(
    col.weights <- if ( exists( "get.col.weights" ) ) get.col.weights( rows, cols, ratios ) else NA
    if ( is.na( col.weights[ 1 ] ) ) rats <- rowMeans( rats, na.rm=T )
    else rats <- apply( abs( rats ), 1, weighted.mean, w=col.weights[ cols ], na.rm=T )
    rats <- log( rats + 1e-99 )
#ifndef PACKAGE
  }
  else if ( method == "pval" ) rats <- get.row.scores.pVals( k, cols, rows, ratios, method, ... ) ## Sam's new func.
  else if ( exists( "get.row.scores.NEW" ) ) rats <- get.row.scores.NEW( k, cols, rows, ratios, method, ... )
#endif
  return( rats )
}

#ifndef PACKAGE
## See minet and infotheo packages for information-based distances and implementations like ARACNE
## TODO: how do we include row.weights with the new scores?
get.row.scores.NEW <- function( k, cols=get.cols( k ), rows, ratios=ratios[[ 1 ]],
                           method=c("cor2","abscor","cor","dist")[1], ... ) {
  rats <- ratios[ ,cols, drop=F ]
  if ( method == "dist" || substr( method, 1, 5 ) == "dist." ) { ## Try a more flexible dist() based function -- mean of distances (slower)
    ## if ( method == "dist" ) method <- "dist.euc" ## Default - euclidean
    ## d <- as.matrix( dist( rats, method=substring( method, 6 ) ) ); diag( d ) <- NA
    ## rats <- log( apply( d[ ,rows, drop=F ], 1, mean, na.rm=T, ... ) + 1e-99 )
    require( proxy ) ## Better, faster!!! (10x)
    if ( method == "dist" ) method <- "dist.Euclidean" ## Default - euclidean    
    d <- proxy::dist( rats, rats[ rows, ,drop=F ], method=substring( method, 6 ) ); d[ cbind( rows, rows ) ] <- NA
    rats <- log( apply( d, 1, mean, na.rm=T, ... ) + 1e-99 )
  }
  else if ( method %in% c( "cor", "cor2", "abscor" ) ) { ## Try a cor^2 (allow anticorr) function -- median of cor^2 's
    rats <- t( rats )
    ##d <- cor( rats, use="pairwise", ... ); diag( d ) <- NA
    ##rats <- apply( log( 1 - d[ ,rows, drop=F ] + 1e-99 ), 1, mean, na.rm=T, ... ) ## + 1 + 1e-99 )
    d <- t( cor( rats[ ,rows, drop=F ], rats, use="pairwise", ... ) ); d[ cbind( rows, rows ) ] <- NA ## faster! (10x)
    if ( method == "cor2" ) d <- d^2
    else if ( method == "abscor" ) d <- abs( d )
    rats <- apply( log( 1 - d + 1e-99 ), 1, mean, na.rm=T, ... ) ## + 1 + 1e-99 )
  }
  ## Try mutual information!!! Options: estimator="mi.empirical", "mi.mm", "mi.shrink", "mi.sg" need a non-default
  ##    'discretizer': disc="none", "equalfreq", "equalwidth" or "globalequalwidth" (see infotheo:::discretize)
  ##    estimator="pearson","spearman","kendall" is not discretized (but use "cor2" option instead - same but faster)
  else if ( method == "mi" || substr( method, 1, 3 ) == "mi." ) { 
    require( minet )
    if ( method == "mi" ) method <- "mi.spearman" ## The default estimator
    d <- build.mim( t( rats ), estimator=substring( method, 4 ), disc="equalwidth" ); diag( d ) <- NA ## slow!!!
    rats <- -apply( d[ ,rows, drop=F ], 1, mean, na.rm=T, ... )
  }
  return( rats )
}
#endif

## get.row.scores.cpp <- function( rats, rats.mn ) { ##, rank=rank.ties ) {
##   if ( TRUE && ! is.loaded( "get_row_scores" ) ) {
##     so.file <- "cpp/get_row_scores"
##     if ( ! exists( "is.slave" ) || ! is.slave ) {
##       require( Rcpp ) ## Dont really need to load the package - but need to make sure it's installed
##       code <-
##         "RcppExport SEXP get_row_scores( SEXP rats, SEXP rats_mn ) {
##        SEXP rl = R_NilValue; // Use this when there is nothing to be returned.
##        char* exceptionMesg = NULL; 
##        try {
##          RcppMatrixView<double> r( rats );
##          RcppVectorView<double> rm( rats_mn );
##          int nr = r.dim1(), nc = r.dim2();
##          vector<double> out( nr );
##          vector<int> ni( nr ); //, nj( nc );

##          for ( int j = 0; j < nc; j ++ ) {
##            double rr = rm( j );
##            if ( ISNA( rr ) ) continue;
##            //nj[ j ] = nj[ j ] + 1;
##            for ( int i = 0; i < nr; i ++ ) {
##              double rrr = r( i, j );
##              if ( ISNA( rrr ) ) continue;
##              out[ i ] = out[ i ] + fabs( rrr - rr ); //pow( rrr - rr, 2 );
##              ni[ i ] = ni[ i ] + 1;
##            }
##          }
##          for ( int i = 0; i < nr; i ++ ) out[ i ] = out[ i ] / ni[ i ];

##          RcppResultSet rs; // Build result set to be returned as a list to R.
##          rs.add( \"out\", out );
##          rl = rs.getReturnList(); // Get the list to be returned to R.
##        } catch( std::exception& ex ) { 
##          exceptionMesg = copyMessageToR( ex.what() ); 
##        } catch( ... ) { 
##          exceptionMesg = copyMessageToR( \"unknown reason\" ); 
##        } 
##        if ( exceptionMesg != NULL ) error( exceptionMesg ); 
##        return rl;
##      }"

##       source( "Rcpp_compile3.R", local=T ) ## Loads compile.cpp.code() (requires Rcpp package)
##       compile.rcpp.code( code, fname=so.file, verbose=T )
##     } else {
##       dyn.load( paste( so.file, ".so", sep="" ) )
##     }
##   }
##   out <- .Call( "get_row_scores", rats, rats.mn )$out
##   out[ out == 0 ] <- NA
##   out
## }

## Default is to normalize (diff)^2 by mean expression level, similar to "index of dispersion"
##    http://en.wikipedia.org/wiki/Index_of_dispersion
get.col.scores <- function( k, for.cols="all", ratios=ratios[[ 1 ]], 
#ifndef PACKAGE
                           method=c("new","orig","ent")[1],                            
#endif
                           norm.method=c("mean","all.colVars","none")[1], ... ) {
  ## Compute scores for ALL cols (over just the rows IN the cluster)
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k ) 
  else rows <- k
  if ( for.cols[ 1 ] == "all" ) for.cols <- colnames( ratios ) ##attr( ratios, "cnames" )
  rows <- rows[ rows %in% rownames( ratios ) ]
  if ( length( rows ) <= 1 ) return( rep( NA, length( for.cols ) ) )

  rats <- ratios[ rows, for.cols, drop=F ]
  row.weights <- if ( exists( "get.row.weights" ) ) get.row.weights( rows, cols, ratios ) else NA

#ifndef PACKAGE
  if ( method == "orig" ) { ## OLD, GOOD EUCLIDEAN DISTANCE METRIC
#endif
    if ( is.na( row.weights[ 1 ] ) ) { ## Default
      rats.mn <- matrix( colMeans( rats, na.rm=T ), nrow=nrow( rats ), ncol=ncol( rats ), byrow=T )
    } else { ## Custom row weights
      rats.mn <- matrix( apply( rats, 2, weighted.mean, w=row.weights[ rows ], na.rm=T ), ncol=ncol( rats ), byrow=T )
    }
    
    rats[,] <- ( rats[,] - rats.mn )^2 ## abs( Multiplying is faster than squaring
    rats <- colMeans( rats, na.rm=T )
#ifndef PACKAGE
  }

  ## else if ( method == "ent" ) {  ## ENTROPY BASED METHOD (do it on abs(rats) to allow for repression too
  ##   require( infotheo )
  ##   rats <- apply( discretize( rats, disc="globalequalwidth", nbins=ncol( ratios ) / 2 ), 2, entropy )
  ## }
  
  else if ( method == "new" ) { ## TRY -- score columns on amount by which they contribute to avg in-row score,
    ###  however that score is calculated! Note - this shouldn't have to worry about the varNorm anymore (???)
    mn <- mean( get.row.scores( k, ratios=ratios, method=row.score.func, ... ), na.rm=T )
    rows <- get.rows( k )
    cols <- get.cols( k )
    rats <- -sapply( for.cols, function( cc ) {
      ##if(row.score.func!="orig")cat(k,cc,"\n")
      if ( is.na( row.weights[ 1 ] ) ) {
        if ( cc %in% cols ) mn / mean( get.row.scores( rows, cols=cols[ cols != cc ], ratios=ratios,
                                                      method=row.score.func, ... ), na.rm=T ) ## - mn
        ##else mean( mn -
        else mean( get.row.scores( k, cols=c( cols, cc ), ratios=ratios, method=row.score.func, ... ), na.rm=T ) / mn
      } else {
        if ( cc %in% cols ) mn / weighted.mean( get.row.scores( rows, cols=cols[ cols != cc ], ratios=ratios,
                                                  method=row.score.func, ... ), w=row.weights, na.rm=T ) ## - mn
        ##else weighted.mean( mn -
        else weighted.mean( get.row.scores( rows, cols=c( cols, cc ), ratios=ratios, method=row.score.func, ... ),
                           w=row.weights, na.rm=T ) / mn
      }
    } )
    return( rats )
  }
#endif
  
  var.norm <- 0.99
  if ( norm.method == "all.colVars" ) {
    all.colVars <- attr( ratios, "all.colVars" )
    if ( ! is.null( all.colVars ) ) var.norm <- all.colVars[ for.cols ]
  } else if ( norm.method == "mean" ) {
#ifndef PACKAGE
    if ( ! exists( "rats.mn" ) ) {
      row.weights <- if ( exists( "get.row.weights" ) ) get.row.weights( rows, cols, ratios ) else NA
      rats.tmp <- ratios[ rows, for.cols, drop=F ]
      if ( is.na( row.weights[ 1 ] ) ) { ## Default
        rats.mn <- matrix( colMeans( rats.tmp, na.rm=T ), nrow=nrow( rats.tmp ), ncol=ncol( rats.tmp ), byrow=T )
      } else { ## Custom row weights
        rats.mn <- matrix( apply( rats.tmp, 2, weighted.mean, w=row.weights[ rows ], na.rm=T ),
                          ncol=ncol( rats.tmp ), byrow=T )
      }
    }
#endif
    var.norm <- abs( rats.mn[ 1, ] ) ##0.99 ## Use the mean expr. level (higher expressed expected to have higher noise)
  }
  
  ##col.weights <- get.col.weights( rows, cols )
  ##if ( is.na( col.weights ) )
  rats <- rats / ( var.norm + 0.01 ) ## default
  ##else rats <- colMeans( rats, na.rm=T ) / ( var.norm * col.weights[ cols ] + 0.01 ) ## customized col. weights

  ##return( log( rats + 1e-99 ) )
  rats
}

get.motif.scores <- function( k, seq.type="upstream meme", for.rows="all" ) { ##m=meme.scores$upstream[[ k ]], for.rows="all" ) { 
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k ) 
  else rows <- k
  if ( for.rows[ 1 ] == "all" ) for.rows <- attr( ratios, "rnames" )
  if ( length( rows ) <= 1 || is.null( meme.scores[[ seq.type ]]$all.pv ) ) ##is.null( m ) || is.null( m$pv.ev ) )
    return( rep( NA, length( for.rows ) ) )
  ##m.scores <- log( m$pv.ev[[ 1 ]][ ,"p.value" ] ); names( m.scores ) <- rownames( m$pv.ev[[ 1 ]] )
  m.scores <- log( meme.scores[[ seq.type ]]$all.pv[ ,k ] )
  m.scores <- m.scores[ for.rows ] ##; names( m.scores ) <- for.rows
  return( m.scores )
}

get.network.scores <- function( k, net=networks$string, for.rows="all", p1.col="protein1", p2.col="protein2", 
                               score.col="combined_score" ) {
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k )
  else rows <- k
  if ( for.rows[ 1 ] == "all" ) for.rows <- attr( ratios, "rnames" )
  if ( length( rows ) < 1 ) return( rep( NA, length( for.rows ) ) )
  cons <- net[ as.character( net[[ p1.col ]] ) %in% rows, c( p2.col, score.col ), drop=F ] ## Note subset is slower!
  if ( is.null( cons ) || nrow( cons ) <= 0 ) return( rep( NA, length( for.rows ) ) )
  cons <- cons[ as.character( cons[[ p2.col ]] ) %in% for.rows, , drop=F ]
  ##cons <- get.cluster.network( rows, net, for.rows )
  if ( is.null( cons ) || nrow( cons ) <= 0 ) return( rep( NA, length( for.rows ) ) )
  tmp <- tapply( as.numeric( cons[[ score.col ]] ), as.character( cons[[ p2.col ]] ), sum, na.rm=T ) / length( rows )
  scores <- rep( NA, length( for.rows ) ); names( scores ) <- for.rows
  scores[ names( tmp ) ] <- tmp
  return( -log( scores + 1 ) )
}

pareto.adjust.weights <- function( kwin=51, diter=21, max.delta=0.05 ) {
  out.scaling <- c( row=row.scaling[ iter - 1 ], mot=mot.scaling[ iter - 1 ], net=net.scaling[ iter - 1 ], set=set.scaling[ iter - 1 ] )
  orig.scaling <- c( row=row.scaling[ iter ], mot=mot.scaling[ iter ], net=net.scaling[ iter ], set=set.scaling[ iter ] )
  if ( iter < diter/2 ) return( out.scaling )
  diter <- min( diter, iter )
  kwin <- min( kwin, iter )
  avs <- NULL
  for ( i in names( out.scaling ) ) {##c( "row", "mot", "net" ) ) {
    if ( i == "row" ) col <- "row.scores" ##"resid"
    else if ( i == "mot" ) col <- "mot.scores" ##"p.clust"
    else col <- paste( i, "scores", sep="." )
    tmp <- stats[ ,col ]; tmp <- tmp[ ! is.na( tmp ) ]
    if ( length( tmp ) > 0 ) {
      av <- runmed( tmp, k=kwin )
      dy <- ( av[ length( av ) ] - av[ max( 1, length( av ) - diter ) ] ) / diff( range( av, na.rm=T ) )
      if ( ! is.na( dy ) ) out.scaling[ i ] <- out.scaling[ i ] + min( 0.1, dy )
      out.scaling[ i ] <- max( out.scaling[ i ], orig.scaling[ i ] )
    }
  }
  out.scaling
}

pareto.adjust.weights.OLD <- function( iter, delta.iter=200, delta.factor=1, n.avg=50, max.delta=0.05 ) { ##0.01
  if ( iter == 1 ) return( c( row=row.scaling[ iter ], mot=mot.scaling[ iter ], net=net.scaling[ iter ] ) )
  out.scaling <- c( row=row.scaling[ iter - 1 ], mot=mot.scaling[ iter - 1 ], net=net.scaling[ iter - 1 ] )
  if ( iter < delta.iter + n.avg + 10 ) return( out.scaling )
  all.diffs <- numeric()
  for ( i in c( "row", "mot", "net" ) ) {
    if ( i == "row" ) col <- "resid"
    else if ( i == "mot" ) col <- "p.clust"
    else col <- paste( i, "scores", sep="." )
    tops <- stats[[ col ]][ (iter-n.avg+1):iter ]
    bots <- stats[[ col ]][ (iter-n.avg+1):iter - delta.iter ]
    all.diffs[ i ] <- mean( tops - bots, na.rm=T ) / abs( diff( range( stats[[ col ]], na.rm=T ) ) )
  }
  for ( i in which( all.diffs > 0 & ! is.na( all.diffs ) & ! is.na( out.scaling ) & out.scaling > 0 ) )
    out.scaling[ names( all.diffs )[ i ] ] <- out.scaling[ names( all.diffs )[ i ] ] +
      min( c( out.scaling[ i ] * max.delta, all.diffs[ i ] * delta.factor ) )
  out.scaling
}

## WARNING: This function relies on A LOT of global variables!!!
get.combined.scores <- function( quantile.normalize=F ) {
  if ( ! exists( "r.scores" ) || is.null( r.scores ) ) {
    r.scores <- row.scores[,]
    r.scores <- matrix.reference( r.scores )
  }
  r.scores[,] <- row.scores[,]
  ##if ( attr( row.scores, "changed" ) == TRUE ) {
#ifndef PACKAGE
  if ( ! quantile.normalize ) {
    tmp <- r.scores[,] < -20; r.scores[,][ tmp ] <- min( r.scores[,][ ! tmp ], na.rm=T )
    rsm <- r.scores[,][ row.memb[,] == 1 ]
    tmp <- mad( rsm, na.rm=T )
    if ( tmp != 0 ) r.scores[,] <- ( r.scores[,] - median( rsm, na.rm=T ) ) / tmp
    else { tmp <- sd( rsm, na.rm=T ); if ( tmp != 0 ) r.scores[,] <- ( r.scores[,] - median( rsm, na.rm=T ) ) / tmp }
    rm( tmp, rsm )
  }
#endif
  
  tmp <- r.scores[,] < -20; r.scores[,][ tmp ] <- min( r.scores[,][ ! tmp ], na.rm=T ); rm( tmp ) ## -220
  r.scores[,][ is.infinite( r.scores[,] ) ] <- NA
  r.scores[,][ is.na( r.scores[,] ) ] <- max( r.scores[,], na.rm=T )
  ## }
  ##cat( "HERE: row", r.scores[1,1], "\n" )
  
#ifndef PACKAGE
  if ( ! quantile.normalize && ! is.null( mot.scores ) || ! is.null( net.scores ) )
    rs.quant <- quantile( r.scores[,], 0.01, na.rm=T )
#endif

  if ( ! is.null( mot.scores ) ) {
    if ( ! exists( "m.scores" ) || is.null( m.scores ) ) {
      m.scores <- mot.scores[,]
      m.scores <- matrix.reference( m.scores ) ## not useful as m.scores (and n.scores) are not saved in the env.
    } else m.scores[,] <- mot.scores[,]
  } else m.scores <- NULL
  if ( ! is.null( mot.scores ) && ! is.null( m.scores ) && ( is.null( attr( mot.scores, "changed" ) ) ||
      attr( mot.scores, "changed" ) == TRUE ) ) { ## Only need to do this if net scores changed!!
    ##cat("HERE: mot", mot.scores[1,1], "\n")
    tmp <- m.scores[,] < -20; m.scores[,][ tmp ] <- min( m.scores[,][ ! tmp ], na.rm=T ); rm( tmp ) ## effective zero is 10^-20
#ifndef PACKAGE
    if ( ! quantile.normalize ) {
      m.scores[,] <- m.scores[,] - quantile( m.scores[,], 0.99, na.rm=T ) ## Make it so no mot -> zero score -> no penalization
      m.scores[,] <- m.scores[,] / abs( quantile( m.scores[,], 0.01, na.rm=T ) ) * abs( rs.quant ) ## Make it so min(mot.scores) == min(row.scores) so no mot.score is too dominant
    }
#endif
  } ##else m.scores <- NULL
  
  if ( ! is.null( net.scores ) ) {
    if ( ! exists( "n.scores" ) || is.null( n.scores ) ) {
      n.scores <- net.scores[,]
      n.scores <- matrix.reference( n.scores ) ## not useful as m.scores (and n.scores) are not saved in the env.
    } else n.scores[,] <- net.scores[,]
  } else n.scores <- NULL
  if ( ! is.null( net.scores ) && ! is.null( n.scores ) && ( is.null( attr( net.scores, "changed" ) ) ||
      attr( net.scores, "changed" ) == TRUE ) ) { ## Only need to do this if net scores changed!!
    ##cat("HERE: net", n.scores[1,1],"\n")
    n.scores[,] <- n.scores[,] - quantile( n.scores[,], 0.99, na.rm=T ) ## Make it so no edge -> zero score -> no penalization
#ifndef PACKAGE
    if ( ! quantile.normalize ) {
      qqq <- abs( quantile( n.scores[,], 0.01, na.rm=T ) )
      if ( qqq == 0 ) qqq <- sort( n.scores[,] )[ 10 ] ## For really sparse networks
      if ( qqq == 0 ) qqq <- min( n.scores[,], na.rm=T ) ## For really sparse networks
      if ( qqq != 0 ) n.scores[,] <- n.scores[,] / qqq * abs( rs.quant ) ## Make it so min(net.scores) == min(row.scores) so no net.score is too dominant
      rm( qqq )
    }
#endif
  }

  if ( ! is.null( set.scores ) ) {
    if ( ! exists( "se.scores" ) || is.null( n.scores ) ) {
      se.scores <- set.scores[,]
      se.scores <- matrix.reference( se.scores )
    } else {
      se.scores[,] <- set.scores[,]
    }
  } else {
    se.scores <- NULL
  }
  if ( ! is.null( set.scores ) && ! is.null( se.scores ) && ( is.null( attr( set.scores, "changed" ) ) || attr( set.scores, "changed" ) == TRUE ) ) {
    se.scores[,] <- se.scores[,] - quantile( se.scores[,], 0.99, na.rm=T ) ## Make it so no edge -> zero score -> no penalization
    if ( ! quantile.normalize ) {
      qqq <- abs( quantile( se.scores[,], 0.01, na.rm=T ) )
      if ( qqq == 0 ) {
        qqq <- sort( se.scores[,] )[ 10 ] ## For really sparse networks
      }
      if ( qqq == 0 ) {
        qqq <- min( se.scores[,], na.rm=T ) ## For really sparse networks
      }
      if ( qqq != 0 ) {
        se.scores[,] <- se.scores[,] / qqq * abs( rs.quant ) ## Make it so min(net.scores) == min(row.scores) so no net.score is too dominant
      }
      rm( qqq )
    }
  }

  if ( ! is.null( col.scores ) ) {
    ## c.scores <- NULL
    if ( ! exists( "c.scores" ) || is.null( c.scores ) ) {
      c.scores <- col.scores[,] * 0
      c.scores <- matrix.reference( c.scores )
    } else c.scores[,] <- col.scores[,]
    if ( is.null( attr( col.scores, "changed" ) ) || attr( col.scores, "changed" ) == TRUE ) {
      ##cat("HERE: col", col.scores[1,1], "\n")
      tmp <- c.scores[,] < -20; c.scores[,][ tmp ] <- min( c.scores[,][ ! tmp ], na.rm=T ); rm( tmp ) ## effective zero ##-220
    }
  } else c.scores <- NULL

  new.weights <- c( row=row.scaling[ iter ], mot=mot.scaling[ iter ], net=net.scaling[ iter ], set=set.scaling[ iter ] ) ##numeric()
#ifndef PACKAGE
  if ( pareto.adjust.scalings && iter > 51 ) {
    new.weights <- pareto.adjust.weights() ## iter ) ##, ... )
    ## TODO: fix pareto.adjusting (use smooth spline?) or runmed (running median) - average over multiple scales
    ## if ( FALSE ) { ## Make sure weights are never smaller than the input schedule provides.. (do we want this?)
    ##   if ( new.weights[ "mot" ] < mot.scaling[ iter ] / row.scaling[ iter ] )
    ##     new.weights[ "mot" ] <- mot.scaling[ iter ] / row.scaling[ iter ]
    ##   else if ( new.weights[ "mot" ] > new.weights[ "row" ] ) new.weights[ "mot" ] <- new.weights[ "row" ]
    ##   if ( new.weights[ "net" ] < net.scaling[ iter ] / row.scaling[ iter ] )
    ##     new.weights[ "net" ] <- net.scaling[ iter ] / row.scaling[ iter ]
    ##   else if ( new.weights[ "net" ] > new.weights[ "row" ] ) new.weights[ "net" ] <- new.weights[ "row" ]
    ## }
    if ( iter %in% stats.iters ) cat( "New weights:", new.weights, "\n" ) 
  }
#endif
  
  ## Hey, instead of just making each distrib. have the same mean and variance, let's make them have the same
  ##   exact distribution!
  if ( quantile.normalize ) {
    tmp <- list( row=r.scores, mot=m.scores, net=n.scores, set=se.scores )
    if ( sum( sapply( tmp, function( i ) ! is.null( i ) ) ) > 1 ) {
      wts <- new.weights[ ! sapply( tmp, is.null ) ]
      tmp <- quantile.normalize.scores( tmp, weights=wts )
      if ( ! is.null( r.scores ) ) r.scores[,] <- tmp$row[,]
      if ( ! is.null( m.scores ) ) m.scores[,] <- tmp$mot[,]
      if ( ! is.null( n.scores ) ) n.scores[,] <- tmp$net[,]
      if ( ! is.null( se.scores ) ) se.scores[,] <- tmp$set[,]
    }
    rm( tmp )
  }
  
  if ( new.weights[ "row" ] != 1 ) r.scores[,] <- r.scores[,] * new.weights[ "row" ] 
  if ( ! is.null( m.scores ) ) {
    tmp <- ! is.na( m.scores[,] )
    r.scores[,][ tmp ] <- r.scores[,][ tmp ] + m.scores[,][ tmp ] * new.weights[ "mot" ]
  }
  if ( ! is.null( n.scores ) ) {
    tmp <- ! is.na( n.scores[,] )
    r.scores[,][ tmp ] <- r.scores[,][ tmp ] + n.scores[,][ tmp ] * new.weights[ "net" ]
  }
  if ( ! is.null( se.scores ) ) {
    tmp <- ! is.na( se.scores[,] )
    r.scores[,][ tmp ] <- r.scores[,][ tmp ] + se.scores[,][ tmp ] * new.weights[ "set" ]
  }
  r.scores <- matrix.reference( r.scores )
  c.scores <- matrix.reference( c.scores )
  if ( ! is.null( n.scores ) ) n.scores <- matrix.reference( n.scores )
  if ( ! is.null( se.scores ) ) se.scores <- matrix.reference( se.scores )
  if ( ! is.null( m.scores ) ) m.scores <- matrix.reference( m.scores )
  invisible( list( r=r.scores, c=c.scores, n=net.scores, se=set.scores, m=mot.scores, scalings=new.weights ) )
}

## Hack to get more genes added to really small (<10 rows) clusters
get.density.scores <- function( ks=1:k.clust, plot="none", bw.scale=function( nr ) exp( -nr / 10 ) * 10 ) { ##r.scores, col.scores, 
  rr <- attr( ratios, "rnames" )
  rs <- r.scores ##[,]
  ## TODO: adjust bw to be function of iter (is this necessary with fuzzy.index turned on???)
  bw.r <- max( diff( range( rs[,], na.rm=T ) ) / 100, 0.001 ) ##"nrd0"
  get.rr.scores <- function( k ) {
    rows <- get.rows( k ) 
    cols <- get.cols( k ) 
    rsk <- rs[ ,k, drop=T ]
    if ( length( rows ) > 0 && length( cols ) > 0 ) {
      bw <- bw.r * bw.scale( length( rows ) )
      d <- density( rsk[ rows ], na.rm=T, bw=bw, adjust=2,
                   from=min( rsk, na.rm=T ) - 1, to=max( rsk, na.rm=T ) + 1, n=256 )
      p <- approx( d$x, rev( cumsum( rev( d$y ) ) ), rsk )$y
      if ( "rows" %in% plot ) { h=hist(rsk,breaks=50,main=NULL,xlab="Combined scores");tmp.scale <- round( attr( ratios, "nrow" ) / length( rows )/4 );hist(rep(rsk[rows],tmp.scale),breaks=h$breaks,col="red",border="red",add=T);hist(rsk,breaks=h$breaks,add=T);lines(d$x,d$y/max(d$y,na.rm=T)*attr(ratios,"nrow")/50,col="blue");lines(sort(rsk),p[order(rsk)]/max(p,na.rm=T)*attr(ratios,"nrow")/50,col="green") }
    } else {
      p <- rep( 1, attr( ratios, "nrow" ) )
    }
    return( p / sum( p, na.rm=T ) ) ##* length( rows )
  }

  if ( ! exists( "rr.scores" ) || is.null( rr.scores ) || any( dim( rr.scores ) != dim( r.scores ) ) ) {
    rr.scores <- row.scores[,] * 0
    rr.scores <- matrix.reference( rr.scores )
  }
  
  mc <- get.parallel( length( ks ) )
  rr.scores[,] <- do.call( cbind, mc$apply( ks, get.rr.scores ) )
  rr.scores[,][ is.infinite( rr.scores[,] ) ] <- NA
  ##rownames( rr.scores ) <- attr( ratios, "rnames" )
  
  cc.scores <- NULL ## this is required for n.clust.per.col == k.clust !!! dont delete it!
  if ( ! is.null( col.scores ) ) {
    cs <- c.scores ##[,]    
    bw.c <- max( diff( range( cs[,], na.rm=T ) ) / 100, 0.001 )
    get.cc.scores <- function( k ) {
      cols <- get.cols( k ) 
      rows <- get.rows( k ) 
      csk <- cs[ ,k, drop=T ]
      if ( length( cols ) > 0 && length( rows ) > 0 &&
          ! all( is.na( csk[ cols ] ) ) && ! all( is.infinite( csk[ cols ] ) ) &
          ! all( csk[ cols ][ ! is.na( csk[ cols ] ) ] == csk[ cols[ ! is.na( csk[ cols ] ) ][ 1 ] ] ) ) {
        d <- density( csk[ cols ], na.rm=T, from=min( csk, na.rm=T ) - 1, to=max( csk, na.rm=T ) + 1,
                     bw=bw.c, adjust=2, n=256 )
        p <- approx( d$x, rev( cumsum( rev( d$y ) ) ), csk )$y
        if ( "cols" %in% plot ) { h=hist(csk,breaks=50,main=NULL,xlab="Combined scores");tmp.scale <- round( attr( ratios, "ncol" ) / length( cols )/4 ) + 1;hist(rep(csk[cols],tmp.scale),breaks=h$breaks,col="red",border="red",add=T);hist(csk,breaks=h$breaks,add=T);lines(d$x,d$y/max(d$y,na.rm=T)*attr(ratios,"ncol")/50,col="blue");lines(sort(csk),p[order(csk)]/max(p,na.rm=T)*attr(ratios,"ncol")/50,col="green") }
      } else {
        p <- rep( 1, attr( ratios, "ncol" ) )
      }
      return( p / sum( p, na.rm=T ) ) ##* length( cols )
    }

    if ( ! exists( "cc.scores" ) || is.null( cc.scores ) || any( dim( cc.scores ) != dim( c.scores ) ) ) {
      cc.scores <- col.scores[,] * 0
      cc.scores <- matrix.reference( cc.scores )
    }

    if ( ! is.null( c.scores ) ) { ##&& ! is.na( c.scores ) ) {
      cc.scores[,] <- do.call( cbind, mc$apply( ks, get.cc.scores ) )
      cc.scores[,][ is.infinite( cc.scores ) ] <- NA
      ##rownames( cc.scores ) <- attr( ratios, "cnames" )
    }
  }
  
  invisible( list( r=rr.scores, c=cc.scores ) )
}  

## This reassigning of all rows over all clusters is too much ... problem is crappy clusters end up losing
##     ALL rows, and okay clusters get TOO MANY. IDEA: allow only 1 change per gene per iter
## TODO: When updating, don't let number of genes in a cluster get > max.cluster.rows or < min.cluster.rows
get.updated.memberships <- function() { ## rr.scores, cc.scores ) {
  ##if ( ! is.null( rr.scores ) ) {
  ##n.rows <- tabulate( row.membership )
  rm <- t( apply( rr.scores[,], 1, order, decreasing=T )[ 1:n.clust.per.row, ,drop=F ] ) ##[ iter ]
  rm <- t( apply( rm, 1, sort ) ); if ( n.clust.per.row == 1 ) rm <- t( rm ) ##[ iter ]
  ##cra <- cluster.rows.allowed

  for ( i in 1:nrow( rm ) ) {
    if ( all( rm[ i, ] %in% row.membership[ i, ] ) ) next
    mc <- max.changes[ "rows" ] ## If 0<max.changes<1, it's a "prob of seeing a change" and choose whether to change randomly
    if ( mc < 1 && mc > 0 && runif( 1 ) > mc ) next ##else mc <- 1
    for ( ii in 1:mc ) { 
    if ( sum( ! rm[ i, ] %in% row.membership[ i, ] ) >= mc ) {
      if ( any( row.membership[ i, ] == 0 ) ) {
        col.change <- which( row.membership[ i, ] == 0 )[ 1 ] ## If any are zero (no cluster assigned) use that index
      } else {
        ttmp <- tabulate( row.membership[ i, ] ) ## Any clusters have this gene assigned more than once?
        if ( any( ttmp > 1 ) ) { ## If so, use that index
          col.change <- which( row.membership[ i, ] %in% which( ttmp > 1 ) )[ 1 ]
        } else { ## Only make change if it improves the rr.score value 
          delta <- rr.scores[ i, rm[ i, ] ] - rr.scores[ i, row.membership[ i, ] ]
          if ( any( row.membership[ i, ] %in% rm[ i, ] ) ) delta[ row.membership[ i, ] %in% rm[ i, ] ] <- 0
          if ( all( is.na( delta ) | delta <= 0 ) ) next
          col.change <- which.max( delta ) 
        }
      }
      if ( exists( "maintain.seed" ) && ! is.null( maintain.seed ) && ! is.null( maintain.seed$rows ) && ## Don't remove row i from cluster
          ! is.null( maintain.seed$rows[[ as.character( row.membership[ i, col.change ] ) ]] ) &&
          rownames( row.membership )[ i ] %in%
          maintain.seed$rows[[ as.character( row.membership[ i, col.change ] ) ]] ) next 
      if ( ! rm[ i, col.change ] %in% row.membership[ i, ] ) row.membership[ i, col.change ] <- rm[ i, col.change ]
    }
  }
  }
  ##} 
  
  if ( ! is.null( cc.scores ) ) {
  ##n.cols <- tabulate( col.membership )
  cm <- t( apply( cc.scores[,], 1, order, decreasing=T )[ 1:n.clust.per.col, ,drop=F ] ) ##[ iter ]

  for ( i in 1:nrow( cm ) ) {
    mc <- max.changes[ "cols" ] ## If 0<max.changes<1, it's a "prob of seeing a change" and choose whether to change randomly
    if ( mc < 1 && mc > 0 && runif( 1 ) > mc ) next ##else mc <- 1
    for ( ii in 1:mc ) {
    if ( sum( ! cm[ i, ] %in% col.membership[ i, ] ) >= mc ) {
      if ( any( col.membership[ i, ] == 0 ) ) {
        col.change <- which( col.membership[ i, ] == 0 )[ 1 ]
      } else {
        ttmp <- tabulate( col.membership[ i, ] )
        if ( any( ttmp > 1 ) ) {
          col.change <- which( col.membership[ i, ] %in% which( ttmp > 1 ) )[ 1 ]
        } else { ## Only make change if it improves the cc.score value 
          delta <- cc.scores[ i, cm[ i, ] ] - cc.scores[ i, col.membership[ i, ] ]
          if ( all( is.na( delta ) | delta <= 0 ) ) next
          col.change <- which.max( delta ) 
        }
      }
      if ( exists( "maintain.seed" ) && ! is.null( maintain.seed ) && ! is.null( maintain.seed$cols ) && ## Don't remove col i from cluster
          ! is.null( maintain.seed$cols[[ as.character( col.membership[ i, col.change ] ) ]] ) &&
          colnames( col.membership )[ i ] %in%
          maintain.seed$cols[[ as.character( col.membership[ i, col.change ] ) ]] ) next 
      col.membership[ i, col.change ] <- cm[ i, col.change ]
    }
  } 
  }
  }
  
  ##cat("HERE:",sum(sapply(1:nrow(rm),function(i)sum(!rm[i,]%in%row.membership[i,]))),
  ##    sum(sapply(1:nrow(cm),function(i)sum(!cm[i,]%in%col.membership[i,]))),"\n")
  invisible( list( r=row.membership, c=col.membership ) )
}

seed.clusters <- function( k.clust, seed.method="rnd", col.method="rnd" ) {
  ## Allow it to be overridden by a custom function if desired (e.g. to seed from a gene list -- no that's a bad
  ## example -- there's the "list=" option below)
  if ( seed.method == "custom" && exists( "seed.clusters.custom" ) )
    return( seed.clusters.custom( k.clust, col.method ) )
  if ( substr( seed.method, 1, 3 ) == "net" && length( networks ) <= 0 ) {
    cat( "Seed method is", seed.method, ", but no networks -- using 'kmeans' instead.\n" )
    seed.method <- "kmeans"
  }
  if ( seed.method == "rnd" ) { ## RANDOM, ALL GENES ASSIGNED TO k.clust CLUSTERS
    row.membership <- t( sapply( 1:attr( ratios, "nrow" ),
                                function( i ) sample( 1:k.clust, n.clust.per.row[ 1 ],
                                                     replace=n.clust.per.row[ 1 ] > attr( ratios, "nrow" ) ) ) )
  } else if ( substr( seed.method, 1, 5 ) == "list=" ) { ## File with sets of genes - one set of genes per line
    row.membership <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=n.clust.per.row[ 1 ] ) ## separated by space,tab,or comma
    rownames( row.membership ) <- attr( ratios, "rnames" )
    fname <- strsplit( seed.method, "=" )[[ 1 ]][ 2 ]
    if ( exists( fname ) ) lists <- get( fname ) ## List of vectors already exists in memory
    else if ( file.exists( fname ) ) lists <- strsplit( readLines( fname ), split="[ \t,]", perl=T )
    for ( k in 1:min( c( k.clust, length( lists ) ) ) ) {
      probes <- unlist( lapply( get.synonyms( lists[[ k ]] ), function( i ) i %in% rownames( row.membership ) ) )
      row.membership[ probes[ row.membership[ probes, 1 ] == 0 ], 1 ] <- k
      row.membership[ probes[ row.membership[ probes, 1 ] != 0 ], 2 ] <- k
    }
    if ( length( lists ) < k.clust ) { ## Fill in additional clusters randomly
      for ( k in ( length( lists ) + 1 ):k.clust ) {
        rnames <- attr( ratios, "rnames" )[ ! attr( ratios, "rnames" ) %in% unlist( lists ) ]
        rows <- sample( rnames, 5 )
        row.membership[ rows[ row.membership[ rows, 1 ] == 0 ], 1 ] <- k
        row.membership[ rows[ row.membership[ rows, 1 ] != 0 & row.membership[ rows, 2 ] == 0 ], 2 ] <- k
      }
    }
  } else if ( substr( seed.method, 1, 4 ) == "rnd=" ) { ## RANDOM SAMPLED N per cluster
    n.samp <- as.integer( strsplit( seed.method, "=" )[[ 1 ]][ 2 ] )
    row.membership <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=n.clust.per.row[ 1 ] )
    rownames( row.membership ) <- attr( ratios, "rnames" )
    for ( i in 1:n.clust.per.row ) {
      sampled <- rep( FALSE, attr( ratios, "nrow" ) ); names( sampled ) <- attr( ratios, "rnames" )
      for ( k in 1:k.clust ) {
        g <- sample( attr( ratios, "rnames" )[ ! sampled ], n.samp )
        row.membership[ g, 1 ] <- k
        sampled[ g ] <- TRUE
      }
    }
  } else if ( seed.method == "kmeans" ) { ## Kmeans seeded
    if ( ! exists( "ratios" ) ) stop( "kmeans seed method but no ratios" )
    tmp.rat <- get.cluster.matrix() ##ratios;
    tmp.rat[ is.na( tmp.rat ) ] <- 0
    ##cat(dim(tmp.rat),k.clust,"\n")
    km <- kmeans
#ifndef PACKAGE
    if ( ! is.na( parallel.cores ) && ( parallel.cores == TRUE || parallel.cores > 1 ) ) {
        ##&& require( biganalytics ) && require( doMC ) )
        get.parallel(); if ( require( biganalytics ) ) km <- bigkmeans }
#endif
    row.membership <- km( tmp.rat, centers=k.clust, iter.max=20, nstart=2 )$cluster
    names( row.membership ) <- attr( ratios, "rnames" )
    if ( n.clust.per.row[ 1 ] > 1 ) row.membership <-
      cbind( row.membership, matrix( rep( 0, attr( ratios, "nrow" ) * ( n.clust.per.row[ 1 ] - 1 ) ),
                                         ncol=n.clust.per.row[ 1 ] - 1 ) )
#ifndef PACKAGE
  } else if ( substr( seed.method, 1, 11 ) == "trimkmeans=" ) { ## Trimmed k-means seeded
    if ( ! exists( "ratios" ) ) stop( "trimkmeans seed method but no ratios" )
    require( trimcluster )
    trim <- as.numeric( strsplit( seed.method, "=" )[[ 1 ]][ 2 ] )
    tmp.rat <- get.cluster.matrix() ##ratios;
    tmp.rat[ is.na( tmp.rat ) ] <- 0
    row.membership <- trimkmeans( tmp.rat, k.clust, trim=trim, maxit=20, runs=2 )$classification
    if ( n.clust.per.row[ 1 ] > 1 ) row.membership <-
      cbind( row.membership, matrix( rep( 0, attr( ratios, "nrow" ) * ( n.clust.per.row[ 1 ] - 1 ),
                                         ncol=n.clust.per.row[ 1 ] - 1 ) ) )
  } else if ( substr( seed.method, 1, 4 ) == "cor=" ) { ## RANDOM GENES AND THEIR N CORRELATED GENES
    if ( ! exists( "ratios" ) ) stop( "cor seed method but no ratios" )
    n.cor <- as.integer( strsplit( seed.method, "=" )[[ 1 ]][ 2 ] )
    rats <- get.cluster.matrix()
    cors <- if ( attr( ratios, "nrow" ) < 6000 ) cor( t( rats ), use="pairwise" ) else NULL
    rm <- rep( 0, attr( ratios, "nrow" ) ); names( rm ) <- attr( ratios, "rnames" )
    sampled <- rep( FALSE, attr( ratios, "nrow" ) ); names( sampled ) <- attr( ratios, "rnames" )
    mc <- get.parallel( n.clust.per.row )
    tmp <- mc$apply( 1:n.clust.per.row, function( i ) {
      for ( k in 1:k.clust ) {
        if ( sum( ! sampled ) < n.cor ) sampled[ sample( 1:length( sampled ) ) ] <- FALSE
        rnames <- attr( ratios, "rnames" )[ ! sampled ]
        g <- sample( rnames, 1 )
        if ( ! is.null( cors ) ) g <- rnames[ order( cors[ g, ! sampled ], decreasing=T )[ 1:n.cor ] ]
        else g <- rnames[ order( apply( rats[ ! sampled, ], 1, cor, rats[ g, ] ), decreasing=T )[ 1:n.cor ] ]
        rm[ g ] <- k
        if ( length( g ) == 1 ) { ## if option was cor=1, try to spread the wealth by preventing genes correlated 
          if ( ! is.null( cors ) ) tmp <- rnames[ order( cors[ g, ! sampled ], decreasing=T )[ 1:10 ] ] ## with this one from getting
          else tmp <- rnames[ order( apply( rats[ ! sampled, ], 1, cor, rats[ g, ] ), decreasing=T )[ 1:10 ] ]
          sampled[ tmp ] <- TRUE ## added to other clusters
        }
        sampled[ g ] <- TRUE
      }
      rm[ attr( ratios, "rnames" ) ]
    } ) 
    row.membership <- do.call( cbind, tmp )
  } else if ( substr( seed.method, 1, 4 ) == "net=" ) { ## Network seeded; e.g. 'net=string:5'
    if ( ! exists( "networks" ) || length( networks ) <= 0 ) stop( "net seed method but no networks" )
    seed.method <- strsplit( seed.method, "=" )[[ 1 ]][ 2 ] ## Name of network to use and # of seed nodes
    net.name <- strsplit( seed.method, ":" )[[ 1 ]][ 1 ]
    net <- networks[[ net.name ]]
    n.seed <- as.integer( strsplit( seed.method, ":" )[[ 1 ]][ 2 ] )
    rm <- rep( 0, attr( ratios, "nrow" ) ); names( rm ) <- attr( ratios, "rnames" )
    sampled <- rep( FALSE, length( unique( as.character( net$protein1 ) ) ) )
    names( sampled ) <- unique( as.character( net$protein1 ) )
    mc <- get.parallel( n.clust.per.row )
    tmp <- mc$apply( 1:n.clust.per.row, function( i ) {
      for ( k in 1:k.clust ) {
        if ( sum( ! sampled ) <= 0 ) sampled[ 1:length( sampled ) ] <- FALSE
        rnames <- names( which( ! sampled ) )
        gs <- sample( rnames, 1 )
        qiter <- 0
        while( length( gs ) < n.seed && qiter < 20 ) {
          ns <- as.character( net$protein2[ as.character( net$protein1 ) %in% gs ] )
          if ( length( ns ) + length( gs ) >= n.seed ) ns <- sample( ns, size=n.seed-length( gs ),
                                      prob=net$combined_score[ as.character( net$protein1 ) %in% gs ] )
          gs <- unique( c( gs, ns ) )
          qiter <- qiter + 1
        }
        rm[ gs ] <- k
        sampled[ gs ] <- TRUE  ## try to spread the wealth by preventing genes connected to these from being added 
        if ( n.seed <= 2 ) sampled[ as.character( net$protein2[ as.character( net$protein1 ) %in% gs ] ) ] <- TRUE ## to other clusts
      }
      rm[ attr( ratios, "rnames" ) ]
    } ) 
    row.membership <- do.call( cbind, tmp )
  } else if ( substr( seed.method, 1, 7 ) == "netcor=" ) { ## Correld net-weighted seeds, e.g. 'netcor=string:5'
    if ( ! exists( "ratios" ) ) stop( "netcor seed method but no ratios" )
    if ( ! exists( "networks" ) || length( networks ) <= 0 ) stop( "netcor seed method but no networks" )
    seed.method <- strsplit( seed.method, "=" )[[ 1 ]][ 2 ] ## Name of network to use and # of seed nodes
    net.name <- strsplit( seed.method, ":" )[[ 1 ]][ 1 ]
    net <- networks[[ net.name ]]
    n.seed <- as.integer( strsplit( seed.method, ":" )[[ 1 ]][ 2 ] )
    rats <- get.cluster.matrix()
    cors <- cor( t( rats ), use="pairwise" )
    tmp.mat <- matrix( 0, nrow=nrow( cors ), ncol=ncol( cors ) ); dimnames( tmp.mat ) <- dimnames( cors )
    tmp.lookup <- 1:attr( ratios, "nrow" ); names( tmp.lookup ) <- attr( ratios, "rnames" )
    net <- net[ as.character( net$protein1 ) %in% attr( ratios, "rnames" ) &
               as.character( net$protein2 ) %in% attr( ratios, "rnames" ), ]
    tmp.mat[ cbind( tmp.lookup[ as.character( net$protein1 ) ], tmp.lookup[ as.character( net$protein2 ) ] ) ] <-
      net$combined_score / 1000
    cors <- cors + tmp.mat; rm( tmp.mat )
    rm <- rep( 0, attr( ratios, "nrow" ) ); names( rm ) <- attr( ratios, "rnames" )
    sampled <- rep( FALSE, attr( ratios, "nrow" ) ); names( sampled ) <- attr( ratios, "rnames" )
    mc <- get.parallel( n.clust.per.row )
    tmp <- mc$apply( 1:n.clust.per.row, function( i ) {
      for ( k in 1:k.clust ) {
        if ( sum( ! sampled ) < n.seed ) sampled[ sample( 1:length( sampled ) ) ] <- FALSE
        rnames <- attr( ratios, "rnames" )[ ! sampled ]
        g <- sample( rnames, 1 )
        g <- rnames[ order( cors[ g, ! sampled ], decreasing=T )[ 1:n.seed ] ]
        rm[ g ] <- k
        if ( length( g ) == 1 ) { ## if option was cor=1, try to spread the wealth by preventing genes correlated 
          tmp <- rnames[ order( cors[ g, ! sampled ], decreasing=T )[ 1:10 ] ] ## with this one from getting 
          sampled[ tmp ] <- TRUE ## added to other clusters
        }
        sampled[ g ] <- TRUE
      }
      rm[ attr( ratios, "rnames" ) ]
    } ) 
    row.membership <- do.call( cbind, tmp )
#endif
  }
  
  if ( is.vector( row.membership ) ) row.membership <- t( row.membership ) ## probably n.clust.per.row == 1
  if ( nrow( row.membership ) == 1 ) row.membership <- t( row.membership ) ## probably n.clust.per.row == 1
 
  if ( col.method == "rnd" ) {
    col.membership <- t( sapply( 1:attr( ratios, "ncol" ), function( i )
                                sample( 1:k.clust, n.clust.per.col[ 1 ],
                                       replace=n.clust.per.col[ 1 ] > k.clust ) ) ) ##attr( ratios, "ncol" ) ) ) )
  } else if ( col.method == "best" ) {
    if ( ! exists( "ratios" ) ) stop( "best col seed method but no ratios" )
    all.rats <- get.cluster.matrix()
    attr( all.rats, "all.colVars" ) <- apply( all.rats, 2, var, use="pair", na.rm=T )
    col.scores <- -sapply( 1:k.clust, function( k )
                          if ( sum( row.membership == k, na.rm=T ) <= 0 ) rep( NA, attr( ratios, "ncol" ) ) else 
                          get.col.scores( k=get.rows( k, row.membership ), ratios=all.rats, method="orig" ) ) 
    col.membership <- t( apply( col.scores, 1, function( i ) order( i, decreasing=T )[ 1:n.clust.per.col[ 1 ] ] ) )
  }

  rownames( row.membership ) <- attr( ratios, "rnames" ) 
  rownames( col.membership ) <- attr( ratios, "cnames" )
  list( row.membership=row.membership, col.membership=col.membership )
}  

## SEEDING -- seed new clusters if no row.membership exists.
##    ALSO re-set the random seed if env$cmonkey.params$rnd.seed doesn't exist (e.g. for multiple runs on the
##       same pre-initialized env.) -- need to this via:
##       rm(list="rnd.seed",envir=env$cmonkey.params);rm(list="rnd.seed",envir=env)
cmonkey.re.seed <- function( env ) {
  if ( ! exists( "rnd.seed", envir=env$cmonkey.params ) ) {
    ##tmp.rnd.seed <- as.integer( Sys.time() ) ## This will only actually change the seed if
    op <- options( digits.secs=10 )
    tmp.time <- as.character( Sys.time() )
    options( op ); rm( op )
    tmp.rnd.seed <- as.integer( substr( gsub( '[-:. ]', "", tmp.time ), 12, 20 ) )
    cat( "RESETTING RANDOM SEED: " ) ## env$rnd.seed and env$cmonkey.params$rnd.seed are GONE
    env$set.param( "date.run", tmp.time, env$cmonkey.params ) 
    env$date.run <- env$cmonkey.params$date.run
    env$set.param( "rnd.seed", tmp.rnd.seed, env$cmonkey.params ) 
    env$rnd.seed <- env$cmonkey.params$rnd.seed
    set.seed( env$rnd.seed )
    rm( tmp.rnd.seed )
    ##env$set.param( "cmonkey.filename", paste( "cmonkey", env$cmonkey.version, env$organism,
    ##                                         paste( sapply( env$ratios, dim ), collapse="x" ),
    ##                                         gsub( " ", "_", env$date.run ), sep="_" ), env$cmonkey.params )
    ##env$cmonkey.filename <- env$cmonkey.params$cmonkey.filename
  }
  if ( ! is.null( env$ratios ) && attr( env$ratios, "ncol" ) > 1 ) {
    cat( "Seeding all clusters using methods:", env$seed.method, "\n" ) ##seed.method.cols, "\n" )
    tmp <- env$seed.clusters( env$k.clust, seed.method=env$seed.method[ "rows" ], col.method=env$seed.method[ "cols" ] ) ##.cols )
  } else {
    cat( "Seeding all clusters using methods: rnd rnd\n" )
    tmp <- env$seed.clusters( env$k.clust, seed.method="rnd", col.method="rnd" )
  }
  env$row.membership <- tmp$row.membership; env$col.membership <- tmp$col.membership; rm( tmp )
  invisible( env )
}

get.cluster.matrix <- function( rows=NULL, cols=NULL, matrices=names( ratios ) ) {
  if ( is.null( rows ) ) rows <- attr( ratios, "rnames" )
  if ( is.null( cols ) ) cols <- attr( ratios, "cnames" )
  cols.b <- attr( ratios, "cnames" )[ attr( ratios, "cnames" ) %in% cols ]

  rats <- matrix( NA, nrow=length( rows ), ncol=length( cols.b ) )
  rownames( rats ) <- rows; colnames( rats ) <- cols.b
  cnames <- character()
  for ( n in matrices ) { ## names( ratios ) ) {
    r.tmp <- ratios[[ n ]][ rows[ rows %in% rownames( ratios[[ n ]] ) ],
                           cols.b[ cols.b %in% colnames( ratios[[ n ]] ) ], drop=F ]
    if ( is.null( r.tmp ) || all( is.na( r.tmp ) ) ) next
    if ( is.vector( r.tmp ) ) { r.tmp <- t( r.tmp ); rownames( r.tmp ) <- rows }
    cnames <- c( cnames, colnames( r.tmp ) )
    rats[ rownames( r.tmp ), colnames( r.tmp ) ] <- r.tmp; rm( r.tmp )
  }
  rats[ ,colnames( rats ) %in% cnames, drop=F ]
}

## Get all potential synonyms from feature.names table and possibly from "translation.tab" if it exists
get.synonyms <- function( gns, ft=genome.info$feature.names, ignore.case=T, verbose=F, fast=F, force=F ) { ##transl.table
  if ( exists( "no.genome.info" ) && no.genome.info ) { out <- as.list( gns ); names( out ) <- gns; return( out ) }
  out <- list()
  if ( ! force && exists( "genome.info" ) && ! is.null( genome.info$synonyms ) ) {
    gns.cached <- gns[ gns %in% names( genome.info$synonyms ) ]
    out <- genome.info$synonyms[ gns.cached ]
    gns <- gns[ ! gns %in% names( genome.info$synonyms ) ]
    if ( length( gns ) <= 0 ) return( out )
  }
  ##cat( "Generating gene synonyms cache... this could take a little while.\n" )
  tmp.out <- as.list( gns ); names( tmp.out ) <- gns
  if ( is.null( ft ) ) return( tmp.out ) ##character() )
  gns.orig <- gns
  gns <- gsub( "m$|_\\d$|\\-\\S$", "", gns, perl=T ) ## specific for Halo, also works for yeast
  gns <- gsub( "([\\[\\]\\(\\)\\{\\}\\.\\+\\-'\"])", "\\\\\\1", gns, perl=T ) ## Escape the regex-y characters
  gns <- gns[ ! is.na( gns ) & gns != "" ]
  ft <- ft[ ,c( "id", "names" ) ]
  if ( exists( "translation.tab" ) && ! is.null( translation.tab ) )
    ft <- rbind( ft, data.frame( id=as.character( translation.tab$V1 ), names=as.character( translation.tab$V2 ) ) )
  ft <- subset( ft, names != "" )

  if ( verbose ) ggggg <- gns[ seq( 100, length( gns ), by=100 ) ]

  mc <- get.parallel( length( gns ), verbose=F )
  tmp <- mc$apply( gns, function( g ) {
    if ( verbose && g %in% ggggg ) cat( " ...", g )
    greg <- paste( "^", g, sep="" )
    tmp <- subset( ft, grepl( greg, id, perl=T, ignore=ignore.case ) |
                  grepl( greg, names, perl=T, ignore=ignore.case ) ) ## ) )
    if ( nrow( tmp ) <= 0 ) return( g ) ##tmp.out ) ##character() )
    tmp2 <- unique( c( g, as.character( tmp$id ), as.character( tmp$names ) ) )
    if ( ! fast ) {
      tmp2 <- subset( ft, id %in% tmp2 | names %in% tmp2 )
      tmp2 <- unique( c( g, as.character( tmp2[ ,1 ] ), as.character( tmp2[ ,2 ] ) ) )
    }
    tmp2 <- gsub( "\\\\([\\[\\]\\(\\)\\{\\}\\.\\+\\-'\"])", "\\1", tmp2, perl=T )
    tmp2
  } ) 

  names( tmp ) <- gns.orig
  if ( verbose ) cat( "\n" )
  c( tmp, out )
}

## Get id's from feature.names table possibly translating via "translation.tab" or other rules if necessary
## get.ids <- function( rows, ft=genome.info$feature.names, ... ) {
##   ids <- subset( ft, grepl( rows, names, ignore=T, perl=T ), c( id, names ) )
##   if ( exists( "translation.tab" ) && ! is.null( translation.tab ) &&
##       nrow( ids ) < length( rows ) ) { ## Use translation table first
##     rows2 <- as.character( subset( translation.tab, V1 %in% rows, V2 )[ ,1 ] )
##     ids <- rbind( ids, subset( ft, names %in% rows2, c( id, names ) ) )
##   }
##   unique( ids )
## }

get.gene.coords <- function( rows, op.shift=T, op.table=genome.info$operons, ... ) { ##, override=F ) {
  rows <- unique( rows )
  syns <- get.synonyms( rows, ... )
  tab <- genome.info$feature.tab
  ##if ( "chrom_position" %in% colnames( tab ) ) tab <- tab[ ! duplicated( tab$chrom_position ), ]
  ids <- lapply( syns, function( s ) s[ s %in% tab$id ] )
  if ( all( sapply( ids, length ) < 1 ) ) {
    warning( "Could not find gene start/stop for any input genes", call.=F ); return( NULL ) }
  if ( any( sapply( ids, length ) < 1 ) ) warning( "Could not find gene start/stop for all input genes", call.=F )
  ids <- ids[ sapply( ids, length ) >= 1 ]
  ids <- sapply( ids, "[", 1 )
  ids <- data.frame( id=ids, names=names( ids ) ) ##, gene=rows )
  
  coos <- NULL
  if ( op.shift ) { ## Replace each gene's coords with the head of its operon, but keep the "locus_tag".
    if ( attr( op.table, "source" ) == "rsat" ) {
      ops <- merge( ids, op.table, by.x="id", by.y="query", all=F )
      ops2 <- ops[ order( ops$lead ), ]
      coos <- merge( ops, tab, by.x="lead", by.y="name", all=F )[ ,c( "id.x", "names", "contig",
                                                                               "strand", "start_pos", "end_pos" ) ]
      ## Due to not-unique mapping of short names in operons table to coords, we need to get rid of dupes
      ##    in an intelligent way:
##       if ( FALSE ) { ## BIG TODO!!!!
##         dupes <- as.character( coos$names[ which( duplicated( coos$names ) ) ] )
##         syns <- get.synonyms( dupes )
##         for ( s in names( syns ) ) {
##           tmp <- subset( genome.info$feature.tab, name %in% syns[[ s ]] )
##           tmp2 <- subset( coos, id %in% tmp$id )
##         }
##       }
    } else if ( attr( op.table, "source" ) == "microbes.online" ) {
      ops <- merge( ids, op.table, by.x="names", by.y="gene", all.x=T )
      if ( any( is.na( ops$head ) ) ) {
        head <- as.character( ops$head ); head[ is.na( head ) ] <- as.character( ops$names[ is.na( head ) ] )
        ops$head <- as.factor( head )
      }
      head.syns <- get.synonyms( unique( as.character( ops$head ) ) )
      head.ids <- lapply( head.syns, function( s ) s[ s %in% tab$id ] )
      head.ids <- head.ids[ sapply( head.ids, length ) >= 1 ]
      head.ids <- data.frame( id=sapply( head.ids, "[", 1 ), names=names( head.ids ) )
      ops2 <- merge( ops, head.ids, by.x="head", by.y="names", all.x=T )
      coos <- merge( ops2, tab, by.x="id.y", by.y="id", all.x=T )[ ,c( "id.x", "names", "contig",
                                                                               "strand", "start_pos", "end_pos" ) ]
    }
  } else { ##if ( ! op.shift )
    coos <- merge( ids, tab, by="id" )[ ,c( "id", "names", "contig",
                                                   "strand", "start_pos", "end_pos" ) ]
  }
  colnames( coos )[ 1 ] <- "id"
  if ( is.factor( coos$start_pos ) ) coos$start_pos <- as.numeric( levels( coos$start_pos ) )[ coos$start_pos ]
  if ( is.factor( coos$end_pos ) ) coos$end_pos <- as.numeric( levels( coos$end_pos ) )[ coos$end_pos ]
  coos[ ! duplicated( coos[ ,1:4 ] ), ]
}  

get.long.names <- function( k, shorter=F ) {
  if ( is.numeric( k[ 1 ] ) ) { rows <- get.rows( k ) }
  else { rows <- k }
  if ( is.null( genome.info$feature.tab ) ) {
    out <- rep( "", length( rows ) ); names( out ) <- rows; return( rows ) }
  ids <- get.synonyms( rows )
  mc <- list( apply=lapply )
  if ( ! shorter ) desc <- mc$apply( ids, function( i ) subset( genome.info$feature.tab, id %in% i, 
                                                               select=c( "id", "description" ) ) )
  else {
    desc <- mc$apply( ids, function( i ) subset( genome.info$feature.tab, id %in% i, select=c( "id", "name", "description" ) ) )
    for ( i in 1:length( desc ) ) if ( length( desc[[ i ]]$name ) > 0 && desc[[ i ]]$name %in% rows ) {
      if ( grepl( "(", desc[[ i ]]$description, fixed=T ) ) ## Try to parse out short name from description
        desc[[ i ]]$name <- strsplit( as.character( desc[[ i ]]$description ), "[()]", perl=T )[[ 1 ]][ 2 ]
    }
  }
  out <- sapply( desc, function( i ) as.character( i[ 1, 2 ] ) )
  out <- out[ rows ]
  names( out ) <- rows
  out[ is.na( out ) | out == names( out ) ] <- "" 
  out
}

extend.vec <- function( v, n=n.iter ) {
  if ( length( v ) < n ) v <- c( v, rep( v[ length( v ) ], n.iter - length( v ) ) ); v }

#ifndef PACKAGE
if ( file.exists( "cmonkey-init.R" ) ) {
  source( "cmonkey.R", local=T )
  source( "cmonkey-init.R", local=T )
  source( "cmonkey-data-load.R", local=T ) ## Functions for loading the data
  source( "cmonkey-motif.R", local=T ) ## Functions for motif finding/scoring
  source( "cmonkey-plotting.R", local=T ) ## Functions for all cmonkey plotting
  source( "cmonkey-postproc.R", local=T ) ## Functions for all post-processing and analysis of cmonkey clusters
  source( "cmonkey-motif-other.R", local=T ) ## Functions for motif finding with other algorithms and blasting
  source( "cmonkey-bigmem.R", local=T ) ## Functions for using on-disk list and matrix storage for big organisms
}
#endif
