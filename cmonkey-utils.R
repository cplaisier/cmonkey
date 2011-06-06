source( "cmonkey.R" )

#!ifndef PACKAGE
cm.package <- function( data=F, bigdata=F, install=T, update.web=F, check=F, version=cm.version, beta=F ) {
  ## Can get halo ratios and translation.tab and other relevant data via data(halo);attach(halo)
  source.files <- list.files( patt="^cmonkey.*\\.R$" )
  source.files <- source.files[ ! source.files %in% c( "cmonkey.code.update.R", "cmonkey-experimental.R",
                                                      "cmonkey-motif-other.R", "cmonkey-utils.R" ) ]
  ##cm.name <- "cMonkey"
  if ( beta ) {
    ##cm.name <- paste( cm.name, "beta", sep="." ) ## Include experimental functions
    ##version <- paste( version, "99", sep="." ) ## Include experimental functions
    update.web <- FALSE ## DO NOT PUT ON THE WEB! bigdata <- 
  }
  
  onLoad <- function( libname, pkgname ) { ##.onAttach
    packageStartupMessage( "Loading ", pkgname, " version ", VERSION, " (", DATE, ")" )
    packageStartupMessage( "Copyright (C) David J Reiss, Institute for Systems Biology; dreiss@systemsbiology.org." )
    packageStartupMessage( "http://baliga.systemsbiology.net/cmonkey" )
    vers <- try( readLines( "http://baliga.systemsbiology.net/cmonkey/VERSION" ), silent=T )
    if ( class( vers ) != "try-error" ) {
      vers <- gsub( " ", "", vers )
      if ( vers != VERSION ) packageStartupMessage( "\nYou are not using the most current version of cMonkey.\nPlease consider upgrading to v", vers, " via:\n\n> download.file( \"http://baliga.systemsbiology.net/cmonkey/cMonkey_", vers, ".tar.gz\", \n\t\t\"cMonkey_", vers, ".tar.gz\" )\n> install.packages( \"cMonkey_", vers, ".tar.gz\", repos=NULL )\n\nOr by following the instructions on the cMonkey website." )
      else packageStartupMessage( "Congratulations! You are using the latest version of cMonkey.\n" )
    } else {
      packageStartupMessage( "Could not check to see if you are using the latest version of cMonkey." )
    }
  }
  
  source( "../generic_scripts/construct.package.R" )

  ## Note: eventually, want to include "ref" package, turn off the "gsubs" option, make "matrix.reference"
  ##    return ref(m), etc.. These are currently all turned off for the official package.
  construct.package( "cMonkey", version=version, source.files=source.files, nocpp=beta,
                    functions.visible=c( "cmonkey", "cmonkey.init" ), ##, "plotClust", "cluster.summary", "plotStats",
                      ##"get.rows", "get.cols", "clusters.w.func", "clusters.w.conds", "clusters.w.genes", 
                      ##"write.project", "update.cmonkey.env", "save.cmonkey.env" ),
                    functions.excluded=if ( beta ) "NONE" else
                    c( "cm.package", "construct.package", "pssm.to.consensus",
                      ### FUNCTIONS NOT YET NEEDED IN OFFICIAL PACKAGE:
                      "get.prolinks.links", "get.predictome.links", "load.ratios.GEO", "load.ratios.MicrobesOnline",
                      ### OLD FUNCTIONS:
                      "get.STRING.links.OLD", "pareto.adjust.weights.OLD",
                      ### EXPERIMENTAL FUNCTIONS:
                      "get.condition.groups", "get.col.weights", "get.row.weights", 
                      "filter.updated.memberships", "pareto.adjust.weights", 
                      "consolidate.duplicate.clusters", "re.seed.empty.clusters",
                      "row.col.membership.from.clusterStack",
                      ##"matrix.reference", 
                      "list.reference", "ffify.env", "un.ffify.env",
                      ## "get.STRING.links.NEW",
                      ### EXPERIMENTAL POST PROCESSING FUNCTIONS:
                      "cluster.GO.annotations", "cluster.GO.pvalues", "cluster.KEGG.pvalues",
                      "motif.similarities.tomtom", "desymmetrize.tomtom.results", "cluster.tomtom.results",
                      "motif.similarities.custom", "cluster.custom.results", "compare.pssms",
                      "process.history", "compare.clusters", "update.cmonkey.env", "adjust.clust.2",
                      "residual.pvalue", "pssm.motif.lines", "cluster.score.pvalues", "bicluster.pvalues",
                      ### NEW (EXPERIMENTAL) MOTIF FUNCTIONS:
                      "get.sequence.psps", "get.row.scores.NEW",
                      "weeder.one.cluster", "spacer.one.cluster",
                      "cosmo.one.cluster", "glam.one.cluster", "gadem.one.cluster", 
                      "blast.align", "parse.blast.out", "blast.match.seqs", "all.dna.seqs",
                      ### NEW OUTPUT FUNCTIONS:
                      "write.bicluster.network"
                      ),
                    ##data=list( halo="halo", hpy="hpy", mpn="mpn" ), ##ecoli="ecoli", yeast="yeast", ath="ath",
                    ##objects.included=c( cm.init="cm.init", cm.main="cm.main" ),
                    required=if ( beta ) c( "ref", "bigmemory", "filehash" ) else NULL,
                    suggested=c( "RCurl", "doMC", "foreach", "multicore", "igraph", "RSVGTipsDevice", "hwriter",
                      "ref", "bigmemory", "filehash" ), ##"fUtilities", , "Cairo", "trimcluster"
                    short.desc="cMonkey intgrated biclustering algorithm",
                    long.desc="cMonkey integrated biclustering algorithm for learning co-regulated gene modules",
                    url="http://baliga.systemsbiology.net/cmonkey",
                    reference='"Integrated biclustering of heterogeneous genome-wide datasets for the inference of global regulatory networks",\nby David J Reiss, Nitin S Baliga and Richard Bonneau: \\url{http://www.biomedcentral.com/1471-2105/7/280}',
                    onLoad=onLoad, gsubs=if ( beta ) list( c( "[,]", "" ), c( "[]", "" ) ) else NULL )

  pkg.name <- sprintf( "lib/cMonkey_%s.tar.gz", version )
  if ( beta ) {
    system( sprintf( "mv %s lib/cMonkey_%s_beta.tar.gz", pkg.name, version ) )
    pkg.name <- sprintf( "lib/cMonkey_%s_beta.tar.gz", version )
  }
  if ( install ) system( sprintf( "R CMD INSTALL %s", pkg.name ) )

  if ( check ) {
    cwd <- setwd( "lib" )
    system( sprintf( "R CMD CHECK %s", pkg.name ) )
    setwd( cwd )
  }
  
  if ( data ) {
    ## Halo data goes in data package
    cat( "Packaging Halo data...\n" )
    halo <- list( organism="hal", cog.org="Hbs", rsat.species="Halobacterium_sp", taxon.id=64091, k.clust=250 )
    load( "data/halo/halo.rats.RData" ); halo$ratios <- ratios
    halo$translation.tab <- read.delim( "data/halo/vng7000_to_vng5000.tsv" )
    colnames( halo$translation.tab ) <- c( "V1", "V2" )
    halo$bat.clust.genes <- c("VNG0654C","VNG0656H","VNG0828H","VNG0829G","VNG0830G","VNG0831G","VNG0832C","VNG1039H",
                              "VNG1370G","VNG1405C","VNG1458G","VNG1459H","VNG1461H","VNG1462G","VNG1463G","VNG1464G",
                              "VNG1465G","VNG1467G","VNG1468H","VNG1628G","VNG1630H","VNG1655H","VNG1656H","VNG1657H",
                              "VNG1755G","VNG1882G","VNG1883G","VNG1884G","VNG2137G","VNG2535H")
    halo$cm.func.each.iter <- function(){if(!iter%in%plot.iters)return();fc<-clusters.w.genes(bat.clust.genes);
                                         cat("HERE:",which.max(fc),max(fc),"\n")}
    halo$favorite.cluster <- function() which.max( clusters.w.genes( bat.clust.genes ) )
    halo$net.weights <- c( string=1 ) ## prolinks=1 ##0.5, operons=0.5 )
    halo <<- halo

    ## Hpy data goes in data package
    cat( "Packaging Hpy data...\n" )
    hpy <- list( organism="hpy", k.clust=75 )
    load( "data/hpy/hpy.rats.RData" ); hpy$ratios <- ratios
    source( "cmonkey-data-load.R" ); hpy$pp.ints <- load.sif.interactions( "data/hpy/HPY-PPINTS.sif.gz" )
    hpy$cm.func.each.iter=function(){if(!iter%in%plot.iters)return();fc<-clusters.w.func("flagell");##flag.clusters();
                                     cat("FLAG CLUSTER:",which.max(fc),max(fc),"\n")}
    hpy$favorite.cluster <- function() which.max( clusters.w.func( "flagell" ) )
    hpy$net.weights <- c( string=1, pp.ints=1 ) ## operons=0.5, ##, `data/hpy/HPY-PPINTS.sif.gz`=1 )
    hpy <<- hpy

    ## Mycoplasma data goes in data package
    cat( "Packaging Myc data...\n" )
    mpn <- list( organism="mpn", k.clust=75 )
    ratios <- read.delim( gzfile( "data/Mycoplasma_pneumoniae/GSE14015_series_matrix.txt.gz" ), comment='!' )
    cnames <- grep( "^!Sample_title", readLines( gzfile( "data/Mycoplasma_pneumoniae/GSE14015_series_matrix.txt.gz" ),
                                                n=100 ), val=T )
    cnames <- gsub( "\302\272", "deg", gsub( "\"", "", strsplit( cnames, "\t" )[[ 1 ]] ) )
    rownames( ratios ) <- sapply( strsplit( as.character( ratios[[ 1 ]] ), "|", fixed=T ), "[", 1 )
    colnames( ratios ) <- cnames; rm( cnames )
    ratios <- as.matrix( ratios[ ,-1 ] )
    mpn$ratios <- ratios
    mpn <<- mpn
    
    source( "../generic_scripts/construct.package.R" )
    construct.package( "cMonkey.data", version=version, source.files=NULL, ##source.files,
                      data=list( halo="halo", hpy="hpy", mpn="mpn" ), ##required="cMonkey",
                      short.desc="Additional example data sets for use with the cMonkey integrated biclustering algorithm",
                      long.desc="Additional example data sets for use with the cMonkey integrated biclustering algorithm" )

    if ( install ) system( sprintf( "R CMD INSTALL lib/cMonkey.data_%s.tar.gz", version ) )
  }
  
  if ( bigdata ) { 
    cat( "Packaging Eco data...\n" ) ## Ecoli data package goes in "bigdata" package
    ecoli <- list( organism="eco", k.clust=300, rsat.species="Escherichia_coli_K12" )
    ecoli$ratios.old <- read.delim( file=gzfile( "data/ecoli/ratios.table.txt-cmpaper.gz" ), sep="\t", as.is=T, header=T )
    ecoli$ratios <- as.matrix( read.delim( file=gzfile( "data/ecoli/ratios.table.txt.gz" ), sep="\t", as.is=T, header=T ) )
    ecoli$cm.func.each.iter=function(){if(!iter%in%plot.iters)return();fc<-clusters.w.func("flagell");##flag.clusters();
                                    cat("FLAG CLUSTER:",which.max(fc),max(fc),"\n")}
    ecoli$favorite.cluster <- function() which.max( clusters.w.func( "flagell" ) )
    ecoli <<- ecoli

    cat( "Packaging Sce data...\n" ) ## Yeast data package goes in "bigdata" package
    yeast <- list( organism="sce", is.eukaryotic=T, k.clust=450, motif.upstream.search=c( -30, 250 ),
                  motif.upstream.scan=c( -30, 500 ) )
    yeast$meme.cmd <- "./progs/meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 12 -mod zoops -nostatus -text -cons $none -pal=non"
    yeast$ratios <- as.matrix( read.delim( file=gzfile( "data/yeast/ratios.table.txt-GOOD.gz" ), sep="\t", as.is=T, header=T ) )
    yeast$preprocess.ratios <- function( ratios ) { ## remove rows with <10% of conditions changing
      cat( "SCE: Filtering out small change (<10% with >2-fold) rows from ratios matrix...\n" )
      nochange <- apply( 2^ratios, 1, function( i ) mean( ! is.na( i ) & ( i < 0.5 | i > 2 ), na.rm=T ) ) < 0.1
      ratios <- ratios[ ! nochange, ]
      super.preprocess.ratios( ratios ) ## Default preprocessing and mean/variance scaling
    }
    source( "cmonkey-data-load.R" )
    yeast$pp.ints <-load.sif.interactions( "data/yeast/yeast-ppInt.sif.gz" )
    yeast$bind.ints <-load.sif.interactions( "data/yeast/BIND_4932.sif.gz" )
    yeast$chip.chip <-load.sif.interactions( "data/yeast/yeast-chipChip.sif.gz" )
    yeast$net.weights <- c( string=0.25, pp.ints=0.5, bind.ints=0.1 ) ##string=0.25
    yeast$grouping.weights <- c( chip.chip=1 )
    yeast$cm.func.each.iter <- function(){if(!iter%in%plot.iters)return();fc<-clusters.w.func('proteasom');
                                       cat("HERE:",which.max(fc),max(fc,na.rm=T),"\n")}
    yeast$favorite.cluster <- function() which.max( clusters.w.func( "proteasom" ) )
    yeast <<- yeast

    cat( "Packaging Ath data...\n" ) ## arabidobsis data "ath" package goes in "bigdata" package
    ath <- list( organism="ath", rsat.species="Arabidopsis_thaliana", is.eukaryotic=T, k.clust=600,
                motif.upstream.search=c( -20, 500 ), motif.upstream.scan=c( -30, 1000 ) )
    ath$meme.cmd <- "./progs/meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 12 -mod zoops -nostatus -text -cons $none -pal=non"
    ratios <- "data/ath/AtGE_Abiotic_SR_2DRnew.txt.gz" ##AtGE_Abiotic_SR_2DR.txt.gz"
    ratios <- read.delim( file=gzfile( ratios ), sep="\t", as.is=T, header=T )
    ##attr( ratios, "rnames" ) <- toupper( ratios[ ,1 ] ); ath$ratios <- as.matrix( ratios[ ,-1 ] ); rm( ratios )
    rownames( ratios ) <- toupper( ratios[ ,1 ] ); ath$ratios <- as.matrix( ratios[ ,-1 ] ); rm( ratios )
    ath$preprocess.ratios <- function( ratios ) { ## remove rows with <50% of conditions changing
      cat( "ATH: Filtering out small change (<10% with >2-fold) rows from ratios matrix...\n" )
      nochange <- apply( 2^ratios, 1, function( i ) mean( ! is.na( i ) & ( i < 0.5 | i > 2 ), na.rm=T ) ) < 0.1
      ratios <- ratios[ ! nochange, ]
      super.preprocess.ratios( ratios ) ## Default preprocessing and mean/variance scaling
    }
    ath <<- ath
    
    source( "../generic_scripts/construct.package.R" )
    construct.package( "cMonkey.bigdata", version=version, source.files=NULL, ##source.files,
                      data=list( ecoli="ecoli", yeast="yeast", ath="ath" ), ##required="cMonkey",
                      short.desc="Additional example (big)data for use with the cMonkey integrated biclustering algorithm",
                      long.desc="Additional example (big)data for use with the cMonkey integrated biclustering algorithm" )

    if ( install ) system( sprintf( "R CMD INSTALL lib/cMonkey.bigdata_%s.tar.gz", version ) )
  }
  
  if ( update.web ) {
    system( sprintf( "cp -fv lib/index.html lib/cMonkey*_%s.tar.gz ~/Sites/cMonkey/", version ) )
    system( sprintf( "rpl VERSION \"%s\" ~/Sites/cMonkey/index.html", version ) )
    system( "cp -fv ~/Sites/cMonkey/index.html ~/Sites/cMonkey/cmonkey.html" )
    if ( install ) {
      cwd <- setwd( "~/Library/R/packages" ) ## This will change - works for pinnacle!
      ##system( sprintf( "zip -r cMonkey_%s.zip cMonkey", version ) )
      ##system( sprintf( "zip -r cMonkey.data_%s.zip cMonkey.data", version ) )
      ##if ( bigdata ) system( sprintf( "zip -r cMonkey.bigdata_%s.zip cMonkey.bigdata", version ) )
      ##system( sprintf( "mv -v cMonkey*.zip %s/lib/", cwd ) )
      setwd( cwd )
    }
    md5sums <- system( sprintf( "md5sum lib/cMonkey*_%s*", version, cwd ), intern=T )
    cat( sprintf( "VERSION %s", version ), md5sums, "\n", sep="\n", file="lib/md5sums.txt" )
    cat( version, "\n", file="lib/VERSION" )
    system( sprintf( "scp lib/VERSION lib/md5sums.txt ~/Sites/cMonkey/cmonkey.html lib/cMonkey*_%s.tar.gz lib/cMonkey*_%s.zip bragi:/local/apache2/htdocs/cmonkey/", version, version ) )
    system( sprintf( "scp lib/cMonkey_%s.tar.gz bragi:/local/apache2/htdocs/cmonkey/cMonkey_latest.tar.gz", version ) )
  }  
}

## Create an ensemble from a single run by setting, e.g.
##   cm.func.each.iter <- function() if ( iter %% 25 == 0 ) save.image( sprintf( "zzz_hpy_%04d.RData", iter ) )
## Then run this with rdata.glob="zzz_hpy_*.RData"
## This can also be used for RData output from multiple runs.
## Then, to get bicluster #24 in RData file from iter=1900, for example:
## which(names(e$fnames.to.cluster)=="./zzz_hpyu_1900.RData"&e$fnames.to.cluster==24)
cmonkey.ensemble <- function( rdata.glob ) {
  files <- list.files( patt=glob2rx( rdata.glob ), full=T )
  all.ratios <- env <- fnames.to.cluster <- NULL
  for ( f in files ) {
    print( f )
    load( f )
    if ( ! exists( "e" ) ) next
    e$clusterStack <- e$get.clusterStack( force=T )
    if ( is.null( env ) ) {
      env <- e
      all.ratios <- e$get.cluster.matrix()
      env$meme.scores <- env$meme.scores[[ 1 ]][ 1:env$k.clust ]
      tmp <- e$get.stats() ##stats[ nrow( e$stats ), ]
      env$stats <- tmp
      next
    }

    env$k.clust <- env$k.clust + e$k.clust
    rats <- e$get.cluster.matrix()
    if ( any( ! rownames( rats ) %in% rownames( all.ratios ) ) ||
        any( ! colnames( rats ) %in% colnames( all.ratios ) ) ||
        any( ! rownames( all.ratios ) %in% rownames( rats ) ) ||
        any( ! colnames( all.ratios ) %in% colnames( rats ) ) ) {
      tmp <- matrix( NA, nrow=length( unique( c( rownames( rats ), rownames( all.ratios ) ) ) ),
                    ncol=length( unique( c( colnames( rats ), colnames( all.ratios ) ) ) ) )
      rownames( tmp ) <- unique( c( rownames( rats ), rownames( all.ratios ) ) )
      colnames( tmp ) <- unique( c( colnames( rats ), colnames( all.ratios ) ) )
      tmp[ rownames( rats ), colnames( rats ) ] <- rats
      tmp[ rownames( all.ratios ), colnames( all.ratios ) ] <- all.ratios
      all.ratios <- tmp
      rm( tmp )
    }

    new.rm <- matrix( 0, nrow=nrow( all.ratios ), ncol=ncol( env$row.membership ) + ncol( e$row.membership ) )
    rownames( new.rm ) <- rownames( all.ratios )
    new.rm[ rownames( env$row.membership ), 1:ncol( env$row.membership ) ] <- env$row.membership
    tmp <- e$row.membership
    tmp[ tmp != 0 ] <- tmp[ tmp != 0 ] + max( env$row.membership )
    new.rm[ rownames( e$row.membership ), ( ncol( env$row.membership ) + 1 ):ncol( new.rm ) ] <- tmp
    env$row.membership <- new.rm; rm( new.rm, tmp )
    ## tmp <- tmp[ rownames( tmp ) %in% rownames( env$row.membership ), ,drop=F ]
    ## tmp2 <- matrix( 0, nrow=nrow( all.ratios ), ncol=ncol( tmp ) ); rownames( tmp2 ) <- rownames( all.ratios )
    ## tmp2[ rownames( tmp ), ] <- tmp
    ## env$row.membership <- cbind( env$row.membership, tmp2 )
    
    tmp <- sort( unique( e$row.membership[ e$row.membership != 0 ] ) ); names( tmp ) <- rep( f, length( tmp ) )
    fnames.to.cluster <- c( fnames.to.cluster, tmp )

    new.cm <- matrix( 0, nrow=ncol( all.ratios ), ncol=ncol( env$col.membership ) + ncol( e$col.membership ) )
    rownames( new.cm ) <- colnames( all.ratios )
    new.cm[ rownames( env$col.membership ), 1:ncol( env$col.membership ) ] <- env$col.membership
    tmp <- e$col.membership
    tmp[ tmp != 0 ] <- tmp[ tmp != 0 ] + max( env$col.membership )
    new.cm[ rownames( e$col.membership ), ( ncol( env$col.membership ) + 1 ):ncol( new.cm ) ] <- tmp
    env$col.membership <- new.cm; rm( new.cm, tmp )
    
    ## tmp <- e$col.membership
    ## tmp[ tmp != 0 ] <- tmp[ tmp != 0 ] + max( env$col.membership )
    ## tmp <- tmp[ rownames( tmp ) %in% rownames( env$col.membership ), ]
    ## tmp2 <- matrix( 0, nrow=ncol( all.ratios ), ncol=ncol( tmp ) )
    ## rownames( tmp2 ) <- colnames( all.ratios )
    ## tmp2[ rownames( tmp ), ] <- tmp
    ## env$col.membership <- cbind( env$col.membership, tmp2 )

    if ( length( e$clusterStack ) < e$k.clust ) e$clusterStack[[ e$k.clust ]] <- ""
    env$clusterStack <- c( env$clusterStack, e$clusterStack[ 1:e$k.clust ] )
    if ( is.null( e$meme.scores ) || is.null( e$meme.scores[[ 1 ]] ) )
      e$meme.scores[[ 1 ]] <- as.list( rep( "", e$k.clust ) )
    if ( length( e$meme.scores[[ 1 ]] ) < e$k.clust ) e$meme.scores[[ 1 ]][[ e$k.clust ]] <- ""
    env$meme.scores <- c( env$meme.scores, e$meme.scores[[ 1 ]][ 1:e$k.clust ] )

    tmp <- e$get.stats() ##stats[ nrow( e$stats ), ]
    if ( nrow( e$stats ) > 0 && ! all( colnames( e$stats ) %in% colnames( tmp ) ) ) {
      tmp2 <- e$stats[ 1, ] * NA
      tmp2[ , colnames( tmp ) ] <- tmp; tmp <- tmp2
    } else if ( nrow( e$stats ) <= 0 || ! all( colnames( tmp ) %in% colnames( e$stats ) ) ) {
      tmp <- tmp[ ,colnames( e$stats ) ]
    }
    env$stats <- rbind( env$stats, tmp ); rm( tmp, tmp2 )
    rm( tmp, tmp2, e )
  }
  env$ratios <- list( ratios=all.ratios )
  attr( env$ratios, "rnames" ) <- sort( unique( unlist( lapply( env$ratios, rownames ) ) ) )
  attr( env$ratios, "cnames" ) <- sort( unique( unlist( lapply( env$ratios, colnames ) ) ) )
  attr( env$ratios, "nrow" ) <- length( attr( env$ratios, "rnames" ) )
  attr( env$ratios, "ncol" ) <- length( attr( env$ratios, "cnames" ) )
  attr( env$ratios$ratios, "maxRowVar" ) <- mean( apply( env$ratios$ratios, 1, var, use="pair" ), na.rm=T ) ##* 1.2
  attr( env$ratios$ratios, "all.colVars" ) <- apply( env$ratios$ratios, 2, var, use="pair", na.rm=T )
  env$meme.scores <- list( `upstream meme`=env$meme.scores )
  env$fnames.to.cluster <- fnames.to.cluster
  rm( n.scores, net.scores, row.scores, r.scores, rr.scores, mot.scores, m.scores, row.memb, col.memb, envir=env )
  env
}

cmonkey.ensemble.analysis <- function( e, cluster.motifs=F ) { ## e from cmonkey.ensemble()
  mclapply <- lapply ## even 2 cores is too much for my mac
  ## Get a matrix of the # of times (total) that each pair of genes is in the same cluster
  all.rm <- lapply( 1:nrow( e$row.membership ), function( i ) { i <- e$row.membership[ i, ]; i[ i != 0 ] } )
  row.ov <- do.call( rbind, mclapply( all.rm, function( i ) sapply( all.rm, function( j ) sum( i %in% j ) ) ) )
  ## Matrix of the total # of (unique) clusters that genes i OR j are in (union)
  all.len <- lapply( all.rm, length )
  n.clust <- do.call( rbind, mclapply( all.len, function( i ) sapply( all.len, function( j ) min( c( i, j ) ) ) ) )
  rm( all.rm, all.len )

  row.ov[ lower.tri( row.ov, diag=T ) ] <- 0
  n.clust[ lower.tri( n.clust, diag=T ) ] <- 1
  rownames( row.ov ) <- rownames( n.clust ) <- rownames( e$row.membership )
  row.ov[ n.clust < 10 ] <- 0; n.clust[ n.clust < 10 ] <- Inf
  r.sif <- which( row.ov / n.clust > 0.2, arr=T )
  r.sif <- data.frame( g1=rownames( r.sif ), g2=rownames( row.ov )[ r.sif[ ,2 ] ], type=rep( "bic", nrow( r.sif ) ),
                      weight=row.ov[ r.sif ] / n.clust[ r.sif ] )
  r.sif <- r.sif[ order( r.sif$weight, decreasing=T ), ]; rownames( r.sif ) <- NULL

  gs <- unique( c( as.character( r.sif$g1 ), as.character( r.sif$g2 ) ) )
  r.na <- data.frame( g=gs, short=e$get.long.names( gs, short=T )[ gs ] )
  sh <- as.character( r.na$short ); sh[ which( sh == "" ) ] <- as.character( r.na$g[ which( sh == "" ) ] )
  r.na$short <- as.factor( sh ); rm( sh, gs )

  out <- list( sif=r.sif, na=r.na )

  if ( cluster.motifs ) {
    tt.out <- e$motif.similarities.tomtom()
    e$parallel.cores <- 1
    tt.out2 <- e$cluster.tomtom.results( tt.out )
    ## TODO: use mcl to do the clustering...? see cmPostProc2_mot_metaclustering2.R
    out$tt.out <- tt.out
    out$tt.out2 <- tt.out2
  }
  out
}

cmonkey.init.ec2 <- function( inst, ... ) {
  source( "~/scratch/halo/generic_scripts/ec2-utils.R" )
  ec2.setenv()
  ##inst <- ec2.spawn.spot.instances( ... )
  ##ec2.setup( inst, ec2.tools=F, full=F )
  f <- tempfile()
  cran.repos <- 'cran.cnr.Berkeley.edu'
  cat( "sudo killall -9 console-kit-daemon\n",
      "sudo mv /usr/sbin/console-kit-daemon /usr/sbin/console-kit-daemon.bkup\n",
      "sudo cp /bin/true /usr/sbin/console-kit-daemon\n",
      sprintf( "tar -xvzf %s\n", basename( options( "ec2.tools.tarball" )$ec2.tools.tarball ) ),
      "mkdir biclust; cd biclust; tar -xvzf ~/z_for_ec2.tgz; cd ~/\n",
      ##"sudo dd if=/dev/zero of=/mnt/2Gb.swap bs=1M count=2048; sudo mkswap /mnt/2Gb.swap\n",
      ##"sudo swapon /mnt/2Gb.swap; sudo sysctl vm.swappiness=10\n",
      "sudo apt-get -f -y --force-yes update\n",
      "sudo apt-get -f -y --force-yes install r-base-core r-base-dev tcsh r-recommended rpl s3cmd libcurl4-openssl-dev\n", ##r-cran-mass r-cran-vr r-cran-matrix ruby subversion libopenssl-ruby littler python-rpy pdftk openjdk-6-jdk openjdk-6-jre emacs23-nox pdftk r-cran-doSNOW heirloom-mailx sendmail
      "echo \"*       soft    nofile  1024\" | sudo tee -a /etc/security/limits.conf\n",
      "echo \"*       hard    nofile  65535\" | sudo tee -a /etc/security/limits.conf\n",
      "echo \"ulimit -n 65535\" >> ~/.bashrc\n", ## set the limit at each login
      "R CMD INSTALL lib/cMonkey_4.8.4.tar.gz\n",
      "mkdir ~/R; echo \".libPaths( c( \\\"~/R/\\\", .libPaths() ) )\" >~/.Rprofile\n",
      "echo \"require( utils ); require( graphics ); require( stats )\" >>~/.Rprofile\n",
      "echo \"my.utils <- attach( NULL, name=\\\"my.utils\\\" )\n",
      "sys.source( \\\"~/my.util.R\\\", env=my.utils )\" >>~/.Rprofile\n",
      sprintf( "Rscript -e \"Sys.setenv(MAKE=\\\"make -j 8\\\");install.packages( c(%s), dep=F, lib=\\\"~/R\\\", repos=\\\"http://%s/\\\" )\"\n", "\\\"cMonkey\\\",\\\"multicore\\\",\\\"igraph\\\",\\\"RSVGTipsDevice\\\",\\\"glmnet\\\",\\\"lars\\\"", cran.repos ),
      sprintf( "Rscript -e \"Sys.setenv(MAKE=\\\"make -j 8\\\");install.packages( c(%s), dep=T, lib=\\\"~/R\\\", repos=\\\"http://%s/\\\" )\"\n", "\\\"doMC\\\",\\\"RCurl\\\"", cran.repos ),
      "mkdir ~/progs; cd ~/progs; wget \"http://meme.nbcr.net/downloads/old_versions/meme_4.3.0.tar.gz\"\n",
      "tar -xvzf meme_4.3.0.tar.gz; cd ~/progs/meme_4.3.0; mkdir local\n",
      "./configure --prefix=`pwd`/local/ --enable-dependency-tracking --enable-opt --disable-shared --disable-fast-install --enable-serial\n",
      "make -j 8; make install; cd ~/progs/; ln -s meme_4.3.0/local/bin/meme\n",
      "ln -s meme_4.3.0/local/bin/mast; ln -s meme_4.3.0/local/bin/dust; ln -s meme_4.3.0/local/bin/tomtom\n",
      "cd ~/progs; tar -xvzf ~/weeder1.4.2.tar.gz; cd Weeder1.4.2; ./compileall\n",
      "cd ~/progs; ln -s Weeder1.4.2/weederlauncher.out; ln -s Weeder1.4.2/weederTFBS.out; ln -s Weeder1.4.2/adviser.out\n",
      sep="", file=f )
  system( "tar cvf - `find . -name '*.R'` lib/cMonkey_4.8.4.tar.gz |gzip -c >z_for_ec2.tgz" )
  ec2.upload( inst, sprintf( "z_for_ec2.tgz '%s' '%s' '%s' '%s'", options( "ec2.tools.tarball" )$ec2.tools.tarball, f,
                            '/Users/dreiss/scratch/halo/generic_scripts/my.util.R', './progs/weeder1.4.2.tar.gz' ),
             wait=T, via.s3=F )
  ec2.exec( inst, sprintf( "source %s", basename( f ) ) )
  ec2.get.from.s3( inst, file="cmonkey_data.tgz", bucket="s3://dreiss-data/cmonkey/" )
  ec2.exec( inst, "cd ~/biclust; tar -xvzf ~/cmonkey_data.tgz; ln -s ~/progs" )
  inst
}
#!endif
