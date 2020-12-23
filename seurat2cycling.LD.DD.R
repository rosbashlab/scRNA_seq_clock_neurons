project.dir <- Sys.getenv( "SC_CLK856" );
source( paste0( project.dir, "/scr/clk856.config.R" ) );

set.seed( l.config[[ "seed" ]] );


library( devtools );

packageVersion( "Seurat" );

library( gtools );

library( MetaCycle );

args <- commandArgs(trailingOnly=TRUE);
if ( length( args ) > 0 ){
    dir.work   <- args[ 1 ];
    fileIn     <- args[ 2 ];
    experiment <- args[ 3 ];
} else{
    stop( "ERROR: not all inputs were defined" );
}


short.experiment <- gsub(".*856_","",experiment);
cat("short.experiment is ",short.experiment,"\n\n");


cat( "dir.work=", dir.work, "\n" );
setwd( dir.work );

cat( "fileIn=", fileIn, "experiment=", experiment, "\n" );


MIN.PERIOD       <- l.config[[ "MIN.PERIOD" ]];
MAX.PERIOD       <- l.config[[ "MAX.PERIOD" ]];
MIN.CLUSTER.MEAN <- l.config[[ "MIN.CLUSTER.MEAN" ]];


l.info <- list();
l.info[[ "MIN.PERIOD" ]]       <- MIN.PERIOD;
l.info[[ "MAX.PERIOD" ]]       <- MAX.PERIOD;
l.info[[ "MIN.CLUSTER.MEAN" ]] <- MIN.CLUSTER.MEAN;
l.info[[ "fileIn" ]]           <- fileIn;
l.info[[ "experiment" ]]       <- experiment;

fileOutCore <- fileIn;
fileOutCore <- gsub( ".rds", "", fileOutCore );
fileOutCore <- gsub( ".*\\/", "", fileOutCore );
fileOutCore <- paste0( fileOutCore, ".", MIN.PERIOD, "-", MAX.PERIOD );
cat( "fileOutCore=", fileOutCore, "\n" );

exper.short <- sub( ".*_", "", experiment );
l.info[[ "exper.short" ]] <- exper.short;


## read seurat file
s.in.orig <- readRDS( file=fileIn );
cat( "ncells in original file =", length( s.in.orig@meta.data$nCount_RNA ), "\n" );


s.in <- subset( s.in.orig, subset = experiment == experiment  );
cat( "ncells in experiment file =", length( s.in@meta.data$nCount_RNA ), "\n" );
l.info[[ "n.cells" ]] <- length( s.in@meta.data$nCount_RNA );
rm( s.in.orig );

v.clusters <- as.character( mixedsort( names( which( table( s.in@active.ident ) > 0 ) ) ) );
cat( "clusters=", v.clusters, "\n" );
l.info[[ "n.clusters" ]] <- length( v.clusters );


## compute mean gene expression per cluster
m.c.mean <-
    sapply( v.clusters, function( cluster ){
        apply( exp( s.in$RNA@data[ , s.in@active.ident == cluster ] ) -1, 1, mean );
    } );
str( m.c.mean );


## data frame containing information about clusters
df.c.info <- data.frame( cluster=v.clusters, row.names=v.clusters, stringsAsFactors=FALSE );
str( df.c.info );


## compute number of genes expressed per cluster 
v.n.expressed <- apply( m.c.mean, 2, function( x ){ sum( x > 0, na.rm=T ); } );
df.c.info[ names( v.n.expressed ), paste0( "n.g.expressed", "_", exper.short) ] <- v.n.expressed;
str( df.c.info );

## compute number of genes above MIN.CLUSTER.MEAN
v.n.mean.th <- apply( m.c.mean, 2, function( x ){ sum( x > MIN.CLUSTER.MEAN, na.rm=T ); } );
str( v.n.mean.th );
df.c.info[ names( v.n.mean.th ),  paste0( "n.g.mean.th", "_", exper.short) ]  <- v.n.mean.th;
str( df.c.info );


write.csv( df.c.info, paste0( fileOutCore, ".cluster.info.csv" ) );

l.info[[ "df.c.info" ]] <- df.c.info;
str( l.info );


## find genes to be kept and used for cycling computations
v.genes.keep <- sort( names( which( apply( m.c.mean > MIN.CLUSTER.MEAN, 1, any ) ) ) );
cat( "length( v.genes.keep ) = ", length( v.genes.keep ), "\n" );
l.info[[ "n.genes.mean.th" ]] <- length( v.genes.keep );
print( l.info );
saveRDS( l.info, paste0( fileOutCore, ".l.info.rds" ) );
cat( file=paste0( fileOutCore, ".v.genes.keep.forcycling.txt" ), v.genes.keep, sep="\n" );


## compute mean expression in clusters per time points
print( date() );
a.c.t <- f.computeTimeMeansSem( s.in, v.clusters, v.genes.keep );
print( date() );
a.c.t[ "tim", , , "mean" ];
a.c.t[ "Clk", , , "mean" ];


###### Run Meta
cat( "computing cycling with meta\n" );
date();
nTimePoints <- l.config[[ "N.TIME.POINTS" ]];
nCoresMeta  <- l.config[[ "N.CORES.META" ]];
l.meta <- f.seurat2meta( s.in, rev( v.clusters ), v.genes.keep, MIN.PERIOD, MAX.PERIOD, nTimePoints, Parallelize = TRUE, nCoresMeta );
#l.meta <- f.seurat2meta( s.in, as.character( seq( 10, 33, by=1) ), v.genes.keep, MIN.PERIOD, MAX.PERIOD, nTimePoints, Parallelize=TRUE, nCoresMeta );
saveRDS( l.meta, file=paste0( fileOutCore, ".l.meta.", "rds" ) );
date();



## process cycling data
cat( "processing cycling data\n" );
                                        #df.meta <- l.meta[[ "cb" ]];
df.meta.rb <- l.meta[[ "rb" ]];
names( df.meta.rb )[ names( df.meta.rb ) == "CycID" ] <- "gene";
str( df.meta.rb );

df.c.res <- data.table::melt( m.c.mean, varnames=c( "gene", "cluster" ), value.name="mean" );
df.c.res$gene <- as.vector( df.c.res$gene );
df.c.res$cluster <- as.character( df.c.res$cluster );
str( df.c.res );


df.meta.rb$gene.cluster <- paste( df.meta.rb$gene, df.meta.rb$cluster, sep=":" );
df.c.res$gene.cluster <- paste( df.c.res$gene, df.c.res$cluster, sep=":" );

v.common.gene.cluster <- intersect( df.c.res$gene.cluster, df.meta.rb$gene.cluster );
df.meat.rb <- df.meta.rb[ df.meta.rb$gene.cluster %in% v.common.gene.cluster, ];
df.c.res <- df.c.res[ df.c.res$gene.cluster %in% v.common.gene.cluster, ];


str( df.meta.rb );
df.c.res <- merge( df.c.res, df.meta.rb, by=c( "gene", "cluster", "gene.cluster" ), all=TRUE )
str( df.c.res );
rm( df.meta.rb );



## compute ts

cat( "Computing df.ts\n" );
df.mean <- data.table::melt( a.c.t[ , , , "mean"], id.vars=c( "gene", "cluster", "time" ), measure.vars=c( "mean" ), value.name="mean" );
df.sem <- data.table::melt( a.c.t[ , , , "sem"], id.vars=c( "gene", "cluster", "time" ), measure.vars=c( "sem" ), value.name="sem" );
df.ts <- merge( df.mean, df.sem, by=c( "gene", "cluster", "time" ) );
rm( df.mean, df.sem );
write.csv( df.ts, file=paste0( fileOutCore, ".df.ts.csv") );

## compute amplitudes
cat( "computing amplitude\n" );
print( date() );
m.amp <- t( apply( df.c.res, 1, function( x ){
    #cat( x, "\n", sep=":" );
    gene <- as.character( x[ "gene" ] );
    cluster <- as.character( x[ "cluster" ] );
    #cat( "gene=", gene, " cluster=", cluster, "\n", sep="" );
    v.m <- a.c.t[ gene, cluster, , "mean" ];
    v.s <- a.c.t[ gene, cluster, , "sem" ];
    
    mx <- max( v.m, na.rm=T );
    mn <- min( v.m, na.rm=T );
    avg <- mean( v.m, na.rm=T );
    med <- median( v.m, na.rm=T );
    
    amp2 <- 0;
    if ( mn > 0 ){ amp2 <- mx / mn; }
    else if ( mx > 0 ){ amp2 <- 1000; }
    
    amp3 <- 0;
    if ( mn + mx > 0 ){ amp3 <- ( mx - mn ) / ( mx + mn ); }
    return( c( "AMP.MXMN" = amp2, "AMP.RC" = amp3, "MX" = mx, "MN" = mn, "MEAN" = avg, "MEDIAN" = med ) );
}
) );
print( date() );

str( m.amp );
df.c.res <- cbind( df.c.res, m.amp );
str( df.c.res );

rm( m.amp );

cat( "formatting df.c.res\n" );
df.c.res$meta2d_pvalue <- signif( df.c.res$meta2d_pvalue, digits=3 );
df.c.res$meta2d_BH.Q <- signif( df.c.res$meta2d_BH.Q, digits=3 );
colnames( df.c.res )[ colnames( df.c.res ) == "AMP" ] <- "AMP.META";
df.c.res$meta2d_AMP <- signif( df.c.res$meta2d_AMP, digits=3);
head( df.c.res );
str( df.c.res );


cat( "Writing long and wide files\n" );
df.c.res$isCycling <- df.c.res$meta2d_BH.Q < 0.05 & df.c.res$AMP.MXMN > 0;
str( df.c.res );
write.csv( df.c.res, file=paste0( fileOutCore, ".cycling.long.META.csv" ), quote=F );

v.genes.cycling <- mixedsort( unique( as.vector( df.c.res$gene[ df.c.res$isCycling ] ) ) );
str( v.genes.cycling );



v.clusters <- mixedsort( unique( df.c.res$cluster ) );
df.meta.wide <- data.frame( gene = mixedsort( unique( df.c.res$gene ) ) );                           
for ( cluster in v.clusters ){
    df.loc <- df.c.res[ df.c.res$cluster == cluster, ];
    colnames( df.loc )[ ! colnames( df.loc ) %in% c( "gene", "cluster" ) ] <- paste( colnames( df.loc )[ ! colnames( df.loc ) %in% c( "gene", "cluster" ) ], cluster, sep=":" );
    df.loc$cluster <- NULL;
    df.meta.wide <- merge( df.meta.wide, df.loc, by="gene", all=T );
}
str( df.meta.wide );

df.meta.form <- df.meta.wide;
rownames( df.meta.form ) <- NULL;
df.meta.form <- df.meta.form[ , c( "gene", mixedsort( setdiff( colnames( df.meta.form ), "gene") ) ) ];
colnames( df.meta.form );
write.csv( df.meta.form, file=paste0( fileOutCore, ".cycling.wide.META.csv" ), quote=F );

sum( df.c.res$meta2d_pvalue < 0.05, na.rm = T );
sum( df.c.res$meta2d_BH.Q < 0.05, na.rm = T );


n.cyc <- length( v.genes.cycling );
cat( "number of cycling genes:", n.cyc, "\n" );

sessionInfo();
sink( "sessionInfo-integrate.LD.txt" );
sessionInfo();
sink();


cat( "seurat2cycling.R done!\n" );


