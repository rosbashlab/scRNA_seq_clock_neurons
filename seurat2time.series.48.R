project.dir <- Sys.getenv( "SC_CLK856" );
source( paste0( project.dir, "/scr/clk856.config.R" ) );

library( gtools );

packageVersion( "Seurat" );


args <- commandArgs(trailingOnly=TRUE);
if ( length( args ) > 0 ){
    fileIn     <- args[ 1 ];
    experiment <- args[ 2 ];
    fileOut    <- args[ 3 ];
} else{
    stop( "ERROR: arguments to seurat2time.series.48.R not defined" );
}

cat( "fileIn=", fileIn, "experiment=", experiment, "\n" );


## read seurat file
s.in.orig <- readRDS( file=fileIn );
cat( "ncells in original file =", length( s.in.orig@meta.data$nCount_RNA ), "\n" );


s.in <- subset( s.in.orig, subset = experiment == experiment  );
cat( "ncells in experiment file =", length( s.in@meta.data$nCount_RNA ), "\n" );

v.times <- mixedsort( unique( s.in@meta.data$time ) );
cat( "v.times=", v.times, "\n" );

v.clusters <- as.character( mixedsort( names( which( table( s.in@active.ident ) > 0 ) ) ) );
cat( "clusters=", v.clusters, "\n" );


df.means <- data.frame();

for ( cluster in v.clusters ){
    cat( "working on cluster=", cluster, "\n" );
    df.loc <- data.frame( gene = rownames( s.in$RNA@data ), cluster = cluster, experiment = experiment );
    for ( time in v.times ) {
        cat( "experiment=", experiment, "cluster=", cluster, "time=", time, "\n" );
        v.cells.all <- s.in@meta.data$cell.names[ s.in@meta.data$time == time & as.vector( s.in@active.ident ) == cluster ];
                                        #cat( "  number of all cells = ", length( v.cells.all ), "\n" );
        v.cells.all.rand <- sample( v.cells.all, length( v.cells.all ), replace = F );
        if ( length( v.cells.all.rand ) > 1 ){ 
            for ( rep in c( "r1", "r2" ) ) {
                v.cells.loc <- c();
                if ( rep == "r1" ) {
                    v.cells.loc <- v.cells.all.rand[ 1:floor( length( v.cells.all.rand ) / 2 ) ];
                }
                else if ( rep == "r2" ) {
                    v.cells.loc <- v.cells.all.rand[ ( floor( length( v.cells.all.rand ) / 2 ) + 1 ):length( v.cells.all.rand ) ];
                }
                                        #cat( "  length( vcells.loc ) = ", length( v.cells.loc ), "\n" );
                m <- as.matrix( exp( s.in$RNA@data[ , v.cells.loc ] ) -1 );

                                        #cat( "ncol( m ) = ", ncol( m ), "\n" );
                                        #print( unique( round( apply( m, 2, sum ), digits=0) ) );
                v.mean <- apply( m, 1, mean, na.rm = T );
                                        #print( str( v.mean ) );
                name <- paste0( time, ".", rep );
                                        #df.loc[ , name ] <- v.mean;
                df.loc[ , name ] <- v.mean;
            }
        } else if ( length( v.cells.all.rand ) == 1 ){
            df.loc[ ,  paste0( time, ".r1" ) ] <- as.vector( s.in$RNA@data[ , v.cells.all.rand ] );
            df.loc[ ,  paste0( time, ".r2" ) ] <- NA;
        } else{
            df.loc[ ,  paste0( time, ".r1" ) ] <- NA;
            df.loc[ ,  paste0( time, ".r2" ) ] <- NA;
        }
    }
    df.means <- rbind( df.means, df.loc ); 
}

df.means$index <- seq( 1:nrow( df.means ) );
                                        #str( df.means );
                                        #apply( df.means[ , names( df.means ) != "gene" ], 2, sum );
v.names.ord <- c( "index", "gene", "cluster", "experiment", "zt02.r1", "zt06.r1", "zt10.r1", "zt14.r1", "zt18.r1", "zt22.r1", "zt02.r2", "zt06.r2", "zt10.r2", "zt14.r2", "zt18.r2", "zt22.r2" );
df.means <- df.means[ , v.names.ord ];

#fileOut <- paste0( "time.series.48.", experiment, ".rand.split.csv" );
cat( "writing output file ", fileOut, "\n" );
write.csv( df.means, file = fileOut, row.names = F  ); 



cat( "seurat2time.series.48.R done!\n" );


