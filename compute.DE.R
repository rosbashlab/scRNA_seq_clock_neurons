project.dir <- Sys.getenv( "SC_CLK856" );
source( paste0( project.dir, "/scr/clk856.config.R" ) );

set.seed( l.config[[ "seed" ]] );

library( gtools );

packageVersion( "Seurat" );


v.args <- commandArgs( trailingOnly = TRUE );
print( "args:\n" );
print( v.args );


if ( length( v.args ) > 0 ){
    file.seurat.in   <- v.args[ 1 ];
    file.markers.out <- v.args[ 2 ];
    cat( "file.seurat.in=", file.seurat.in, "file.markers.out=", file.markers.out, "\n" );
}else{
    file.seurat.in   <- l.config[[ "file.seurat.inte.both" ]];
    file.markers.out <- l.config[[ "file.markers.inte.LD" ]];
    print( "arguments to compute.DE.R not defined, using default file locations" );
    #stop( "ERROR: arguments to compute.DE.R not defined" );
}

cat( file.seurat.in, file.markers.out, "\n" );

if ( ! exists( "file.seurat.in" ) | ! exists( "file.markers.out" ) ){
  stop( "ERROR: not all arguments to compute.DE.R were defined" );
}



cat( "reading", file.seurat.in, "\n" );
s.in <- readRDS( file = file.seurat.in );

library( future );
plan( "multiprocess", workers = 4 );
plan();


cat( "computing all markers\n" );
date();
df.m <- FindAllMarkers( s.in, test.use="negbinom", only.pos=FALSE, assay="RNA", latent.vars=c("time.date"), verbost=TRUE  );
#df.m <- FindAllMarkers( s.in, test.use="negbinom", only.pos=TRUE, min.pct = 0.5, logfc.threshold = 1.0, assay="RNA", latent.vars=c("time.date"), verboset=TRUE  );
date();

## compute mean per cluster
m.c.mean <- sapply( mixedsort( levels( s.in@active.ident ) ), function( cluster ){ apply( exp( s.in$RNA@data[ , s.in@active.ident == cluster ] ) -1, 1, mean );  } );

df.m$TP10K.mean <- NA;
for ( ir in 1:nrow( df.m ) ){
    df.m$TP10K.mean[ ir ] <- m.c.mean[ as.character( df.m$gene[ ir ] ), as.character( df.m$cluster[ ir ] ) ];
}
df.m <- df.m[ , c( "gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2", "TP10K.mean" ) ];

if ( 0 ){
    df.out <- df.m;
    df.out <- format( df.m, scientific=TRUE, digits=3 );
    rownames( df.out ) <- NULL;
}
write.csv( df.m, file = file.markers.out );

sessionInfo();

cat( "compute.DE.R done!\n" );

