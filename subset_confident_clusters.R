project.dir <- Sys.getenv( "SC_CLK856" );
source( paste0( project.dir, "/scr/clk856.config.R" ) );

set.seed( l.config[[ "seed" ]] );

library( gtools );

packageVersion( "Seurat" );

library( tidyverse );
library( cowplot );
library( magrittr );
library( ComplexHeatmap );
library( circlize );

args <- commandArgs(trailingOnly=TRUE);
if ( length( args ) > 0 ){
  dir.in        <- args[ 1 ];
  file.seurat.in  <- args[ 2 ];
}else{
  stop( "ERROR: arguments to plot.DE.R not defined" );
}

if ( ! exists( "dir.in" ) | ! exists( "file.seurat.in" ) ){
  stop( "ERROR: not all arguments to compute.DE.R were defined" );
}

setwd( dir.in );

s.in <- readRDS( file = file.seurat.in );

if( !exists("Subset_cluster_results")){
  dir.create("Subset_cluster_results",showWarnings = FALSE);
  cat("Creating Subset_cluster_results directory \n change work directory to Subset_cluster_results... \n");
  setwd("Subset_cluster_results");
}


nDims <- l.config[[ "N.PCA.INTE" ]];
v.clk.genes <- l.config[[ "clock.genes" ]];

l.confident.clusters <- c("1:DN1p_CNMa","2:s_LNv","3:DN1a","4:DN1p","5:LNd_Trissin","6:DN1p","7:DN1p","8:LN_ITP","9:LNd_NPF","12:LNd","14:DN1p","15:DN1p_CNMa","18:DN1p",
                          "19:DN2","20:DN3","25:l_LNv","29:LPN")


s.in.hcf <- subset(s.in,cells = s.in@meta.data$cell.names[s.in@active.ident %in% l.confident.clusters])

s.in.hcf <- RunPCA( s.in.hcf, assay="integrated", npcs = nDims, verbose = FALSE ); cat( "RunPCA done\n" );
cat( "Running RunTSNE\n" );
s.in.hcf <- RunTSNE( s.in.hcf, reduction = "pca", dims = 1:nDims ); cat( "RunTSNE done\n" );

cat( "PCA and TSNE ran\n" );

cat( "Finding Neighbors\n" );
s.in.hcf <- FindNeighbors(s.in.hcf, assay="integrated", reduction = "pca", dims = 1:nDims, force.recalc=T );

cat( "Finding Clusters\n" );
s.in.hcf <- FindClusters(s.in.hcf, assay="integrated", resolution = 1.0 )

cat( "table of clusters:", table( s.in.hcf@active.ident ), "\n" );
cat( "table of time: ", table( s.in.hcf@meta.data$time ), "\n" );
cat( "done with high confident cluster subseting!\n" );


s.in.hcf@meta.data$Idents <- s.in@active.ident[colnames(s.in.hcf)]

Idents(s.in.hcf) <- s.in.hcf$Idents;


saveRDS(s.in.hcf,file = "s.both.ld.dd.inte.50.subset.clusters.rds")
saveRDS(s.in.hcf,file = "../s.both.ld.dd.inte.50.subset.clusters.rds")

cat("Jobs are done. \n")




