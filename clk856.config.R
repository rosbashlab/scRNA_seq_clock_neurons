## set seed

l.config                       <- list();
l.config[[ "seed" ]]           <- 1313;
l.config[[ "seed.read.data" ]] <- 138;


set.seed( l.config[[ "seed" ]] );


## directories
l.config[[ "dir.project" ]]    <- Sys.getenv( "SC_CLK856" );
l.config[[ "dir.analyses" ]]   <- paste0( l.config[[ "dir.project" ]], "analyses/" );
l.config[[ "dir.scripts" ]]    <- paste0( l.config[[ "dir.project" ]], "scr/" );

l.config[[ "dir.dbs" ]]        <- paste0( l.config[[ "dir.project" ]], "other_data/dbs/" );
l.config[[ "dir.gene.sets" ]]  <- paste0( l.config[[ "dir.project" ]], "other_data/gene_sets/" );


l.config[[ "dir.parsed.data" ]]      <- paste0( l.config[[ "dir.project" ]], "data_sc_parsed/" );

l.config[[ "dir.integrated.LD.DD" ]] <- paste0( l.config[[ "dir.analyses" ]], "integrated_LD_DD/" );


## functions
l.config[[ "file.functions" ]] <- paste0( l.config[[ "dir.scripts" ]] , "functions.main.R" );
source( l.config[[ "file.functions" ]] );

library( Seurat, lib.loc= paste0( l.config[[ "dir.project" ]], "/scr/R_loc/3.6/Resources/library_seurat_3.0.2" ) );


l.config[[ "dgeinfo.in" ]] <- list();
l.config[[ "dgeinfo.in" ]][[ "data" ]] <- list();
l.config[[ "dgeinfo.in" ]][[ "data" ]][[ "CLK856_LD" ]]  <- paste0( l.config[[ "dir.project" ]], "data/clk856_sc_ld/clk856.sc.ld.dgecounts.info.csv" );
l.config[[ "dgeinfo.in" ]][[ "data" ]][[ "CLK856_DD" ]]  <- paste0( l.config[[ "dir.project" ]], "data/clk856_sc_dd/clk856.sc.dd.dgecounts.info.new.csv" );

l.config[[ "dgeinfo.in" ]][[ "empty" ]] <- list();
l.config[[ "dgeinfo.in" ]][[ "empty" ]][[ "CLK856_DD" ]]  <- paste0( l.config[[ "dir.project" ]], "data/clk856_sc_dd//negative.controls.CLK856_DD.csv" );

l.config[[ "dgeinfo.in" ]][[ "empty" ]][[ "CLK856_LD" ]] <- "AGACTC,AGCTAG";

l.config[[ "dgeinfo.in" ]][[ "dirs" ]] <- list();
l.config[[ "dgeinfo.in" ]][[ "dirs" ]][[ "CLK856_LD" ]]  <- paste0( l.config[[ "dir.project" ]], "data/clk856_sc_ld/" );
l.config[[ "dgeinfo.in" ]][[ "dirs" ]][[ "CLK856_DD" ]]  <- paste0( l.config[[ "dir.project" ]], "data/clk856_sc_dd/" );


l.config[[ "gene2smbl.file" ]]         <- paste0( l.config[[ "dir.dbs" ]], "ensemblToSymbol.csv" ); 



## original seruat files
l.config[[ "file.seurat.orig.LD" ]] <- paste0( l.config[[ "dir.parsed.data" ]], "s.dat.in.CLK856_LD.new.ex.rds" );
l.config[[ "file.seurat.orig.DD" ]] <- paste0( l.config[[ "dir.parsed.data" ]], "s.dat.in.CLK856_DD.new.ex.rds" );



## seurat files with clustering
l.config[[ "file.seurat.inte.LD" ]] <- paste0( l.config[[ "dir.integrated.LD.DD" ]], "s.both.inte.LD.rds" );
l.config[[ "file.seurat.inte.DD" ]] <- paste0( l.config[[ "dir.integrated.LD.DD" ]], "s.both.inte.DD.rds" );
l.config[[ "file.seurat.inte.both" ]] <- paste0( l.config[[ "dir.integrated.LD.DD" ]], "s.both.inte.rds" );

## DE marker files
l.config[[ "file.markers.inte.LD" ]] <- paste0( l.config[[ "dir.integrated.LD.DD" ]], "LD_DD.markers.negbinom.csv" );


##gene sets
l.config[[ "gene2smbl.file" ]]         <- paste0( l.config[[ "dir.dbs" ]], "ensemblToSymbol.csv" ); 


## parameters
l.config[[ "SCALE.FACTOR" ]]     <- 10000;
l.config[[ "EXON.ONLY" ]]        <- TRUE;

l.config[[ "MIN.nCount.RNA" ]]   <-  6000;
l.config[[ "MAX.nCount.RNA" ]]   <- 75000;

l.config[[ "MIN.nFeature.RNA" ]] <- 1000;
l.config[[ "MAX.nFeature.RNA" ]] <- 6000;

l.config[[ "MIN.ENTROPY" ]]      <- 5.5;


l.config[[ "N.PCA.INTE" ]]       <- 50;

l.config[[ "N.PCA.PROJ" ]]       <- 30;


## cycling parameters
l.config[[ "MIN.PERIOD" ]]        <- 24;
l.config[[ "MAX.PERIOD" ]]        <- 24;
l.config[[ "MIN.CLUSTER.MEAN" ]]  <- 0.5;
l.config[[ "N.TIME.POINTS" ]]     <-  6;
l.config[[ "N.CORES.META" ]]      <-  6;


## colors
l.config[[ "v.rand.colors" ]] <- sample( rainbow( 40 ), 40, replace=F );
