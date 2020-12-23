project.dir <- Sys.getenv( "SC_CLK856" );
source( paste0( project.dir, "/scr/clk856.config.R" ) );


set.seed( l.config[[ "seed" ]] );


library( devtools );


packageVersion( "Seurat" );

library( gtools );
library( Hmisc );


dir.work <- l.config[[ "dir.integrated.LD.DD" ]];
cat( "dir.work=", dir.work, "\n" );
dir.create( dir.work, recursive=TRUE, showWarnings=TRUE );
setwd( dir.work );


dir.dbs <- l.config[[ "dir.dbs" ]];
dir.gene.sets <-  l.config[[ "dir.gene.sets" ]];

v.rand.colors <- l.config[[ "v.rand.colors" ]];


#s.dat.in.LD <- readRDS( l.config[[ "file.seurat.orig.LD" ]] );
s.dat.in.LD <- readRDS( "/Users/dariusz/clk856_elife/data_sc_parsed/fromDylan/s.dat.in.CLK856_LD.new.ex.rds" );
unique( s.dat.in.LD@meta.data$experiment );

#s.dat.in.DD <- readRDS( l.config[[ "file.seurat.orig.DD" ]] );
s.dat.in.DD <- readRDS( "/Users/dariusz/clk856_elife/data_sc_parsed/fromDylan/s.dat.in.CLK856_DD.new.ex.rds" );
unique( s.dat.in.DD@meta.data$experiment );


s.dat.in <- merge( x = s.dat.in.LD, y = s.dat.in.DD );
s.dat.in@meta.data$time.date.exper <- paste( s.dat.in@meta.data$time.date, s.dat.in@meta.data$experiment, sep = ":" );
table( s.dat.in@meta.data$time.date.exper );

s.dat.in@meta.data$time.exper <- paste( s.dat.in@meta.data$time, s.dat.in@meta.data$experiment, sep = ":" );
table( s.dat.in@meta.data$time.exper );



s.dat.in <- subset( s.dat.in, subset = experiment %in% c( "CLK856_LD", "CLK856_DD" ) );
table( s.dat.in@meta.data$time.date.exper );

s.dat.in@meta.data$log2nCount_RNA <- log2( s.dat.in@meta.data$nCount_RNA + 1 );
s.dat.in@meta.data$log2nFeature_RNA <- log2( s.dat.in@meta.data$nFeature_RNA + 1 );
s.dat.in@meta.data$entropy <- f.entropy( s.dat.in$RNA@counts );


## structure containing information about data
l.info <- list();

## compute genes not to be used in clustering
l.info[[ "geneSets" ]][[ "mt" ]]      <- v.mt.genes <- grep( "^mt:", rownames( s.dat.in@assays$RNA@data ), value=T, ignore.case=T );
l.info[[ "geneSets" ]][[ "rpls" ]]    <- v.ribo.genes <- grep( "^rp[ls]", rownames( s.dat.in@assays$RNA@data ), value=T, ignore.case=T );
l.info[[ "geneSets" ]][[ "rRNA" ]]    <- v.rRNA.genes <- grep( "rRNA", rownames( s.dat.in@assays$RNA@data ), value=T, ignore.case=T );
l.info[[ "geneSets" ]][[ "tRNA" ]]    <- v.tRNA.genes <- grep( "tRNA", rownames( s.dat.in@assays$RNA@data ), value=T, ignore.case=T );
l.info[[ "geneSets" ]][[ "ERCC" ]]    <- v.ercc.genes <- grep( "^ERCC", rownames( s.dat.in@assays$RNA@data ), value=T, ignore.case=F);
l.info[[ "geneSets" ]][[ "exclude" ]] <- v.genes.exclude <- mixedsort( unique( c( v.mt.genes, v.ribo.genes, v.rRNA.genes, v.tRNA.genes, v.ercc.genes, "EGFP" ) ) );

cat( "length( v.genes.exclude ) = ", length( v.genes.exclude ), "\n" );
write( v.genes.exclude, file="genes.exclude.txt" );

v.experiments <- unique( s.dat.in@meta.data$experiment );
cat( "v.experiments:", v.experiments, "\n" );
for ( experiment in v.experiments ){
    l.info[[ "nCells_condition" ]][[ experiment ]]  <- sum( s.dat.in@meta.data$experiment == experiment );
}
print( l.info );

for ( experiment in mixedsort( unique( s.dat.in@meta.data$experiment ) ) ){
    for ( zt in mixedsort( unique( s.dat.in@meta.data$time ) ) ){
        n.tmp <- sum( s.dat.in@meta.data$experiment == experiment & s.dat.in@meta.data$time == zt );
        if ( n.tmp > 0 ){
            l.info[[ "nCells_experiment_time" ]][[ paste(  experiment, zt, sep="_" ) ]] <- n.tmp;
        }
    }
}
print( l.info );

for ( experiment in mixedsort( unique( s.dat.in@meta.data$experiment ) ) ){
    for ( ar in mixedsort( unique( s.dat.in@meta.data$ar ) ) ){
        for ( date in mixedsort( unique( s.dat.in@meta.data$date ) ) ){
            n.tmp <- sum( s.dat.in@meta.data$experiment == experiment & s.dat.in@meta.data$ar == ar & s.dat.in@meta.data$date == date );
            if ( n.tmp > 0 ){
                l.info[[ "nCells_experiment_plate_date" ]][[ paste( experiment, ar, date, sep="_" ) ]] <- n.tmp;
            }
        }
    }
}
print( l.info );



v.times <- unique( mixedsort( s.dat.in@meta.data$time ) );
cat( "v.times: ", v.times, "\n" );
v.plates <- unique( mixedsort( s.dat.in@meta.data$ar ) );
cat( "v.plates: ", v.plates, "\n" );
v.dates <- unique( mixedsort( s.dat.in@meta.data$date ) );
cat( "v.dates: ", v.dates, "\n" );

l.info[[ "date2times" ]] <- sapply( v.dates, function( x ){ mixedsort( unique( s.dat.in@meta.data$time[ s.dat.in@meta.data$date == x ] ) ) } );
saveRDS( l.info, "l.info.rds" );

 

########### DO INTEGRATION
s.both <- subset( s.dat.in, cells = s.dat.in@meta.data$cell.names[
                                   s.dat.in@meta.data$experiment %in% c( "CLK856_LD", "CLK856_DD" ) &
                                   s.dat.in@meta.data$nCount_RNA   > l.config[[ "MIN.nCount.RNA" ]]  &
                                   s.dat.in@meta.data$nCount_RNA   < l.config[[ "MAX.nCount.RNA" ]] &
                                   s.dat.in@meta.data$nFeature_RNA > l.config[[ "MIN.nFeature.RNA" ]] &
                                   s.dat.in@meta.data$nFeature_RNA < l.config[[ "MAX.nFeature.RNA" ]] &
                                   s.dat.in@meta.data$entropy      > l.config[[ "MIN.ENTROPY" ]]
                                   ] );

l.info[[ "n.cells.filt" ]] <- length( s.both@meta.data$cell.names );
cat( "ncells in s.both=", l.info[[ "n.cells.filt" ]], "\n" );
cat( "experiments in s.both=", unique( s.both@meta.data$experiment ), "\n" );


s.both <- PercentageFeatureSet( s.both, pattern = "^mt:", col.name = "percent.mito" );
s.both@meta.data$time2 <- as.numeric( gsub( "zt", "", s.both@meta.data$time ) );



date();
nDims <- l.config[[ "N.PCA.INTE" ]];
s.both.inte <- f.cluster.seuratv3( s.both, nDims, l.info[[ "geneSets" ]][[ "exclude" ]], split.by.var = "time.exper" );
date();



# rename clusters

l.c2new <- list();
l.c2new[[ "0" ]] <- "0";
l.c2new[[ "1" ]] <- "1:DN1p_CNMa";
l.c2new[[ "2" ]] <- "2:s_LNv";
l.c2new[[ "3" ]] <- "3:DN1a";
l.c2new[[ "4" ]] <- "4:DN1p";
l.c2new[[ "5" ]] <- "5:LNd_Trissin";
l.c2new[[ "6" ]] <- "6:DN1p";
l.c2new[[ "7" ]] <- "7:DN1p";
l.c2new[[ "8" ]] <- "8:LN_ITP";
l.c2new[[ "9" ]] <- "9:LNd_NPF";
l.c2new[[ "10" ]] <- "10";
l.c2new[[ "11" ]] <- "11";
l.c2new[[ "12" ]] <- "12:LNd";
l.c2new[[ "13" ]] <- "13";
l.c2new[[ "14" ]] <- "14:DN3";
l.c2new[[ "15" ]] <- "15:DN1p_CNMa";
l.c2new[[ "16" ]] <- "16";
l.c2new[[ "17" ]] <- "17";
l.c2new[[ "18" ]] <- "18:DN1p";
l.c2new[[ "19" ]] <- "19:DN2";
l.c2new[[ "20" ]] <- "20:DN3";
l.c2new[[ "21" ]] <- "21";
l.c2new[[ "22" ]] <- "22";
l.c2new[[ "23" ]] <- "23";
l.c2new[[ "24" ]] <- "24";
l.c2new[[ "25" ]] <- "25:l_LNv";
l.c2new[[ "26" ]] <- "26";
l.c2new[[ "27" ]] <- "27";
l.c2new[[ "28" ]] <- "28";
l.c2new[[ "29" ]] <- "29:LPN";
l.c2new[[ "30" ]] <- "30";
l.c2new[[ "31" ]] <- "31";
l.c2new[[ "32" ]] <- "32";
l.c2new[[ "33" ]] <- "33";
l.c2new[[ "34" ]] <- "34";
l.c2new[[ "35" ]] <- "35";
l.c2new[[ "36" ]] <- "36";
l.c2new[[ "37" ]] <- "37";
l.c2new[[ "38" ]] <- "38";

s.tmp <- s.both.inte;
for ( c in levels( s.tmp@active.ident ) ){
  v.cells <- WhichCells( s.tmp, idents=c(c) );
  s.tmp <- SetIdent( s.tmp, cells=v.cells, value=l.c2new[[ c ]] );
}
levels( s.tmp@active.ident );
s.both.inte <- s.tmp;
rm( s.tmp );


s.both <- s.both.inte;

s.both@meta.data$Repeats <- " ";

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt02"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt02"), ]$date == "20181231","LD_1","LD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt06"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt06"), ]$date == "20190228","LD_1","LD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt10"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt10"), ]$date == "20190228","LD_1","LD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt14"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt14"), ]$date == "20181231","LD_1","LD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt18"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt18"), ]$date == "20190301","LD_1","LD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt22"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_LD" & s.both@meta.data$time == "zt22"), ]$date == "20190219","LD_1","LD_2")

# DD repeats information

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt02"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt02"), ]$date == "20190528","DD_1","DD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt06"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt06"), ]$date == "20190710","DD_1","DD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt10"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt10"), ]$date == "20190704","DD_1","DD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt14"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt14"), ]$date == "20190001","DD_1","DD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt18"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt18"), ]$date == "20190613","DD_1","DD_2")

s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt22"), ]$Repeats <- 
  ifelse(s.both@meta.data[(s.both@meta.data$experiment == "CLK856_DD" & s.both@meta.data$time == "zt22"), ]$date == "20190613","DD_1","DD_2")

s.both.inte <- s.both;
rm(s.both);


saveRDS( s.both.inte, file = l.config[[ "file.seurat.inte.both" ]] );


pdf( "elbow.plot.pdf" );
ElbowPlot( s.both.inte, ndims=nDims ); 
dev.off();
table( s.both.inte@meta.data$seurat_clusters );

write( s.both.inte$integrated@var.features, file= paste0( "LD.inte.var.features.nd.", nDims, ".txt" ) );

pdf( paste0( "LD.DD.inte.clusters.tsne.nd.", nDims, ".pdf" ) );
DimPlot(s.both.inte, reduction = "tsne", assay= "integrated", label=TRUE, label.size=4, repel=TRUE ) + theme(legend.position = "none");
dev.off();

pdf( paste0( "LD.DD.inte.clusters.Repeats.tsne.nd.", nDims, ".pdf" ) );
DimPlot(s.both.inte, reduction = "tsne", assay= "integrated", label = F,group.by = "Repeats", label.size=4, repel=TRUE );
dev.off();

pdf( paste0( "LD.DD.inte.time.tsne.nd.", nDims, ".pdf" ) );
DimPlot(s.both.inte, reduction = "tsne", assay= "integrated", label=F, group.by="time" );
dev.off();

pdf( paste0( "LD.DD.inte.experiment.tsne.nd.", nDims, ".pdf" ) );
DimPlot(s.both.inte, reduction = "tsne", assay= "integrated", label=F, group.by="experiment" );
dev.off();

pdf( paste0( "LD.DD.inte.nCount_RNA.tsne.nd.", nDims, ".pdf" ) );
FeaturePlot( s.both.inte, reduction="tsne", features=c( "nCount_RNA" ), min.cutoff="q1", max.cutoff="q99", order=T, pt.size=1  );
dev.off();

pdf( paste0( "LD.DD.inte.frac.mito.tsne.nd.", nDims, ".pdf" ) );
FeaturePlot( s.both.inte, reduction="tsne", features=c( "frac.mito" ), min.cutoff="q1", max.cutoff="q99", order=T, pt.size=1  );
dev.off();


# split the rds to LD set and DD set

dir.create("JTK_F24")

s.both.inte.ld <- subset( s.both.inte, cells=s.both.inte@meta.data[s.both.inte@meta.data$experiment == "CLK856_LD",]$cell.names) 
s.both.inte.dd <- subset( s.both.inte, cells=s.both.inte@meta.data[s.both.inte@meta.data$experiment == "CLK856_DD",]$cell.names) 

saveRDS(s.both.inte.ld,file = "s.both.inte.LD.rds");
saveRDS(s.both.inte.dd,file = "s.both.inte.DD.rds");


#system( "ln -s s.both.inte.LD.rds JTK_F24/s.both.inte.LD.rds" );
#system( "ln -s s.both.inte.DD.rds JTK_F24/s.both.inte.DD.rds" );


sessionInfo();
sink( "sessionInfo-integrate.LD.txt" );
sessionInfo();
sink();

cat( "\n", "integrate.LD.DD.R done!\n" );
