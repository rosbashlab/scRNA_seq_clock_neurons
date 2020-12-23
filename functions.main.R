

##### cluster seurat using integration
f.cluster.seuratv3 <- function( s.in, nDims, v.var.genes.exclude, k.anchor.use = 5, k.filter.use = 199, k.score.use = 30, k.weight.use = 100, split.by.var = "time" ){
                                        
    cat( "f.cluster.seurat.v3:\n" );
    cat( "k.anchor.use=", k.anchor.use, "k.filter.use=", k.filter.use, "k.score.use=", k.score.use, "k.weight.use=", k.weight.use, "split.by.var=", split.by.var, "\n" );
    
    cat( "splitting data by ", split.by.var, "\n" );
    l.in <- SplitObject( s.in, split.by = split.by.var );
    cat( "split into", length( l.in ), "objects\n" );

    cat( "running SCTransform\n" );
    for (i in 1:length( l.in ) ) {
        #l.in[[i]] <- SCTransform( l.in[[i]], vars.to.regress = c( "date", "nCount_RNA", "nFeature_RNA", "percent.mito", "dissoc.score"  ), assay="RNA", verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = NULL );
        l.in[[i]] <- SCTransform( l.in[[i]], vars.to.regress = c( "date", "nCount_RNA", "nFeature_RNA", "percent.mito"  ), assay="RNA", verbose = TRUE, return.only.var.genes = FALSE, variable.features.n = 3000 );
        l.in[[i]]@assays$SCT@var.features <- setdiff( l.in[[i]]@assays$SCT@var.features, v.var.genes.exclude );
    }

    #v.inte.features <- SelectIntegrationFeatures( object.list = l.in, nfeatures = 3000 );
    #cat( "length( v.inte.features ) = ", length( v.inte.features ), "\n" );
    
    
                                        # find features/genes common to all
    v.common <- rownames( l.in[[1]]$SCT@data );
    for (  i in 2:length( l.in ) ){
        v.common <- intersect( v.common, rownames( l.in[[i]]$SCT@data ) );
    } 
    length( v.common );
    cat( "length(v.common)=", length(v.common), "\n" );

    v.inte.features <- vector();
    if ( 1 ){                                     # find variable features/genes common to all
        v.var.common <- l.in[[1]]@assays$SCT@var.features;
        for (  i in 2:length( l.in ) ){
            v.var.common <- c( v.var.common, l.in[[i]]@assays$SCT@var.features );
        }
        length( unique( v.var.common ) );
       v.var.common <- names( which( table( v.var.common ) == length( l.in ) ) );
        # v.var.common <- names( which( table( v.var.common ) == 8 ) ); # Change this just for integration LD and DD
        
        v.var.common <- setdiff( v.var.common, v.var.genes.exclude );
        length( v.var.common );
        cat( "length(v.var.common)=", length(v.var.common), "\n" );
        v.inte.features <- v.var.common;
    } else {
        v.inte.features <- SelectIntegrationFeatures( object.list = l.in, nfeatures = 3000 );
    }
    
    cat( "length( v.inte.features ) = ", length( v.inte.features ), "\n" );
    
    #cat( "running PrepSCTIntegration\n" );
    #l.in <- PrepSCTIntegration( l.in, assay = rep( "SCT", length( l.in ) ), anchor.features = v.inte.features, sct.clip.range = NULL, verbose = TRUE );

    print( date() );
    cat( "running FindIntegrationAnchors\n" );
    cat( "k.anchor.use=", k.anchor.use, "k.filter.use=", k.filter.use, "k.score.use=", k.score.use, "\n" );
    l.in.anchors <- FindIntegrationAnchors( object.list = l.in, dims = 1:nDims, assay=rep( "SCT", length( l.in )), anchor.features=v.inte.features, k.anchor=k.anchor.use, k.filter=k.filter.use, k.score=k.score.use, verbose=T );
    #l.in.anchors <- FindIntegrationAnchors( object.list = l.in, dims = 1:nDims, assay=rep( "SCT", length( l.in )), k.anchor=k.anchor.use, k.filter=k.filter.use, k.score=k.score.use, verbose=F );
    print( date() );
    cat( "integrating data:\n\n" );
    cat( "k.weight.use=", k.weight.use, "\n" );

    cat( "running IntegrateData\n" );
    s.in.inte <- IntegrateData( anchorset = l.in.anchors, dims = 1:nDims, features.to.integrate=v.common, k.weight = k.weight.use, verbose = T );
    cat( "Data integrated\n" );
    print( date() );
    VariableFeatures( s.in.inte, assay="integrated" ) <- v.inte.features;
    cat( "length(s.in.inte$integrated@var.features)=", length(s.in.inte$integrated@var.features), "\n" );

    DefaultAssay( s.in.inte ) <- "integrated";
    s.in.inte <- ScaleData( s.in.inte, verbose=TRUE, assay="integrated", features=rownames( s.in.inte$integrated@data ) );
    cat( "Data scaled\n" );
    #cat( "Finding variable features\n" );
    #s.in.inte <- FindVariableFeatures( s.in.inte, assay="integrated", nfeatures = 2000 );
    #cat( "Variable features found\n" );
    
    
    cat( "Running PCA\n" );
    s.in.inte <- RunPCA( s.in.inte, assay="integrated", npcs = nDims, verbose = FALSE ); cat( "RunPCA done\n" );
    cat( "Running RunTSNE\n" );
    s.in.inte <- RunTSNE( s.in.inte, reduction = "pca", dims = 1:nDims ); cat( "RunTSNE done\n" );
    #cat( "Running UMAP\n" );
    #s.in.inte <- RunUMAP( s.in.inte, assay="integrated", reduction = "pca", dims = 1:nDims ); cat( "RunUMAP done\n" );
    cat( "PCA and TSNE ran\n" );

    cat( "Finding Neighbors\n" );
    s.in.inte <- FindNeighbors(s.in.inte, assay="integrated", reduction = "pca", dims = 1:nDims, force.recalc=T );
    cat( "Finding Clusters\n" );
    s.in.inte <- FindClusters(s.in.inte, assay="integrated", resolution = 1.0 )
    
    cat( "table of clusters:", table( s.in.inte@active.ident ), "\n" );;
    cat( "table of time: ", table( s.in.inte@meta.data$time ), "\n" );
    cat( "done with f.cluster.seuratv3!\n" );
    return( s.in.inte );
}

#####
f.entropy <- function( m.in ){
    m.in <- apply( m.in, 2, function( x ){ x / sum(x, na.rm=T ); } );
    v.entropy <- apply( m.in, 2, function(x){ y <- unname(x); y <- y[ y > 0 ]; sum( -y * log(y) ); } );
    return( v.entropy );
}


f.compute.time.in.clusters <- function( s.in ){
    df.tmp <- data.frame( cluster=s.in@active.ident, time=s.in@meta.data$time );
    df.tmp$cluster <- factor( df.tmp$cluster, levels=mixedsort( levels( df.tmp$cluster ) ) );
    df.tmp$cluster.perc <- as.vector( 100 / table( df.tmp$cluster )[ df.tmp$cluster ] );
    print( "sorted" );
    #head( df.tmp )
    p.frac <- ggplot() + geom_bar( data=df.tmp, aes( x=cluster, fill=time, weight=cluster.perc ) );
    p.frac <- p.frac + ylab( "percentage of cells" ) + ggtitle( "time points in clusters (%)" );
    p.frac <- p.frac + coord_flip();
        
    p.sizes <- ggplot() + geom_bar( data=df.tmp, aes( x=cluster, fill=cluster ) ) + coord_flip();
    p.sizes <- p.sizes + theme( panel.grid.major='element_line'(colour="grey40", size=0.3 ), panel.grid.minor='element_line'(colour="grey50", size=0.1 ) );
    
    l.out <- list();
    l.out[[ "time.in.clusters.plot" ]] <- p.frac;
    l.out[[ "cluster.sizes.plot" ]] <- p.sizes; 
    return( l.out );
}


f.seurat2meta <- function( s.dat3, v.clusters, v.genes.keep, MinPeriod, MaxPeriod, nTimePoints, Parallelize, nCoresMeta ){
    cat( "f.seurat2meta:\n" );
    dirMETAout <- "meta_out";
    dir.create( dirMETAout, showWarnings=FALSE );
    nSampGlob <- 1000;
    df.meta <- data.frame();
    df.meta.rb <- data.frame();
    v.times <- mixedsort( unique( s.dat3@meta.data$time ) );
    
    for ( cluster in v.clusters ){
    #for ( cluster in c( "6", "5" ) ){
        cat( "Working on cluster", cluster, "\n" );
        s.loc <- subset( s.dat3, idents = c(cluster) );
        
        df.data.loc <- data.frame( list( geneSymbol=rownames( s.loc$RNA@data ) ), stringsAsFactors=FALSE );
        
        for ( zt in v.times ){
            cat( "\t", zt, ":", sep="" );
            v.cn <- colnames( s.loc$RNA@data )[ s.loc@meta.data$time == zt ];
            cat( "length( v.cn )=", length( v.cn ), "\n" );
            v.rep <- paste( s.in@meta.data$date, s.in@meta.data$ar, s.in@meta.data$time, sep=":" )
            if ( length( v.cn ) == 0 ){
                #m.mean <- matrix( NA, nrow=nrow( s.loc$RNA@data ), ncol = 1 );
                #colnames( m.mean ) <- paste0( zt, ".Rep", "0" );
                m.mean <- matrix( NA, nrow=nrow( s.loc$RNA@data ), ncol = 0 );
            }else {
                nSamp <- ifelse( nSampGlob > length( v.cn ), length( v.cn ), nSampGlob );
                cat( "\t\tnSamp=", nSamp, "\n" );
                if ( length( v.cn ) < nSamp ){ print( "ERROR" ); stop(); }
                v.rs <- sample( 1:nSamp, nSamp, replace=F ); 
                v.labels <- rep( v.rs, ceiling( length( v.cn )/nSamp ) )[ 1:length( v.cn ) ];
                v.labels <- as.character( v.labels );
                names( v.labels ) <- v.cn;
                f.cn <- factor( v.labels );
                
                l.split <- split( v.cn, f.cn );
                m.mean <- sapply( names( l.split ),
                                 function( ns ){
                                     m.loc <- exp( s.loc$RNA@data[ , l.split[[ ns ]] ] ) - 1;
                                     if ( is.vector( m.loc ) ){ return( m.loc ); }
                                     else{ apply( m.loc, 1, mean ); }
                                 } );
                
                colnames( m.mean ) <- paste0( zt, ".Rep", names( l.split ) );
            }
            df.data.loc <- cbind( df.data.loc, m.mean );
        }
        df.data.loc <- df.data.loc[ df.data.loc$gene %in% v.genes.keep, ];
        cat( "nrow( df.data.loc ) = ", nrow( df.data.loc ), "\t", "ncol( df.data.loc ) = ", ncol( df.data.loc ),"\n" );
        
        rownames( df.data.loc ) <- df.data.loc$gene;
        fileLoc <- paste0( cluster, ".meta.input.csv" );
        cat( "fileLoc=", fileLoc, "\n" );
        write.csv( df.data.loc, file= paste0( dirMETAout, "/", fileLoc) , row.names=FALSE, quote=TRUE );
        #print( head( df.data.loc, n=2 ) );
        
        cat( "computing cycling\n" );
        fileOutLoc <- paste0( dirMETAout, "/", "meta2d_", fileLoc );
        cat( "fileOutLoc=", fileOutLoc, "\n" );
        unlink( fileOutLoc );
        v.timepoints <- as.numeric( gsub( "^zt(\\d+).Rep.*", "\\1", colnames( df.data.loc )[ 2:ncol(df.data.loc) ] ) );
        cat( "v.timepoints=", v.timepoints, "\n" );
        if ( length( unique( v.timepoints ) ) != nTimePoints ){
            cat( "WARNING: skipping cluster", cluster, "because cells are not present at every time point\n" );
            next;
        }
        MetaCycle::meta2d( infile=paste0( dirMETAout, "/", fileLoc), outdir=dirMETAout, filestyle="csv", timepoints=v.timepoints, minper=MinPeriod, maxper=MaxPeriod, cycMethod = c( "JTK", "LS" ), parallelize=Parallelize, nCores=nCoresMeta );
        cat( "MetaCycle has finished for cluster", cluster, "\n" );
        ##MetaCycle::meta2d( infile="cycMouseLiverProtein.txt", filestyle="txt", outdir="example", timepoints=rep(seq(0, 45, by=3), each=3), cycMethod=c("JTK","LS","ARS"), outIntegration="noIntegration" );
        df.out.loc <- read.csv( fileOutLoc );
        df.out.loc.rb <- df.out.loc;
        df.out.loc.rb$cluster <- cluster;
        colnames( df.out.loc )[ colnames( df.out.loc ) != "CycID" ] <- paste( colnames( df.out.loc )[ colnames( df.out.loc ) != "CycID" ], cluster, sep=":" );
        if ( ncol( df.meta ) == 0 ){ df.meta <- df.out.loc; df.meta.rb <- df.out.loc.rb;}
        else{
            df.meta <- merge( df.meta, df.out.loc, by=c( "CycID" ), all=TRUE );
            df.meta.rb <- rbind( df.meta.rb, df.out.loc.rb );
        }
        cat( "cluster", cluster, "meta results have been merged\n");
    }
    cat( "Writing output file\n" );
    fileMETAout <- paste0( dirMETAout, "/", "META_merged_results_", MinPeriod, "-", MaxPeriod, ".csv" );
    write.csv( df.meta, file=fileMETAout, row.names=FALSE, quote=FALSE );
    rownames( df.meta ) <- df.meta$CycID;
    fileMETAout.rb <- paste0( dirMETAout, "/", "META_appended_results_", MinPeriod, "-", MaxPeriod, ".csv" );
    write.csv( df.meta.rb, file=fileMETAout.rb, row.names=FALSE, quote=FALSE );
    return( list( cb=df.meta, rb=df.meta.rb ) );
}

f.computeTimeMeansSem <- function( s.in, v.clusters, v.genes.keep ){
    if ( ! "time" %in% names( s.in@meta.data ) ){ stop( "ERROR: time not defined" ); }

    s.in <- subset( s.in, features = v.genes.keep );
    
    v.times <- mixedsort( unique( s.in@meta.data$time ) );
    cat( "v.times=", v.times, "\n" );
    v.clusters <- unique( c( "all", v.clusters ) );
    a.d <- array( NA, c( nrow( s.in$RNA@counts), length( v.clusters ), length( v.times ), 2 ), dimnames=(list( "gene"=rownames( s.in$RNA@counts), "cluster"=v.clusters, "time"=v.times, "type"=c( "mean", "sem") ) ) );  

    df.c.data <- data.frame( list( geneSymbol=rownames( s.in$RNA@data ) ), stringsAsFactors=FALSE );
    for ( cluster in v.clusters ){
        print( cluster );
        v.cn <- vector();
        if ( cluster == "all" ){
            s.loc <- s.in;
        }
        else{
            s.loc <- subset( s.in, idents = c( cluster ) );
        }
        cat( "cluster=", cluster, "ncol( s.loc$RNA@data )=", ncol( s.loc$RNA@data ), "\n" );
        

        for ( zt in v.times ){
            #v.cn <- WhichCells( s.loc, expression='time == zt' );
            v.cn <- s.loc@meta.data$cell.names[ s.loc@meta.data$time == zt ];
            cat( "zt=", zt, "length( v.cn )=", length( v.cn ), "\n" );
            
            if ( length( v.cn ) > 0 ){
                a.d[ , cluster, zt, "mean" ] <- apply( as.matrix( exp( s.loc$RNA@data[ , v.cn ] ) -1 ), 1, mean, na.rm=T );
                a.d[ , cluster, zt, "sem" ] <-  apply( as.matrix( exp( s.loc$RNA@data[ , v.cn ] ) -1 ), 1, function( x ){ x <- x[ ! is.na( x ) ]; sd( x ) / sqrt( length( x ) ); } );
            }
            else{
                a.d[ , cluster, zt, "mean" ] <- NA;
                a.d[ , cluster, zt, "sem" ] <- NA;
            }
            
        }
            
     }
    return( a.d );
}



# The default seurat object is s.both.inte

f.plot.GE.6.timepoints <- function(Gene){
  
  df <- GetAssayData(s.in,assay = "RNA",slot = "data");
  df <- as.data.frame(df[Gene,]);
  df$experiment <- s.in@meta.data[rownames(df),]$experiment;
  df$time <- s.in@meta.data[rownames(df),]$time;
  df$cluster <- s.in@active.ident;
  colnames(df) <- c("Gene","experiment","time","cluster");
  df$Gene <- exp(df$Gene)-1;
  
  df <- dplyr::group_by(df,experiment,time,cluster) %>% dplyr::summarise("avg" = mean(Gene),"SEM" = sd(Gene)/sqrt(dplyr::n()));
  df <- as.data.frame(df);
  
  # df$cluster <- factor(df$cluster,levels = mixedsort(as.character(unique(df$cluster))))
  df$cluster <- factor(df$cluster,levels = c("2:s_LNv","25:l_LNv","8:LN_ITP","5:LNd_Trissin","9:LNd","12:LNd"))
  
  pdf(paste0(Gene," expression in all clusters.pdf"),width = 10,height = 8);
  
  # scaleFUN <- function(x) sprintf("%.1f", x) scale_y_continuous(labels=scaleFUN)
  # scales="free_y"
  p <- ggplot2::ggplot(df,aes(x=time,y=avg,group = experiment,color=experiment)) + geom_point() + geom_line(size=1.5) + theme_cowplot() + facet_wrap(~cluster,scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45,size = 12,vjust = 0.5),legend.text = element_text(size = 8),legend.title = element_blank()) + 
    geom_errorbar(aes(ymin=avg-SEM, ymax=avg+SEM), width=.2) +  labs(title = paste0(Gene," expression (TP10K)"),y="Mean(TP10K)") + xlab("") +
    scale_color_manual(values=c("#8F8F8F","#00C5CD"))
  print( p );
  dev.off();
  
}

