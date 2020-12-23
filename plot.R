setwd("~/Desktop/clk856_figures/");


library( magrittr )
library( ComplexHeatmap )
library( circlize )
library( ggpubr )
library( corrplot )
library( tidyverse )
library( cowplot )
library( gtools )
library( Seurat  )

options(stringsAsFactors = FALSE)

getwd()

rm( list=ls())

ls()

s.raw <- readRDS("~/Desktop/clk856_figures/raw_rds/s.both.inte.rds")
s.17 <- readRDS("~/Desktop/clk856_figures/raw_rds/s.both.ld.dd.inte.50.subset.clusters.newnames.2.rds")


# -------------------------------------------------- Dimplot of 17 clusters --------------------------------------------------


pdf("tnse_plot_17_clusters_renamed.pdf",width = 8,height = 8)
p <- DimPlot(s.17,label.size = 6,pt.size = 0.5,label = TRUE) + NoLegend() + 
  scale_y_continuous(sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis())
print( p ); dev.off();


# ------------------------------------- plot detected number of genes and UMIs in empty wells and positive wells -------------------------------------

s.dd <- readRDS("~/Desktop/sc_clk856/data_in_merged/s.dat.in.CLK856_DD.new.ex.rds")
df.dd <- s.dd@meta.data %>% select(nCount_RNA,nFeature_RNA,empty)
colnames(df.dd) <- c("nCount_RNA","nFeature_RNA","Wells")
table( df.dd$Wells )
df.dd$Wells <- ifelse(df.dd$Wells == "1","Control","Positive")

s.ld <- readRDS("~/Desktop/sc_clk856/data_in_merged/s.dat.in.CLK856_LD.new.ex.rds")
df.ld <- s.ld@meta.data %>% select(nCount_RNA,nFeature_RNA,empty)
colnames(df.ld) <- c("nCount_RNA","nFeature_RNA","Wells")
table( df.ld$Wells )
df.ld$Wells <- ifelse(df.ld$Wells == "1","Control","Positive")

df.all <- rbind(df.ld,df.dd) %>% arrange(desc(Wells))


pdf("umi and ngene correlation.pdf")
p <- ggplot(df.all,aes(x=nCount_RNA,y=nFeature_RNA,color=Wells)) + geom_point(alpha=0.8) + scale_x_log10(labels = scales::comma) +
  labs(x="Number of UMIs",y="Number of genes") + theme_cowplot()
print( p );  dev.off()


# ------------------------------------- Pearson gene expression correlation -------------------------------------

df <- data.frame(genes = rownames(s.raw))
cat("Repeat info:",unique(s.raw$Repeats))

for(Repeat in unique(s.raw$Repeats)){
  v.cells <- s.raw@meta.data$cell.names[s.raw$Repeats == Repeat];
  cat("Repeats:",Repeat," ncells:", length(v.cells),"\n")
  df[,Repeat] <- apply(expm1(s.raw@assays$RNA@data[rownames(s.raw),v.cells]),1,mean)
}

rownames(df) <- df$genes; df$genes <- NULL;


pdf("pearson.correlation.of.LD.data.pdf");

p <-  ggplot(df,aes(x=LD_1,y=LD_2)) + geom_point()  + stat_cor(method = "pearson") + theme_cowplot() + 
  theme() + labs(x="LD Repeat-1",y="LD Repeat-2") + scale_y_log10(labels = scales::comma) + scale_x_log10(labels = scales::comma)
print( p ); dev.off();


pdf("pearson.correlation.of.DD.data.pdf");

p <-  ggplot(df,aes(x=DD_1,y=DD_2)) + geom_point()  + stat_cor(method = "pearson") + theme_cowplot() + 
  theme() + labs(x="DD Repeat-1",y="DD Repeat-2") + scale_y_log10(labels = scales::comma) + scale_x_log10(labels = scales::comma)

print( p ); dev.off();


df <- data.frame(genes = rownames(s.raw))

cat("Experiments info:",unique(s.raw$experiment))

for(experiment in unique(s.raw$experiment)){
  v.cells <- s.raw@meta.data$cell.names[s.raw$experiment == experiment];
  cat("experiment:",experiment," ncells:", length(v.cells),"\n")
  df[,experiment] <- apply(expm1(s.raw@assays$RNA@data[rownames(s.raw),v.cells]),1,mean)
}

rownames(df) <- df$genes; df$genes <- NULL;

pdf("pearson.correlation.of LD and DD.data.pdf");

p <-  ggplot(df,aes(x=CLK856_LD,y=CLK856_DD)) + geom_point()  + stat_cor(method = "pearson") + theme_cowplot() + 
  theme() + labs(x="CLK856_LD",y="CLK856_DD") + scale_y_log10(labels = scales::comma) + scale_x_log10(labels = scales::comma)

print( p ); dev.off();


# --------------------------------------------- tSNE plot of raw clustering result, group_by repeats  ---------------------------------------------


pdf("dimplot.repeats.pdf")
p <- DimPlot(s.raw,group.by = "Repeats",pt.size = 0.5) ;
print( p ); dev.off();


# ----------------------------------- plot cell distribution from LD and DD -----------------------------------

d.tmp <- as.data.frame(s.raw@meta.data) %>% select(experiment,time,cell.names,seurat_clusters);
d.tmp <- d.tmp %>% group_by(experiment,time,seurat_clusters) %>% summarise(ncells = dplyr::n())
d.tmp$experiment <- factor(d.tmp$experiment,levels = rev(unique(d.tmp$experiment)))


# scale_fill_brewer(palette="Set2") to change the color palette, but doesn't improve too much.

pdf("clusters_cell_distribution.pdf");
p <- ggplot(d.tmp,aes(seurat_clusters,ncells,fill=time)) + geom_bar(stat="identity",position = "fill") + coord_flip() + xlab("") + ylab("") +
  scale_y_continuous(labels = scales::percent) + theme_cowplot() + theme(legend.position = "none") + facet_wrap(~experiment) 

print( p ); dev.off();

# Detected genes and UMI in each cluster ----------------------------------



pdf("~/Desktop/clk856_figures/figures/Number_UMI_raw_cluster.pdf",width = 15)
p <- ggplot(s.raw@meta.data,aes(x = seurat_clusters, y= nCount_RNA, fill=seurat_clusters)) + geom_violin() + theme_cowplot() + scale_y_log10() +
  expand_limits(x = 0, y = 0) +
  theme( legend.position = "none") + geom_boxplot(width=0.2, fill="white") + labs(x="", y="Number of UMI")
print( p ); dev.off();

pdf("~/Desktop/clk856_figures/figures/Number_nGene_raw_cluster.pdf",width = 15)
p <- ggplot(s.raw@meta.data,aes(x = seurat_clusters, y= nFeature_RNA, fill=seurat_clusters)) + geom_violin() + theme_cowplot() + scale_y_log10() +
  theme( legend.position = "none") + geom_boxplot(width=0.2, fill="white") + labs(x="", y="Number of Genes") + expand_limits(x = 0, y = 0)
print( p ); dev.off();


# --------------------- Plot core clock gene expression in all 39 clusters ----------------------

f.plot.GE.6.timepoints <- function(Gene, s.in = s.raw){
  
  df <- GetAssayData(s.in,assay = "RNA",slot = "data");
  df <- as.data.frame(df[Gene,]);
  df$experiment <- s.in@meta.data[rownames(df),]$experiment;
  df$time <- s.in@meta.data[rownames(df),]$time;
  df$cluster <- s.in@active.ident;
  colnames(df) <- c("Gene","experiment","time","cluster");
  df$Gene <- exp(df$Gene)-1;
  
  df <- df %>% as.data.frame()
  
  df <- dplyr::group_by(df,experiment,time,cluster) %>% dplyr::summarise("avg" = mean(Gene),"SEM" = sd(Gene)/sqrt(dplyr::n()));
  
  df <- as.data.frame(df);
  df$cluster <- factor(df$cluster,levels = mixedsort(as.character(unique(df$cluster))))
  
  df$experiment <- gsub("CLK856_","",df$experiment);
  df$time <- gsub("zt","",df$time);
  df[is.na(df)] = 0;
  
  pdf(paste0(Gene," expression in all clusters_fixed_y.pdf"),width = 8 );
  
  # scales="free_y"
  p <- ggplot2::ggplot(df,aes(x=time,y=avg,group = experiment,color=experiment)) + geom_line() + theme_bw() +
    facet_wrap(~cluster,ncol = 8) +
    theme(axis.text.x = element_text(angle = 90,size = 9,vjust = 0.5,hjust = 1),legend.text = element_text(size = 6),
          legend.position="bottom", legend.box = "horizontal",legend.title = element_blank()) +
    geom_errorbar(aes(ymin=avg-SEM, ymax=avg+SEM), width=.2) +  labs(title = paste0(Gene," expression (TP10K)"),x="",y="Mean") +
    scale_color_manual(values=c( "#8F8F8F", "#00C5CD")) +
    theme(strip.text.x = element_text(size = 10),strip.background = element_blank())
  
  print( p );
  dev.off();
  
}


f.plot.GE.6.timepoints("tim")
f.plot.GE.6.timepoints("Clk")
f.plot.GE.6.timepoints("vri")


#  --------------------------------------- Custom tSNE plot for marker genes ---------------------------------------


l.markers.1 <- c( "Pdf","opa","ITP","sNPF","CCHa1","AstC-R2","VGlut","gl","DIP-beta");
l.markers.2 <- c("AstC","CNMa","Dh31","AstA","Rh7","Trissin");

v.select.genes <- c( l.markers.1,l.markers.2 );
selectDir <- "scatterplots.TP10K";
dir.create( selectDir );
s.in <-  s.17;

for ( gene in v.select.genes ){
  fileOutpdf <- paste0( selectDir, "/", paste0( gene, ".tsne.TP10K.pdf" ) );
  pdf( fileOutpdf );
  df <- as.data.frame( s.in@reductions$tsne@cell.embeddings );
  df$TP10K <- exp( s.in$RNA@data[ gene, ] ) -1;
  pf <- ggplot( df, aes(x=tSNE_1, y=tSNE_2));
  pf <- pf + geom_point( data=df, aes( x=tSNE_1, y=tSNE_2 ), colour="grey75", size=1.0, alpha=0.5, shape=16 )
  
  pf <- pf + geom_point( data=df[ df$TP10K > 0, ], aes( x=tSNE_1, y=tSNE_2, color=TP10K ), size=1.0, alpha=0.5, shape=16 ) + 
    scale_color_gradientn( colors=c("grey75", "red1", "red2", "red3", "red4", "black" ) );
  pf <- pf + labs( title=gene ) + guides(fill=guide_legend(title="TP10K")) + theme_cowplot();
  print( pf );
  dev.off();
}


# --------------------------------------- tSNE plot of NPF expression in LD ---------------------------------------


s.in <- subset(s.17, cells = s.17$cell.names[s.17$experiment == "CLK856_LD"])

df <- as.data.frame( s.in@reductions$tsne@cell.embeddings );
df$TP10K <- exp( s.in$RNA@data[ "NPF", ] ) -1;
pf <- ggplot( df, aes(x=tSNE_1, y=tSNE_2));
pf <- pf + geom_point( data=df, aes( x=tSNE_1, y=tSNE_2 ), colour="grey75", size=1.5, alpha=0.5, shape=16 )
pf <- pf + geom_point( data=df[ df$TP10K > 0, ], aes( x=tSNE_1, y=tSNE_2, color=TP10K ), size=1.5, alpha=0.5, shape=16 ) + scale_color_gradientn( colors=c("grey75", "red1", "red2", "red3", "red4", "black" ) );
pf <- pf + labs( title="NPF" ) + guides(fill=guide_legend(title="TP10K")) + theme_cowplot();

pdf("tsne.NPF.LD.pdf")
print( pf ); dev.off();



# -------------------------------- Plot gene expression in specific cluster, curves --------------------------------


dir <- "cycler_clusters_curves"
dir.create("cycler_clusters_curves")

f.plot.GE.6.timepoints.cluster <- function(Gene,cluster,s.in){
  
  if(!dir.exists(dir)){
    cat("There is no folder to store results\n")
    cat("Creat a dirctory \n")
    dir.create(dir,showWarnings = FALSE)
  }
  
  s.loc <- subset(s.in,cells = s.in@meta.data[s.in$Idents == cluster,]$cell.names);
  df <- GetAssayData(s.loc,assay = "RNA",slot = "data");
  
  cat("cluster:",cluster,"ncells:",ncol(s.loc),"\n")
  df <- as.data.frame(df[Gene,]);
  df$experiment <- s.loc@meta.data[rownames(df),]$experiment;
  df$time <- s.loc@meta.data[rownames(df),]$time;
  df$cluster <- s.loc@active.ident;
  colnames(df) <- c("Gene","experiment","time","cluster");
  df$Gene <- exp(df$Gene)-1;
  
  df <- dplyr::group_by(df,experiment,time,cluster) %>% dplyr::summarise("avg" = mean(Gene),"SEM" = sd(Gene)/sqrt(dplyr::n()));
  df <- as.data.frame(df);
  df$cluster <- factor(df$cluster,levels = mixedsort(as.character(unique(df$cluster))))
  
  df$experiment <- gsub("CLK856_","",df$experiment);
  df$time <- gsub("zt","",df$time);
  
  pdf(paste0(dir,"/",Gene,"_",cluster, "_expression.pdf"));
  
  # scales="free_y"
  p <- ggplot2::ggplot(df,aes(x=time,y=avg,group = experiment,color=experiment)) + geom_point() + geom_line(size=1.5) + theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45,size = 12,vjust = -0.01),legend.text = element_text(size = 8),legend.title = element_blank()) +
    geom_errorbar(aes(ymin=avg-SEM, ymax=avg+SEM), width=.2) +  labs(title = paste0(Gene,"_",cluster," expression (TP10K)"),x="",y="Mean") +
    scale_color_manual(values=c("#8F8F8F","#00C5CD")) + expand_limits(x = 0, y = 0)
  print( p );
  dev.off();
  
}

f.plot.GE.6.timepoints.cluster("tim","19:DN2",s.in = s.17)
f.plot.GE.6.timepoints.cluster("DIP-beta","2:s_LNv",s.in = s.17)
f.plot.GE.6.timepoints.cluster("pdm3","1:DN1p_CNMa",s.in = s.17)
f.plot.GE.6.timepoints.cluster("CrebA","2:s_LNv",s.in = s.17)
f.plot.GE.6.timepoints.cluster("dpr8","15:DN1p_CNMa",s.in = s.17)
f.plot.GE.6.timepoints.cluster(Gene = "beat-Ic",cluster = "2:s_LNv",s.in = s.17)
f.plot.GE.6.timepoints.cluster(Gene = "beat-Va",cluster = "1:DN1p_CNMa",s.in = s.17)


# -------- Compute differentially expressed genes and gene expression corre --------


df.marker <- FindAllMarkers( s.17, test.use="negbinom", only.pos=FALSE, assay="RNA", latent.vars=c("time.date"), verbose=TRUE  );

write.csv( df.all, file = "marker.genes.all.17.cluster_computed_0702.csv" );

s.17 <- ScaleData(s.17,features = rownames(s.17@assays$RNA@data),assay = "RNA")


df.all <- df.marker %>% filter(p_val_adj < 0.05,avg_logFC > log(1.5));


df <- sapply(split(s.17$cell.names,s.17$Idents),function(cells){
    apply(s.17@assays$RNA@scale.data[unique(df.all$gene),cells],1,mean)
  })
  
  
df.cor <- cor(df);diag(df.cor) <- NA;

pdf("Marker.genes.correlation.17.clusters.pdf");
p <- ComplexHeatmap::Heatmap(df.cor,name = "DE.cor",row_dend_width = unit(2, "cm"),column_dend_height = unit(2,"cm"),cluster_columns = TRUE);
print( p ); dev.off(); 


# ----------------------------------- plot secretory signals expression heatmap -----------------------------------

v.features.npp <- scan("~/Documents/scRNA_seq/resource/neuropeptides_55.txt",what = character())
v.singal <- c(c("ChAT","VGlut","CG5549"),v.features.npp)
v.singal <- v.singal[v.singal %in% rownames(s.17)]

df.marker <- df.marker %>% filter(p_val_adj <0.05, avg_logFC > log(1.5));
v.genes.plot <- intersect(unique(df.marker$gene),v.singal)
  
df <- sapply(split(s.17$cell.names,s.17$Idents),function(cells){
    apply(s.17@assays$RNA@scale.data[v.genes.plot,cells],1,mean)
  })
    

pdf("secretory_signals_heatmap.pdf");
p <- ComplexHeatmap::Heatmap(df,cluster_columns = FALSE,cluster_rows = TRUE,name = "Scaled",row_dend_width = unit(2, "cm"),
                               column_title = "Signal");
print( p ); dev.off();
  

# ---------------------------------- Heatmap of Marker genes in 17 clusters ----------------------------------

df.markers <- read.csv("~/Desktop/clk856_figures/figures/marker.genes.all.17.cluster.renamed.csv",stringsAsFactors = FALSE)

all( unique(df.markers$cluster) == unique( Idents( s.17 ) ))

df.m.o <- df.markers %>% filter(p_val_adj < 0.05, avg_logFC > log(1.5)) 

v.cluster.level <- c("2:s_LNv","25:l_LNv","8:LN_ITP","5:LNd_Trissin","9:LNd_NPF","12:LNd","1:DN1p_CNMa","15:DN1p_CNMa",
                     "3:DN1a","4:DN1p","20:DN3","14:DN3","29:LPN","7:DN1p","6:DN1p","19:DN2","18:DN1p")

df.m.o$cluster <- factor(df.m.o$cluster,levels = mixedsort(v.cluster.level))


tmp <- sapply( mixedsort(v.cluster.level), function(x){ v <- df.m.o$gene[ df.m.o$cluster == x ]; v <- v[ 1:min(3,length(v) )]; return(v); } )

str( tmp );

v.features <- unique( unlist( as.vector( tmp ) ) ) %>% as.character();

df <- sapply(split(s.17$cell.names,s.17$Idents),function(cells){
  apply(s.17@assays$RNA@scale.data[v.features,cells],1,mean)
})

df <- df[,mixedsort(v.cluster.level)]

pdf("Top_3_Markers_in_17_cluster_heatmaps.pdf")
p <- Heatmap(df,cluster_rows = FALSE,cluster_columns = FALSE,name = "Expression",row_names_gp = gpar(fontsize = 8),border = F,
             col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
print( p ); dev.off();


# ------------------------------------------TFs, heatmap --------------------------------------------------------------------------


v.dna.binding.tfs <- read.csv("~/Documents/scRNA_seq/resource/DNA_binding_TFs.csv") %>% pull(name) %>% unique()
v.tfs.kca <- read.csv("~/Documents/scRNA_seq/resource/fllybase_transcription_factors_KCA.csv") %>% pull( SYMBOL ) %>% unique()

v.intersect.tfs <- intersect(v.dna.binding.tfs,v.tfs.kca)

df.m.o <- df.markers 

df.m.o$cluster <- factor(df.m.o$cluster,levels = c("2:s_LNv","25:l_LNv","8:LN_ITP","5:LNd_Trissin","9:LNd_NPF","12:LNd","1:DN1p_CNMa","15:DN1p_CNMa",
                                                   "3:DN1a","4:DN1p","20:DN3","14:DN3","29:LPN","7:DN1p","6:DN1p","19:DN2","18:DN1p"))


df.tfs <- df.m.o[ df.m.o$p_val_adj < 0.01 & df.m.o$avg_logFC > log(1.25) & df.m.o$gene %in% v.intersect.tfs , ];
head(df.tfs)



tmp <- sapply( c("2:s_LNv","25:l_LNv","8:LN_ITP","5:LNd_Trissin","9:LNd_NPF","12:LNd","1:DN1p_CNMa","15:DN1p_CNMa","3:DN1a","4:DN1p","20:DN3","14:DN3",
                 "29:LPN","7:DN1p","6:DN1p","19:DN2","18:DN1p"), function(x){ v <- df.tfs$gene[ df.tfs$cluster == x ]; v <- v[ 1:min(5,length(v) )]; return(v); } )
str( tmp );

v.features <- unique( unlist( as.vector( tmp ) ) ) %>% as.character();

df.markers.loc.grouped <- df.markers %>% filter(gene %in% v.intersect.tfs, p_val_adj < 0.01, avg_logFC > log(1.25))  %>% group_by(gene) %>% 
  summarise(ncluster = n()) %>% arrange(desc(ncluster))

v.multi.tfs <- df.markers.loc.grouped %>% filter(ncluster >= 3) %>% pull(gene) %>% unique()

v.tfs <- c(v.multi.tfs,v.features) %>% unique()

df <- sapply(split(s.in$cell.names,s.in$Idents),function(cells){
  apply(s.in@assays$RNA@scale.data[v.tfs,cells],1,mean)
})

pdf("~/Desktop/clk856_figures/figures/top5.tfs.heatmaps.pdf")
p <- Heatmap(df,cluster_rows = FALSE,cluster_columns = FALSE,name = "Expression",row_names_gp = gpar(fontsize = 8),border = F,
             col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
print( p ); dev.off();


# pie chart for cyclers ---------------------------------------------------

df.ld <- read_excel("~/Documents/sc_mannuscript/table_for_submission/Table_S2_Identified_rhythmic_genes_in_each_cluster.xlsx",
                    sheet = 1)

df.dd <- read_excel("~/Documents/sc_mannuscript/table_for_submission/Table_S2_Identified_rhythmic_genes_in_each_cluster.xlsx",
                    sheet = 2)


df.cyclers <- data.frame("ID"=c("Non_cycling","LD cycling","LD & DD cycling","DD cycling"),
                         "number" = c(4197,
                                      setdiff( unique(df.ld$Gene),unique(df.dd$Gene) ) %>% n_distinct(),
                                      intersect(unique(df.ld$Gene),unique(df.dd$Gene)) %>% n_distinct(),
                                      setdiff( unique(df.dd$Gene), unique(df.ld$Gene)) %>% n_distinct()))


# Create Data
Prop <- df.cyclers$number
myPalette <- brewer.pal(4, "Paired")

pdf( "~/Desktop/clk856_figures/figures/pie_chart_for_cyclers.pdf")
p <- pie(Prop, 
    labels = c("Non_cycling Transcripts","LD cyclers","LD & DD cyclers","DD cyclers"),
    clockwise = TRUE,
    border="white", col=myPalette )

print( p ); dev.off();



                                                                       
