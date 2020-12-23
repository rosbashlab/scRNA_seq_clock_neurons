# define and set environment variable SC_CLK856. This is a directory having scr/ as subdirectory (where this project's scripts are located), e.g. export SC_CLK856="this_project_directory"

#conda activate scClock

## analyse LD and DD data together
time Rscript --vanilla $SC_CLK856/scr/integrate.LD.DD.R 2>&1 | tee Out-integrate.LD.DD.R.txt || { echo 'integrate.LD.DD.R failed' ; exit 1; }

## compute DE
time Rscript --vanilla $SC_CLK856/scr/compute.DE.R $SC_CLK856/analyses/integrated_LD_DD/s.both.inte.rds $SC_CLK856/analyses/integrated_LD_DD/LD.DD.markers.negbinom.csv 2>&1 | tee Out-compute.DE.R.LD.DD.txt || { echo 'compute.DE.R failed' ; exit 1; }


# subset confident clock neuron clusters

time Rscript --vanilla $SC_CLK856/scr/subset_confident_clusters.v3.R $SC_CLK856/analyses/integrated_LD_DD/ s.both.inte.rds $SC_CLK856/analyses/integrated_LD_DD/LD_DD.markers.negbinom.csv >& Out-subset.confident.clock.clusters.v3.txt || { echo 'subset_confident_clusters.v3.R failed' ; exit 1; }



# Cycler analysis

time Rscript --vanilla $SC_CLK856/scr/seurat2cycling.LD.DD.R $SC_CLK856/analyses/integrated_LD_DD/JTK_F24/ $SC_CLK856/analyses/integrated_LD_DD/s.both.inte.LD.rds CLK856_LD 2>&1 | tee Out-seurat2cycling.R.LD.txt || { echo 'seurat2cycling.LD.DD.R on LD failed' ; exit 1; }

time Rscript --vanilla $SC_CLK856/scr/seurat2cycling.LD.DD.R $SC_CLK856/analyses/integrated_LD_DD/JTK_F24/ $SC_CLK856/analyses/integrated_LD_DD/s.both.inte.DD.rds CLK856_DD 2>&1 | tee Out-seurat2cycling.R.DD.txt || { echo 'seurat2cycling.LD.DD.R on DD failed' ; exit 1; }

time Rscript --vanilla $SC_CLK856/scr/seurat2time.series.48.R $SC_CLK856/analyses/integrated_LD_DD/s.both.inte.rds CLK856_LD $SC_CLK856/analyses/integrated_LD_DD/JTK_F24/time.series.48.LD.rand.split.csv 2>&1 | tee Out-seurat2time.series.48.R.both.LD.txt || { echo 'seurat2time.series.48.R on LD failed' ; exit 1; }

time Rscript --vanilla $SC_CLK856/scr/seurat2time.series.48.R $SC_CLK856/analyses/integrated_LD_DD/s.both.inte.rds CLK856_DD $SC_CLK856/analyses/integrated_LD_DD/JTK_F24/time.series.48.DD.rand.split.csv 2>&1 | tee Out-seurat2time.series.48.R.both.DD.txt || { echo 'seurat2time.series.48.R on DD failed' ; exit 1; }

