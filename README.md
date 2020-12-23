# scRNA_seq_clock_neurons
Analysis of scRNA_seq of Drosophila clock neurons in LD and DD conditions around the clock.

The configurations are parameters used in this analysis are in clk856.config.R. The clustering of the cells collected in LD and DD conditions around the clock was done with integrate.LD.DD.R, confident clusters were subsequently filtered by subset_confident_clusters.R. seurat2cycling.LD.DD.R and seurat2time.series.48.R were used for rhythmic gene expression analysis. Differentia gene expression analysis was done with compute.DE.R. The codes for results visulization are in plot.R. Command lines to run the scripts on shell are deposited in RUN.sh

For more details, please see:

## A transcriptomic taxonomy of Drosophila circadian neurons around the clock

Dingbang Ma, Dariusz Przybylski, Katharine C. Abruzzi, Matthias Schlichting, Qunlong Li, Xi Long, Michael Rosbash

doi: https://doi.org/10.1101/2020.09.15.297051
