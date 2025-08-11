This repository contains R and Python scripts for the analyses of the study “Human-Specific Spatiotemporal Transcriptomic and Regulatory Divergence in Postnatal Prefrontal Cortex Development”.

1. Transcriptional_development_matching.R : R script for matching the timeline of development between human and macaque.
2. Scrublet.py : Python script for removing doublets from single-cell data using Scrublet.
3. RNAdata_processing.R : R script for snRNA data processing following Seurat workflow.
4. ATACdata_processing.R : R script for scATAC data processing following Signac workflow.
5. Xeniumdata_processing.R : R script for spatial data processing following Seurat workflow.
6. RNAATACSpatialdataset_integration.R : R script for integrating multiple datasets using FindIntegrationAnchors and IntegrateData functions in Seurat. 
7. markersHeatmap.R : R script for plotting Heatmap of the expression of marker genes.
8. HumanMacaca_neuron_integration.R : R script for integrating neuron datasets from different species using FindIntegrationAnchors and IntegrateData functions in Seurat.
9. ligand_receptor_module.R : R script for analyzing and plotting ligand-receptor module features using CellChat and ggplot2. 
