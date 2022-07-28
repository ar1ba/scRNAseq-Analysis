# scRNAseq-Analysis

This repository is a notebook of various scRNA sequencing analysis techniques used to process raw data from 10x CellRanger results.

The workflow is as follows:
1. Quality Control + Preprocessing
2. Detect and remove potential doublets (softwares used: DoubletFinder, DoubletDetection, Scrublet)
3. Cluster remaining cells
4. Apply cell type annotations
5. Pseudotime Trajectory analysis
