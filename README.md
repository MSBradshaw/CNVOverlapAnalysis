# CNVOverlapAnalysis

The purpose of this repo is to document a series of experiements for determining criteria for saying if two CNV calls are the same. It accompanies the tool [SeeNV](https://github.com/MSBradshaw/SeeNV) as backgorund research. 
The calls processed in this experiment are derived from the 1000 genomes project data using whole exome sequencing data for all samples processed in the Broad Institute.
Calls were made using [GATK gCNV](https://github.com/MSBradshaw/GATK_CNV_caller), [CNVkit](https://github.com/MSBradshaw/CNVkitSnakeMake) and [Savvy CNV](https://github.com/MSBradshaw/SavvySnakemake).

All results here can be reproduced using the provided snakemake pipeline `do_all.smk`. 
The pipeline assumes there is a `Data/` with two files `all_calls.bed` and `all_SVs.bed`, which are too big for github. You can download them from google drive [here](https://drive.google.com/file/d/1s5jjXGS_dkLZu21HdpcwQ2ECEESPLIL5/view?usp=sharing).
