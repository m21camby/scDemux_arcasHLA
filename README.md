# scDemux_arcasHLA


The package provides a streamlined workflow for `HLA genotyping` of 5' single cell RNA-seq data specifically designed to validate sample demultiplexing results. By seamlessly integrating UMI-tools and arcasHLA, the tool offers an intuitive pipeline to cross-verify and confirm your multiplexed sample assignments.

```
python scDemux_arcasHLA.py -i /well/ovg-lambe/users/jdj468/data/processed/P4A/P4A_processed.h5ad -o /well/ovg-lambe/users/jdj468/scripts/scDemux_arcasHLA_results/P4A_HTO/ -c Classification -l '["HTO_1","HTO_2","HTO_3","HTO_4"]' -r P4A_ -s P4A -f /well/ovg-lambe/projects/IMMPROVE_CoP/FASTQ/GEX/GEX4A/

```
