# scDemux_arcasHLA


The package provides a streamlined workflow for `HLA genotyping` of 5' single cell RNA-seq data specifically designed to validate sample demultiplexing results. By seamlessly integrating UMI-tools and arcasHLA, the tool offers an intuitive pipeline to cross-verify and confirm your multiplexed sample assignments.

***
## How to install

```
git clone https://github.com/m21camby/scDemux_arcasHLA.git
cd scDemux_arcasHLA
conda env create -f environment.yml
```

This is C++ standard library version mismatch issue 
```
ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.29' not found
```

It occurs when C++ extension requires a newer version than the one installed in your operating system's global `lib64` directory. You can solve by below
```
export LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6
```

***

## How to use

```
-i = input of h5ad or csv or tsv file
-o = output path
-c = demultiplex column name
-l = list of hashtag or cluster 
-r = remove prefix from h5ad index
-s = sample name
-f = path to fastq file

e.g.
python scDemux_arcasHLA.py -i /data/processed/P4A/P4A_processed.h5ad -o /scripts/scDemux_arcasHLA_results/P4A_HTO/ -c Classification -l '["HTO_1","HTO_2","HTO_3","HTO_4"]' -r P4A_ -s P4A -f /projects/IMMPROVE_CoP/FASTQ/GEX/GEX4A/

python scDemux_arcasHLA.py -i /data/soc/P2A/clusters.tsv -o /scripts/scDemux_arcasHLA_results/P2A/ -c assignment -l '["0","1","2","3"]' -r P2A_ -s P2A -f /projects/IMMPROVE_CoP/FASTQ/GEX/GEX2A/

```
