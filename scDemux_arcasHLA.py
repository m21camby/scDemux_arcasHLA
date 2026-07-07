
import os
from argparse import ArgumentParser
import sys
import pandas as pd
import anndata as an
import ast
import glob
import re


from src.functions import arcasHLA_pipeline


if __name__ == '__main__':

    parser = ArgumentParser(add_help = False)    
    parser.add_argument('-h', '--help', action='store_true') #print help
    parser.add_argument('-i', '--input_adata', default=None)
    parser.add_argument('-o', '--output_CB', default=None)
    parser.add_argument('-c', '--cluster_name', default= None)
    parser.add_argument('-l', '--cluster_list', type=str, default= None)
    parser.add_argument('-r', '--remove_prefix', default= None)
    parser.add_argument('-s', '--sample_id', default= None)
    parser.add_argument('-f', '--fastq_input_path', default= None)
    args = parser.parse_args()

    if args.help:
        print("####################################################################################################")
        print("--input_adata '/well/ovg-lambe/users/jdj468/data/processed/P2A/P2A_processed.h5ad'")
        print("--output_CB /well/ovg-lambe/users/jdj468/scripts/arcasHLA/analysis_results/P2A/' ")
        print("####################################################################################################")
        sys.exit(0)

        
# assign arguments
    PATH_adata = args.input_adata
    PATH_CB = args.output_CB
    cluster_name = args.cluster_name
    cluster_list = ast.literal_eval(args.cluster_list)
    Remove_Prefix = args.remove_prefix
    Sample_ID = args.sample_id
    PATH_FASTQ = args.fastq_input_path

# extract CB for each cluster
arcasHLA_pipeline.extract_Barcode(PATH_adata = PATH_adata,
                                 PATH_CB = PATH_CB,
                                 cluster_name = cluster_name,
                                 cluster_list = cluster_list,
                                 Remove_Prefix = Remove_Prefix,
                                 Sample_ID = Sample_ID)


# extract fastq from CB
for i in cluster_list:
    print(f"start extracting cluster : {i}")
    arcasHLA_pipeline.extract_fastq(PATH_FASTQ = PATH_FASTQ,
                                   PATH_FASTQ_out = PATH_CB,
                                   Sample_ID = Sample_ID,
                                   cluster_each = i,
                                   PATH_EACH_CB = PATH_CB + Sample_ID + "_" + i + "_CB.txt")

# merged extract fastqs    
    arcasHLA_pipeline.merge_fastq(PATH_FASTQ_out = PATH_CB,
                                  cluster_each = i,
                                  Sample_ID = Sample_ID)



# Run arcasHLA from merged fastq
    arcasHLA_pipeline.run_arcasHLA(PATH_FASTQ_out = PATH_CB,
                     Sample_ID = Sample_ID,
                     cluster_each = i)









