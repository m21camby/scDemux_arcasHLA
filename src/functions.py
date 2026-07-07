import anndata as an
import os
import sys
import pandas as pd
import glob
import re

class arcasHLA_pipeline():

    def __init__(self, PATH_adata):
        self.PATH_adata = PATH_adata
        self.PATH_CB = PATH_CB
        self.cluster_name = cluster_name
        self.cluster_list = cluster_list
        self.Remove_Prefix = Remove_Prefix
        self.Sample_ID = Sample_ID

        self.PATH_FASTQ = PATH_FASTQ
        self.PATH_FASTQ_out = PATH_FASTQ_out
        self.cluster_each = cluster_each
        self.PATH_EACH_CB = PATH_EACH_CB

    # This function is to extract Cell barcode for each cluster from demultiplexing (e.g. souporcell or Hashtag)
    def extract_Barcode(PATH_adata = '/PATH/P1A_processed.h5ad',
                        PATH_CB = "/PATH/",
                        cluster_name = "assignment",
                        cluster_list = ['0','1','2','3'],
                        Remove_Prefix = "P1A_",
                        Sample_ID = "P1A"):
        
        # check if output path exist, if not, then create
        isExist = os.path.exists(PATH_CB)

        if not isExist:
        # Create a new directory because it does not exist
            os.makedirs(PATH_CB)
            print("------------- The new output directory is created! -------------")
        else:
            print("------------- Directory already exists! ------------------")


        # load adata 
        if "h5ad" in PATH_adata:
            adata = an.read_h5ad(PATH_adata)
            print("----------- adata load ---------------")

            # extract cell barcode within the cluster
            for i in cluster_list:
                print(f'{i} has been load') 
                soc_list = adata.obs[adata.obs[cluster_name].isin([i])].index.tolist()
            
                # remove suffix -1
                new_soc_list = [x.replace('-1', '') for x in soc_list]
            
                # remove prefix if exist
                if Remove_Prefix is None:
                    pass
                else:    
                    new_soc_list2 = [x.replace(Remove_Prefix, '') for x in new_soc_list]

        

                # save        
                with open(PATH_CB + Sample_ID + "_" + i + "_CB.txt", "a") as f:
                    for item in new_soc_list2:
                        f.write(str(item) + "\n")
        
        # load csv
        elif "csv" in PATH_adata:
            csv = pd.read_csv(PATH_adata, index_col = 0)
            print("----------- csv load ---------------")

            # extract cell barcode within the cluster
            for i in cluster_list:
                print(f'{i} has been load')
                soc_list = csv[csv[cluster_name].isin([i])].index.tolist()

                # remove suffix -1
                new_soc_list = [x.replace('-1', '') for x in soc_list]

                # remove prefix if exist
                if Remove_Prefix is None:
                    pass
                else:
                    new_soc_list2 = [x.replace(Remove_Prefix, '') for x in new_soc_list]



                # save        
                with open(PATH_CB + Sample_ID + "_" + i + "_CB.txt", "a") as f:
                    for item in new_soc_list2:
                        f.write(str(item) + "\n")
        

        # load tsv
        elif "tsv" in PATH_adata:
            tsv = pd.read_csv(PATH_adata, index_col = 0, sep = "\t")
            print("----------- csv load ---------------")

            # extract cell barcode within the cluster
            for i in cluster_list:
                print(f'{i} has been load')
                soc_list = tsv[tsv[cluster_name].isin([i])].index.tolist()

                # remove suffix -1
                new_soc_list = [x.replace('-1', '') for x in soc_list]

                # remove prefix if exist
                if Remove_Prefix is None:
                    pass
                else:
                    new_soc_list2 = [x.replace(Remove_Prefix, '') for x in new_soc_list]



                # save
                with open(PATH_CB + Sample_ID + "_" + i + "_CB.txt", "a") as f:
                    for item in new_soc_list2:
                        f.write(str(item) + "\n")




        # if neither h5ad or csv or tsv, then exit  
        else:
            print(" --------------- wrong type of file --------------- ")
            sys.exit(1)


    # This function is to extract fastq files from above cell barcode file
    def extract_fastq(PATH_FASTQ = "/well/ovg-lambe/projects/IMMPROVE_CoP/FASTQ/GEX/GEX10A/",
                      PATH_FASTQ_out = "/well/ovg-lambe/users/jdj468/scripts/arcasHLA_tmp/analysis_results/",
                      Sample_ID = "P1A",
                      cluster_each = "C0",
                      PATH_EACH_CB = "/well/ovg-lambe/users/jdj468/scripts/arcasHLA_tmp/analysis_results/P14A_6_CB.txt"):

        
        # get list of fastqs
        R1_list = glob.glob(PATH_FASTQ + "*R1*")
        R2_list = glob.glob(PATH_FASTQ + "*R2*")

        # sorted based on lane
        R1_list_sorted = sorted(R1_list, key=lambda x: int(re.search(r"L(\d+)_R1", x).group(1)))
        R2_list_sorted = sorted(R2_list, key=lambda x: int(re.search(r"L(\d+)_R2", x).group(1)))
   
        print(f"R1 list: {R1_list_sorted}")
        print(f"R2 list: {R2_list_sorted}")

        # Check output path exit
        isExist = os.path.exists(PATH_FASTQ_out)
        if not isExist:
        # Create a new directory because it does not exist
            os.makedirs(PATH_FASTQ_out)
            print("The new directory is created!")
        else:
            print("Directory already exists!")
        
        # extract fastqs by each cluster
        for i in range(len(R1_list_sorted)):
            print(R1_list_sorted[i])
            print(R2_list_sorted[i])
            cmd = "umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN " + "--stdin=" + R1_list_sorted[i] + " --stdout=" + PATH_FASTQ_out + Sample_ID + "_" + cluster_each  + "_Sub_" + R1_list_sorted[i][-20:] + " --read2-in=" + R2_list_sorted[i] + " --read2-out=" + PATH_FASTQ_out + Sample_ID + "_" + cluster_each + "_Sub_" + R2_list_sorted[i][-20:] + " --whitelist=" + PATH_EACH_CB
            os.system(cmd)

    # This function is to merge above extract fastq files into one merged file
    def merge_fastq(PATH_FASTQ_out = "/well/ovg-lambe/projects/IMMPROVE_CoP/FASTQ/GEX/GEX10A/",
                    cluster_each = "0",
                    Sample_ID = "P1A"):

        # R1 merge
        Sub_list = glob.glob(PATH_FASTQ_out + "*Sub*")
        Sub_list = [x for x in Sub_list if '_' + cluster_each + "_" in x]

        R1_list = [x for x in Sub_list if 'R1' in x]
        R1_list_sorted = sorted(R1_list, key=lambda x: int(re.search(r"L(\d+)_R1", x).group(1)))
        cmd = "cat " + " ".join(R1_list_sorted) + " > " + PATH_FASTQ_out + Sample_ID + "_" + cluster_each + "_merged_R1.fastq.gz"
        os.system(cmd)

        # R2 merge
        Sub_list = glob.glob(PATH_FASTQ_out + "*Sub*")
        Sub_list = [x for x in Sub_list if '_' + cluster_each + "_" in x]

        R2_list = [x for x in Sub_list if 'R2' in x]
        R2_list_sorted = sorted(R2_list, key=lambda x: int(re.search(r"L(\d+)_R2", x).group(1)))
        cmd = "cat " + " ".join(R2_list_sorted) + " > " + PATH_FASTQ_out + Sample_ID + "_" + cluster_each + "_merged_R2.fastq.gz"
        os.system(cmd)



    # Run arcasHLA from merged fastq file
    def run_arcasHLA(PATH_FASTQ_out = "/well/ovg-lambe/users/jdj468/scripts/arcasHLA_tmp/analysis_results/",
                     Sample_ID = "P1A",
                     cluster_each = "C0"):


        # Check output path exit
        isExist = os.path.exists(PATH_FASTQ_out + "C" + cluster_each )
        if not isExist:
        # Create a new directory because it does not exist
            os.makedirs(PATH_FASTQ_out + "C" + cluster_each)
            print("The new directory is created!")
        else:
            print("Directory already exists!")

        # Run arcasHLA
        cmd = "./arcasHLA genotype " + PATH_FASTQ_out + Sample_ID + "_" + cluster_each + "_merged_R1.fastq.gz " + PATH_FASTQ_out + Sample_ID + "_" + cluster_each + "_merged_R2.fastq.gz " + "-g A,B,C,DPB1,DQB1,DQA1,DRB1 -o " + PATH_FASTQ_out + "C" + cluster_each + " -t 8 -v" 
        os.system(cmd)




