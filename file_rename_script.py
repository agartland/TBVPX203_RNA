# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:06:04 2024

@author: hsrivast
"""




#################### renaming the filenames from ptid to pubid#############################################################
import os
import pandas as pd

def extract_file_names(directory):
    file_data = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.startswith("203_"):  
                file_name = os.path.splitext(file)[0]  
                file_info = file_name.split('_S')[0]  
                file_data.append({'File_Name': file, 'Info_Before_S': file_info})
    return file_data

directory = "/fh/scratch/delete30/gilbert_p/data/"
file_data = extract_file_names(directory)

# Create a DataFrame with the file names and the extracted info
df = pd.DataFrame(file_data)
df['ptid_day'] = df['Info_Before_S'].str.replace('203_', '')
counts = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/salmon.merged.gene_counts_corrected.csv")


meta_data = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/salmon_metadata.csv", index_col=0)
treatment_file = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/trt_pubid_2023-NOV-28.csv")
counts = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/salmon.merged.gene_counts_corrected.csv")


# Modify meta_data 'ptid' column
meta_data['ptid'] = meta_data['ptid'].str.replace('203-', '')
meta_data['ptid'] = meta_data['ptid'].str.replace('-', '_')

# Merge meta_data and treatment_file based on the common column "ptid"
merged_data = pd.merge(meta_data, treatment_file, on='ptid', how='inner')

# Modify 'pubid' column
merged_data["pubid"] = merged_data['pubid'].str.replace('-', '_')

# Create 'ptid_day' column
merged_data['ptid_day'] = merged_data['ptid'] + '_' + merged_data['day'].astype(str)
df = pd.merge(df, merged_data[['ptid_day', 'pubid']], on='ptid_day', how='left')

df['New_File_Name'] = df.apply(lambda row: row['File_Name'].replace(row['File_Name'][:12], row['pubid']), axis=1)


#### renaming the files in the folder ####################################
for root, dirs, files in os.walk(directory):
    for file in files:
        if file in df['File_Name'].tolist():
            original_file_path = os.path.join(root, file)
            new_file_path = os.path.join(root, df.loc[df['File_Name'] == file, 'New_File_Name'].iloc[0])
            os.rename(original_file_path, new_file_path)
            
            
################################ creating bio attributes file ################################################



############################### creating a SRA META data file #################################################
sample_acess = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/SRA_Data_upload/attributes.tsv",sep="\t")
df3 = pd.DataFrame()
df3["Biosample accesion"] = sample_acess["accession"]
df3["Library Id"] = sample_acess["sample_name"]
df3["Title"] = "Homo sapien RNAseq from blood after ID93+GLA-SE vaccination: " + sample_acess["sample_name"] + " at day " + sample_acess["day"].astype(str)
df3["Library strategy"] = "RNA-Seq"
df3["Library source"] = "TRANSCRIPTOMIC"
df3["Library selection"] = "PolyA"
df3["Library layout"] = "paired"
df3["Platform"] = "ILLUMINA"
df3["Instrument model"] = "Illumina NovaSeq 6000"
df3["Design description"] = "RNA extracted from PAXgene tubes of whole blood. Globin mRNA removed using bead-based hybridization (Thermofisher GLOBINclear). Libraries created from polyadenylated RNA transcripts (Ilumina TruSeq, stranded). Sequencing on Illumina NovaSeq 6000 with 100 bp paired-end reads."
df3["Filetype"] = "fastq"
matching_filenames_list = []

# Iterate over rows in df3
for index, row in df3.iterrows():
    library_id = row['Library Id']
    matching_filenames = [filename for filename in df['New_File_Name'] if library_id in filename]
    
    # Store matching filenames in the list
    matching_filenames_list.append(matching_filenames)

# Convert the list to DataFrame
matching_filenames_df = pd.DataFrame(matching_filenames_list)

# Rename the columns to filename_1, filename_2, ...
matching_filenames_df.columns = [f'filename_{i+1}' for i in range(len(matching_filenames_df.columns))]

# Concatenate df3 and matching_filenames_df
df3 = pd.concat([df3, matching_filenames_df], axis=1)
df3= df3.drop[""]


############################################# replacing ptid with pubid ###############################
counts = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/salmon.merged.gene_counts_corrected.tsv")
merged_data['ptid_day'] = merged_data['ptid'] + '_' + merged_data['day'].astype(str)
counts.columns = counts.columns.str.replace('X', 'P')

name_mapping = dict(zip(merged_data["samplename"], merged_data["pubid_day"]))

# Replace the column names in counts using the mapping
counts.rename(columns=name_mapping, inplace=True)
counts.to_csv('/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/salmon.merged.gene_counts_corrected.csv', index=False)


elisa = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/elisa_adata.csv")
ics = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/ics_adata.csv")
wbdata = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/wbics_adata.csv")
module = pd.read_csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/module_adata.csv")


name_mapping = dict(zip(merged_data["ptid"], merged_data["pubid"]))


ics['ptid'] = ics['ptid'].replace(name_mapping)
elisa['ptid'] = elisa['ptid'].replace(name_mapping)
wbdata['ptid'] = wbdata['ptid'].replace(name_mapping)
module['ptid'] = module['ptid'].replace(name_mapping)

module.to_csv('/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/module_adata.csv', index=False)
wbdata.to_csv('/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/wbics_adata.csv', index=False)
ics.to_csv('/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/ics_adata.csv', index=False)
elisa.to_csv('/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/data/pub_repo/elisa_adata.csv', index=False)