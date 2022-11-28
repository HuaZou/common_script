import argparse as ap
import os
import re
import pandas as pd
import numpy as np
import glob


parser = ap.ArgumentParser(description='obtaining reads\' status from pipeline')
parser.add_argument('-f', '--folder', 
                    type=str,
                    help='folder of metaphlan/humann result',
                    default="none", 
                    required=True)                   
parser.add_argument('-p', '--prefix', 
                    type=str,
                    help='prefix of output file', 
                    required=True)                                        
parser.add_argument('-o', '--out', 
                    type=str, 
                    default="./",
                    help='output path', 
                    required=False)
args = parser.parse_args()


def file_to_path(foldername, foldertype):

    file2path = {}
    folder_list = glob.glob(foldername + "/*/" + foldertype + "/", recursive=True)
    for folders in folder_list:
        for file in os.listdir(folders):
            filename = folders + file
            if re.findall(r'kneaddata_read_counts.txt|fastp', filename):
                if os.path.exists(filename):
                    if foldertype == "fastp":
                        seqID = str(re.match(r'\S+fastp'
                                                r'\/fastp\_(\d+|\d+\_\d+|\d+\_\d+\_\d+)'
                                                r'\.log', filename)[1])
                    elif foldertype == "kneaddata":
                        seqID = str(re.match(r'\S+kneaddata'
                                                r'\/(\d+|\d+\_\d+|\d+\_\d+\_\d+)'
                                                r'\.kneaddata_read_counts.txt', filename)[1])
                file2path[seqID] = filename
        
    return(file2path)


def Read_KneaddataFile(filenames, samplenames, typenames):

    dat = pd.read_table(filenames, 
                        sep='\t', 
                        skip_blank_lines=True,
                        index_col=0)

    # drop columns with orphan
    dat = dat[dat.columns.drop(list(dat.filter(regex='orphan')))]

    if typenames == 1:
        # extract number row
        res = dat.iloc[0:1, 2:6]
        # Change the column names
        res.columns = ["trimmomatic_trim1", "trimmomatic_trim2", 
                       "kneaddata_rmhost1", "kneaddata_rmhost2"]
    elif typenames == 2:
        # extract number row
        res = dat.iloc[0:1, :6]
        res.columns = ["fastp_trim1", "fastp_trim2",
                       "trimmomatic_trim1", "trimmomatic_trim2", 
                       "kneaddata_rmhost1", "kneaddata_rmhost2"]                
    
    # Change the row indexes
    res.index = [samplenames]

    return(res)


def Read_FastpFile(filenames, samplenames):

    reads_count = []
    with open(filenames, 'r') as inf:
        for line in inf.readlines():
            line = line.strip()
            if re.findall(r'reads:\s+\d+', line):
                read_number = line.split(r': ')[1]
                reads_count.append(read_number)

    read_list = ['Raw_read1', 'Raw_read2', 'fastp_trim1', 'fastp_trim2']
    res = pd.DataFrame([reads_count], columns=read_list)
    
    # Change the row indexes
    res.index = [samplenames]

    return(res)


def Merge_table(filedict, foldertype, types = 1): 

    df = pd.DataFrame()
    for key in filedict:
        if foldertype == "fastp":
            temp_df = Read_FastpFile(filedict[key], key)
        elif foldertype == "kneaddata":
            temp_df = Read_KneaddataFile(filedict[key], key, typenames=types)
        if df.empty:
            df = temp_df
        else:
            frames = [df, temp_df]
            df = pd.concat(frames, axis=0)

    res = df.copy()

    # dropping ALL duplicate values
    res.drop_duplicates()
    
    return(res)


def Make_dir(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)


def main():
    dir = args.folder
    prefix = args.prefix
    out = args.out

    # fastp
    fastp_file2dic = file_to_path(dir, "fastp")
    if len(fastp_file2dic.keys()) != 0:
        fastp_table = Merge_table(fastp_file2dic, "fastp")

    # knead 
    knead_file2dic = file_to_path(dir, "kneaddata")
    if len(fastp_file2dic.keys()) != 0:
        knead_table = Merge_table(knead_file2dic, "kneaddata")
    else:
        knead_table = Merge_table(knead_file2dic, "kneaddata", types=2)
    
    # merge results
    if len(fastp_file2dic.keys()) != 0:
        df_res = pd.merge(fastp_table, knead_table, left_index=True, right_index=True)
    else:
        df_res = knead_table
    df_out = df_res.rename_axis('SeqID')

    Make_dir(out)
    outfile_name = out + "/" + prefix + ".tsv"
    df_out.to_csv(outfile_name, sep='\t', encoding='utf-8', index=True)

    print('Congratulations, Program Ended Without Problem')

main()
