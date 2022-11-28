import argparse as ap
import os
import re
import pandas as pd
import numpy as np


parser = ap.ArgumentParser(description='combining results from humann')
parser.add_argument('-f', '--file', type=str,
                    help='file contains humann results', required=True)
parser.add_argument('-p', '--prefix', type=str,
                    help='prefix of output', required=True)
parser.add_argument('-d', '--datatype', type=str,
                    help='datatype of column', required=True)                                        
parser.add_argument('-o', '--out', type=str, default="./",
                    help='output path', required=False)
args = parser.parse_args()


def Filename_to_filepath(filelist):

    file2path = {}
    with open(filelist, 'r') as f:
        files = f.readlines()
    for File in files:
        File = File.strip()
        if os.path.exists(File):
            FileName = str(re.match(r'\S+[humann2_card|humann2_kegg|humann2_vfdb|humann2|humann3]'
                                    r'\/(\d+|\d+\_\d+|\d+\_\d+\_\d+)'
                                    r'\_[genefamilies.tsv|genefamilies_des.tsv|'
                                    r'pathway_abundance_mean.tsv|pathway_abundance_mean_des.tsv｜'
                                    r'pathway_abundance_median.tsv|pathway_abundance_median_des.tsv|'
                                    r'pathabundance.tsv]', File)[1])
            file2path[FileName] = File
    return(file2path)


def Read_File(filenames, samplenames, type):
    # skip row with '#' and the 1st row
    dat = pd.read_table(filenames, 
                        sep='\t', 
                        skip_blank_lines=True, 
                        comment='#',
                        #skiprows=1,
                        header=None)
    res = dat.iloc[:, 0:2]

    # adding column name to the respective columns
    res.columns = [type, samplenames]
    return(res)


def Merge_humann(filedict, type):    
    df = pd.DataFrame()
    for key in filedict:
        temp_df = Read_File(filedict[key], key, type)
        if df.empty:
            df = temp_df
        else:
            df = pd.merge(df,
                          temp_df,
                          on=type,
                          how='outer')
    df_na = df.replace(np.nan, 0)

    res = df_na.copy()

    # sorting by type
    res.sort_values(type, inplace=True)

    # dropping ALL duplicate values
    res.drop_duplicates(
                        subset=type, 
                        keep=False, 
                        inplace=True)

    return(res)


def Make_dir(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)


def main():
    file = args.file
    out = args.out
    type = args.datatype
    prefix = args.prefix

    file2dic = Filename_to_filepath(file)

    df_res = Merge_humann(file2dic, type)

    Make_dir(out)
    outfile_name = out + "/" + prefix + ".csv"
    df_res.to_csv(outfile_name, sep='\t', encoding='utf-8', index=False)

    print('Congratulations, Program Ended Without Problem')

main()
