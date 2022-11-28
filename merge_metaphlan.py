import argparse as ap
import os
import re
import pandas as pd
import numpy as np


parser = ap.ArgumentParser(description='combining results from metaphlan')
parser.add_argument('-f', '--file', type=str,
                    help='file contains metaphlan results', required=True)
parser.add_argument('-p', '--prefix', type=str,
                    help='prefix of output', required=True)
parser.add_argument('-t', '--taxa', type=str, default="species",
                    help='specific taxa level', required=False)
parser.add_argument('-k', '--kind', type=str, default="metaphlan2",
                    help='version of metaphlan', required=False) 
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
            FileName = str(re.match(r'\S+[humann2|humann3]\/'
                                    r'(\d+|\d+\_\d+|\d+\_\d+\_\d+)'
                                    r'\_metaphlan_\S+_list.tsv', 
                                    File)[1])
            file2path[FileName] = File
    return(file2path)


def Read_File(filenames, samplenames, taxalevels, filetype):
  
    if filetype == 'metaphlan2':
        # skip row with '#' and the 1st row
        dat = pd.read_table(filenames, 
                            sep='\t', 
                            skip_blank_lines=True, 
                            comment='#',
                            skiprows=1,
                            header=None)
        dat1 = dat.iloc[:, 0:2]
    else:
        # skip row with '#' 
        dat = pd.read_table(filenames, 
                            sep='\t', 
                            skip_blank_lines=True, 
                            comment='#',
                            #skiprows=1,
                            header=None)
        dat1 = dat.iloc[:, [0, 2]]  

    if taxalevels[0] == taxalevels[1]:
        df_temp = dat1[dat1.iloc[:, 0].str.contains(taxalevels[0])]
    else:
        # filtering the rows where taxa doesn't contains the specific taxa level
        df = dat1[dat1.iloc[:, 0].str.contains(taxalevels[0])]
        df_temp = df[df.iloc[:, 0].str.contains(taxalevels[1]) == False]
    # adding column name to the respective columns
    df_temp.columns = ['taxaid', samplenames]

    res = df_temp.copy()

    # sorting by taxaid
    res.sort_values('taxaid', inplace=True)

    # dropping ALL duplicate values
    res.drop_duplicates(
                        subset='taxaid', 
                        keep=False, 
                        inplace=True)

    return(res)


def Merge_Metaphlan(filedict, taxa, filetype):
    if taxa == "species":
        taxalevels = 's__|t__'
    elif taxa == "genus":
        taxalevels = 'g__|s__'
    elif taxa == "family":
        taxalevels = 'f__|g__'
    elif taxa == "order":
        taxalevels = 'o__|f__'
    elif taxa == "class":
        taxalevels = 'c__|o__'
    elif taxa == "phylum":
        taxalevels = 'p__|c__' 

    taxa_levels = taxalevels.split('|')
    if len(taxa_levels) != 2:
        print("please check your taxa")
        exit()
    
    df = pd.DataFrame()
    for key in filedict:
        temp_df = Read_File(filedict[key], key, taxa_levels, filetype)
        if df.empty:
            df = temp_df
        else:
            df = pd.merge(df,
                          temp_df,
                          on="taxaid",
                          how='outer')
    res = df.replace(np.nan, 0)
    return(res)


def Make_dir(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)


def main():
    file = args.file
    out = args.out
    taxa = args.taxa
    kind = args.kind
    prefix = args.prefix

    file2dic = Filename_to_filepath(file)
    df_res = Merge_Metaphlan(file2dic, taxa, kind)
    
    Make_dir(out)
    outfile_name = out + "/" + prefix + ".csv"
    df_res.to_csv(outfile_name, sep='\t', encoding='utf-8', index=False)

    print('Congratulations, Program Ended Without Problem')

main()
