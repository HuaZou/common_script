# Script for daily routine

The scripts are used to deal with data in my daily work

- [Script for daily routine](#script-for-daily-routine)
  - [Merging metaphlan](#merging-metaphlan)
  - [Merging humann](#merging-humann)
  - [Calculating md5sum of csv file](#calculating-md5sum-of-csv-file)
  - [Generating phyloseq object in BJ GI cancer](#generating-phyloseq-object-in-bj-gi-cancer)
  - [Summarizing accounts of patients and samples in BJ GI cancer](#summarizing-accounts-of-patients-and-samples-in-bj-gi-cancer)
  - [Plotting swimming plot in in BJ GI cancer](#plotting-swimming-plot-in-in-bj-gi-cancer)
  - [Obtaining Reads' status from metaphlan pipeline](#obtaining-reads-status-from-metaphlan-pipeline)
  - [phyloseq for BJ GI for three types' dataset](#phyloseq-for-bj-gi-for-three-types-dataset)
  - [Integrating phyloseq object per taxa levels into list object](#integrating-phyloseq-object-per-taxa-levels-into-list-object)
  - [Summarizing accounts of patients and samples in BJ GI cancer (version 2)](#summarizing-accounts-of-patients-and-samples-in-bj-gi-cancer-version-2)
  - [checking whether the sra files were successfully downloaded](#checking-whether-the-sra-files-were-successfully-downloaded)
  - [converting sra into fastq](#converting-sra-into-fastq)
  - [scp file between local and remote](#scp-file-between-local-and-remote)
  - [rsync files between source and destination directory](#rsync-files-between-source-and-destination-directory)
  - [batch run for scripts](#batch-run-for-scripts)

## Merging metaphlan
```bash
# python merge_metaphlan.py \
#   -f dataset/metaphlan/metaphlan_test_filepath.tsv \
#   -p merge_metaphlan \
#   -t s__|t__ \
#   -o result
  
python merge_metaphlan.py \
  -f dataset/metaphlan/metaphlan_test_filepath.tsv \
  -p merge_metaphlan \
  -t species \
  -o result  
```

## Merging humann
```bash
python merge_humann.py \
  -f dataset/humann/Humann2_pathabundance_test_filepath.tsv \
  -p merge_humann \
  -d humann2 \
  -o result
```

## Calculating md5sum of csv file
```bash
Rscript calculate_file_md5sum.R \
  -f dataset/phenotype \
  -p csv \
  -o result
```

## Generating phyloseq object in BJ GI cancer
```bash
Rscript generate_phyloseq_for_BJ.R \
  -m dataset/phenotype/bjch_sample_sequence_metadata_20220923.csv \
  -p dataset/metaphlan/metaphlan2_BJ_RoundB_merge.csv \
  -r RoundB \
  -t metaphlan \
  -n metaphlan2 \
  -o result
```

## Summarizing accounts of patients and samples in BJ GI cancer
```bash
Rscript Summarize_Patients_Samples.R \
  -a dataset/phyloseq/metaphlan2_BJ_all_ps.RDS \
  -p dataset/phyloseq/metaphlan2_BJ_all_ps_baseline_R6.RDS \
  -o result 
```

## Plotting swimming plot in in BJ GI cancer
```bash
Rscript BJ_GI_Swimming_plot.R \
  -p dataset/phenotype/bjch_sample_sequence_metadata_20220923.csv \
  -c colon_cancer \
  -n Swimming_Plot \
  -w 10 \
  -g 20 \
  -o result
```

## Obtaining Reads' status from metaphlan pipeline
```bash
# fastp & kneaddata
python obtain_metaphlan_ReadsQC.py \
  -f dataset/pipeline_result/BJ_RoundG \
  -p BJ_RoundG_ReadsQC \
  -o result

# only kneaddata
python obtain_metaphlan_ReadsQC.py \
  -f dataset/pipeline_result/BJ_RoundB \
  -p BJ_RoundB_ReadsQC \
  -o result
```

## phyloseq for BJ GI for three types' dataset
```bash
Rscript generate_phyloseq_for_BJ_v2.R \
  -m dataset/phenotype/bjch_sample_sequence_metadata_20221008.csv \
  -p dataset/profile/metaphlan2_BJ_RoundB-G_merge.csv \
  -r all \
  -t metaphlan \
  -n metaphlan2 \
  -o result
```

## Integrating phyloseq object per taxa levels into list object
```bash
Rscript integrate_taxa_phyloseq_for_BJ.R \
  -p ./dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_total.RDS \
  -c ./dataset/phyloseq/metaphlan2_class_BJ_all_ps_total.RDS \
  -o ./dataset/phyloseq/metaphlan2_order_BJ_all_ps_total.RDS \
  -f ./dataset/phyloseq/metaphlan2_family_BJ_all_ps_total.RDS \
  -g ./dataset/phyloseq/metaphlan2_genus_BJ_all_ps_total.RDS \
  -s ./dataset/phyloseq/metaphlan2_species_BJ_all_ps_total.RDS \
  -n total \
  -d ./result/
```

## Summarizing accounts of patients and samples in BJ GI cancer (version 2)
```bash
Rscript Summarize_Patients_Samples_v2.R \
  -a dataset/phyloseq/metaphlan2_BJ_all_metadata_total.csv \
  -b dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_filter.RDS \
  -c dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_baseR6.RDS \
  -d dataset/phyloseq/16s_BJ_all_metadata_total.csv \
  -e dataset/phyloseq/16s_BJ_all_metadata_filter.csv \
  -f dataset/phyloseq/16s_BJ_all_metadata_baseR6.csv \
  -n Summarize_Patients_Samples \
  -o result
```

## checking whether the sra files were successfully downloaded
```bash
perl get_update_sra_NCBI.pl \
  -f SRR_Acc_List.txt \
  -d /share/projects/Analytics/IO/zouhua/pipeline/GSE164150/SRA/sra \
  -o SRR_Acc_List_redownload.txt
prefetch --option-file SRR_Acc_List_redownload.txt
```

## converting sra into fastq
```bash
find /share/projects/Analytics/IO/zouhua/pipeline/GSE164150/SRA/sra/ -name "*.sra" | \
  perl -e 'print"SampleID\tMode\tPath\n"; while(<>){chomp; $name=(split("\/", $_))[-1]; $name=~s/.sra//g; $type = "PAIRED"; print "$name\t$type\t$_\n";}' \
  > samples.path.tsv
perl sra2fq.pl \
  -f samples.path.tsv \
  -d fastq \
  -p sra_into_fastq
sh sra_into_fastq.sh
```

## scp file between local and remote
```bash
perl get_each_scp.pl \
  -f filepath.tsv \
  -s xbiome_102 \
  -d /share/projects/Analytics/IO/zouhua/database/ \
  -o scp_file
```

## rsync files between source and destination directory
```bash
python rsync.py \
  -s /share/projects/Analytics/IO/zouhua/pipeline/GSE164150/SRA/sra/ \
  -d /share/projects/Analytics/IO/zouhua/database \
  -p sra \
  -l sra_log.tsv
```

## batch run for scripts
```bash
python run_batch_shell.py \
  -m sra_into_fastq.sh \
  -o sra_into_fastq_batch.sh
```