# csv file_md5sum
Rscript calculate_file_md5sum.R \
  -f /share/projects/Analytics/IO/zouhua/Script/dataset/phenotype \
  -m NULL \
  -p csv \
  -o result
# #Rscript calculate_file_md5sum.R \
#   -f dataset/phenotype \
#   -m result/Gmetadata_from_Doctor_csv_md5sum.csv \
#   -p csv \
#   -o result

# metaphlan
python merge_metaphlan.py \
  -f dataset/metaphlan/metaphlan_test_filepath.tsv \
  -p merge_metaphlan \
  -t species \
  -o result
  
# Humann2_card
python merge_humann.py \
  -f dataset/humann/Humann2_card_test_filepath.tsv \
  -p merge_humann_card \
  -d humann2_card \
  -o result
  
# Humann2_keggMean
python merge_humann.py \
  -f dataset/humann/Humann2_keggMean_test_filepath.tsv \
  -p merge_humann_keggMean \
  -d humann2_keggMean \
  -o result

# Humann2_keggMedian
python merge_humann.py \
  -f dataset/humann/Humann2_keggMedian_test_filepath.tsv \
  -p merge_humann_keggMedian \
  -d humann2_keggMedian \
  -o result
  
# humann2 & metacyc
python merge_humann.py \
  -f dataset/humann/Humann2_pathabundance_test_filepath.tsv \
  -p merge_humann \
  -d humann2 \
  -o result

# Humann2_vfdb
python merge_humann.py \
  -f dataset/humann/Humann2_vfdb_test_filepath.tsv \
  -p merge_humann_vfdb \
  -d humann2_vfdb \
  -o result
  
# humann3 & metacyc
python merge_humann.py \
  -f dataset/humann/Humann3_pathabundance_test_filepath.tsv \
  -p merge_humann3 \
  -d humann3 \
  -o result


# phyloseq for BJ GI
Rscript generate_phyloseq_for_BJ.R \
  -m dataset/phenotype/bjch_sample_sequence_metadata_20220923.csv \
  -p dataset/metaphlan/metaphlan2_BJ_RoundB_merge.csv \
  -r RoundB \
  -t metaphlan \
  -n metaphlan2 \
  -o result
  

# Swimming plot for BJ GI cancer
Rscript BJ_GI_Swimming_plot.R \
  -p dataset/phenotype/bjch_sample_sequence_metadata_20220923.csv \
  -c colon_cancer \
  -n Swimming_Plot -w 10 -g 20 -o result
  

# Reads QC from pipeline
python obtain_metaphlan_ReadsQC.py \
  -f dataset/pipeline_result/BJ_RoundG \
  -p BJ_RoundG_ReadsQC \
  -o result
python obtain_metaphlan_ReadsQC.py \
  -f dataset/pipeline_result/BJ_RoundB \
  -p BJ_RoundB_ReadsQC \
  -o result
  

# phyloseq for BJ GI for three types' dataset
Rscript generate_phyloseq_for_BJ_v2.R \
  -m dataset/phenotype/bjch_sample_sequence_metadata_20221008.csv \
  -p dataset/profile/metaphlan2_BJ_RoundB-G_merge.csv \
  -r all \
  -t metaphlan \
  -n metaphlan2 \
  -o result


# Summarizing accounts of patients and samples in BJ GI cancer
Rscript Summarize_Patients_Samples.R \
  -a dataset/phyloseq/metaphlan2_BJ_all_ps.RDS \
  -p dataset/phyloseq/metaphlan2_BJ_all_ps_baseline_R6.RDS \
  -o result 


# Integrating phyloseq object per taxa levels into list object
Rscript integrate_taxa_phyloseq_for_BJ.R \
  -p ./dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_total.RDS \
  -c ./dataset/phyloseq/metaphlan2_class_BJ_all_ps_total.RDS \
  -o ./dataset/phyloseq/metaphlan2_order_BJ_all_ps_total.RDS \
  -f ./dataset/phyloseq/metaphlan2_family_BJ_all_ps_total.RDS \
  -g ./dataset/phyloseq/metaphlan2_genus_BJ_all_ps_total.RDS \
  -s ./dataset/phyloseq/metaphlan2_species_BJ_all_ps_total.RDS \
  -n total \
  -d ./result/


# Summarizing accounts of patients and samples in BJ GI cancer (version 2)
Rscript Summarize_Patients_Samples_v2.R \
  -a dataset/phyloseq/metaphlan2_BJ_all_metadata_total.csv \
  -b dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_filter.RDS \
  -c dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_baseR6.RDS \
  -d dataset/phyloseq/16s_BJ_all_metadata_total.csv \
  -e dataset/phyloseq/16s_BJ_all_metadata_filter.csv \
  -f dataset/phyloseq/16s_BJ_all_metadata_baseR6.csv \
  -n Summarize_Patients_Samples \
  -o result 


# checking whether the sra files were successfully downloaded
perl get_update_sra_NCBI.pl \
  -f SRR_Acc_List.txt \
  -d /share/projects/Analytics/IO/zouhua/pipeline/GSE164150/SRA/sra \
  -o SRR_Acc_List_redownload.txt
prefetch --option-file SRR_Acc_List_redownload.txt


# converting sra into fastq
find /share/projects/Analytics/IO/zouhua/pipeline/GSE164150/SRA/sra/ -name "*.sra" | \
  perl -e 'print"SampleID\tMode\tPath\n"; while(<>){chomp; $name=(split("\/", $_))[-1]; $name=~s/.sra//g; $type = "PAIRED"; print "$name\t$type\t$_\n";}' \
  > samples.path.tsv
perl sra2fq.pl \
  -f samples.path.tsv \
  -d fastq \
  -p sra_into_fastq
sh sra_into_fastq.sh


# scp file between local and remote
perl get_each_scp.pl \
  -f filepath.tsv \
  -s xbiome_102 \
  -d /share/projects/Analytics/IO/zouhua/database/ \
  -o scp_file


# rsync files between source and destination directory
python rsync.py \
  -s /share/projects/Analytics/IO/zouhua/pipeline/GSE164150/SRA/sra/ \
  -d /share/projects/Analytics/IO/zouhua/database \
  -p sra \
  -l sra_log.tsv


# batch run for scripts
python run_batch_shell.py \
  -m sra_into_fastq.sh \
  -o sra_into_fastq_batch.sh


# APAlyzer Expression
Rscript APAlyzer_Expression.R \
  -b bam_file.tsv \
  -r RData \
  -g mm9_REF.RData \
  -c chr19 \
  -e 3UTR \
  -o result
