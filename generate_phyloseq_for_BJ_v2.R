# generate the phyloseq object of WGS for BJ GI Cancer
## Fecal samples on sequencing: 1. WGS; 2. 16s V4
## Dataset types:
##  1. Total samples (WGS for phyloseq object; 16s for metadata table)
##  2. Filtering samples according to the following criterion
##    a. Treatment_class are not "excluded" or "other";
##    b. Response_6 are "R", "CR" or "NR";
##    c. Batch are in Round B to G;
##    d. Diagnosis are "gastric cancer", "esophageal carcinoma" and "colon cancer" (no execute)
##    e. Features' occurrence are more than 0.
##  3. Another dataset according to the unique Baseline samples and Response_6 and corresponding after samples
##    a. the unique baseline samples and after samples; 
##    b. the corresponding after samples based on the unique baseline sample; 
##    c. Response_6 are "R" , "CR", and NR".
## 

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(XMAS2))

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

option_list <- list(
  make_option(c("-m", "--metadata"), 
              type = "character",
              help = "metadata", 
              metavar = "character"),
  make_option(c("-p", "--profile"), 
              type = "character",
              help = "profile", 
              metavar = "character"),
  make_option(c("-r", "--round"), 
              type = "character",
              help = "specific round", 
              metavar = "character"), 
  make_option(c("-t", "--type"), 
              type = "character",
              help = "metaphlan or humann", 
              metavar = "character"),
  make_option(c("-l", "--taxalevel"), 
              type = "character",
              help = "Phylum/Class/Order/Family/Genus/Species",
              default = "Species",
              metavar = "character"), 
  make_option(c("-n", "--name"), 
              type = "character",
              help = "name of file", 
              metavar = "character"),  
  make_option(c("-o", "--out"), 
              type = "character",
              help = "output file path", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# input parameters
metadata <- read.csv(opt$metadata)
profile <- data.table::fread(opt$profile, header = TRUE)
round <- opt$round
type <- opt$type
taxalevel <- opt$taxalevel
name <- opt$name
out <- opt$out

# Script dir
# metadata <- read.csv("./dataset/phenotype/bjch_sample_sequence_metadata_20221008.csv")
# #profile <- data.table::fread("./dataset/profile/metaphlan2_BJ_RoundB-G_merge_phylum.csv", header = TRUE)
# profile <- data.table::fread("./dataset/profile/metaphlan2_BJ_RoundB-G_merge_species.csv", header = TRUE)
# round <- "all"
# type <- "metaphlan"
# taxalevel <- "Phylum"
# taxalevel <- "Species"
# name <- "metaphlan2"
# out <- "result"

# combined_round dir
# metadata <- read.csv("../Metadata/bjch_sample_sequence_metadata_20221008.csv")
# profile <- data.table::fread("./profile/metaphlan2_BJ_RoundB-G_merge_phylum.csv", header = TRUE)
# round <- "all"
# type <- "metaphlan"
# taxalevel <- "Phylum"
# name <- "metaphlan2_phylum"
# out <- "phyloseq_v2"

# n250 irAE 
# metadata <- read.csv("./Metadata/n250_irAE_Sequence_Metadata.csv")
# profile <- data.table::fread("./profile/metaphlan2_BJ_irAE_merge_phylum.csv", header = TRUE)
# round <- "all"
# type <- "metaphlan"
# taxalevel <- "Phylum"
# name <- "metaphlan2_phylum"
# out <- "phyloseq_v2"

############## 20220928: 2nd baseline definition###################
# https://confluence.xbiome.com/pages/viewpage.action?pageId=17039367
# 基线样本规则如下：0天表示第一次治疗
#  1. 优先选择0天作为基线样本
#  2. 若没有, 则考虑-30天以内的样本（离0最近的样本）作为基线样本
#  3. 若没有-30天以内样本, 则考虑+21天以内的样本（离0最近的样本）作为基线样本
#  4. 若以上样本都没有，该病人没有基线样本
#############################################
get_baseline <- function(dat) {
  
  # dat = phen_WGS
  
  uniq_PID <- unique(dat$SubjectID)
  dat$Stage_v2 <- NA
  phen_res <- data.frame()  
  for (i in 1:length(uniq_PID)) {
    # i = 27
    phen_pid <- dat[which(dat$SubjectID == uniq_PID[i]), , F]
    phen_baseline <- phen_pid[which(phen_pid$Stage == "baseline"), , F]
    phen_after <- phen_pid[which(phen_pid$Stage == "after"), , F]
    
    if (nrow(phen_baseline) != 0) {
      # 1. 0 days as baseline
      day0 <- which(phen_baseline$Day == 0)
      if (length(day0) != 0) {
        phen_base_final <- phen_baseline[day0, , F] %>%
          dplyr::slice(1)
      } else {
        # 2. close to day0 before day0
        before_day0 <- which(phen_baseline$Day < 0)
        phen_before_day0 <- phen_baseline[before_day0, , F]
        if (length(before_day0) != 0) {
          before_day0_close <- max(phen_before_day0$Day)
          phen_base_final <- phen_baseline[which(phen_baseline$Day == before_day0_close), , F] %>%
            dplyr::slice(1)
        } else {
          # 3. close to day0 after day0
          after_day0 <- which(phen_baseline$Day > 0)
          phen_after_day0 <- phen_baseline[after_day0, , F]
          after_day0_close <- min(phen_after_day0$Day)
          phen_base_final <- phen_baseline[which(phen_baseline$Day == after_day0_close), , F] %>%
            dplyr::slice(1)
        }
      }
      # before and after baseline samples
      phen_baseline$Stage_v2 <- ifelse(phen_baseline$SeqID == phen_base_final$SeqID,
                                       "baseline",
                                       ifelse(phen_baseline$Day == phen_base_final$Day,
                                              "same_basline",
                                              ifelse(phen_baseline$Day < phen_base_final$Day,
                                                     "pre_baseline",
                                                     "post_baseline")
                                       )
      )
    }
    if (nrow(phen_after) != 0) {
      phen_after$Stage_v2 <- "after"
    }
    
    phen_temp <- rbind(phen_baseline, phen_after)
    phen_res <- rbind(phen_res, phen_temp)
    
  }
  # phen_res %>%
  # dplyr::select(SubjectID, Day, Stage, Stage_v2) %>%
  # arrange(SubjectID, Day) -> a1
  columns <- c("SubjectID", "Sub_SubjectID", "SeqID", 
               "Gender", "Age", "Height", "Weight", "BMI", "BMI_category",
               "Diagnosis", "Diagnosis_new", 
               "Treatment_Origin", "Treatment", "Treatment_class", 
               "Therapy_Start_Date", "Sample.Collection.Date", 
               "Response_3", "Response_6",
               "MSI_Status", 
               "Analysis.Project", "Batch",
               "Day", "Stage", "Stage_v2",
               grep("irAE", colnames(dat), value = T))
  phen_final <- phen_res %>%
    dplyr::select(all_of(columns), everything())
  
  # reorder  Treatment_class (20221012)
  phen_final$Treatment_class[which(phen_final$Treatment_class == "ICI|chemo|target")] <- "ICI|target|chemo"
  phen_final$Treatment_class[which(phen_final$Treatment_class == "ICI|targe")] <- "ICI|target"
  phen_final$Treatment_class[which(phen_final$Treatment_class == "XELOX")] <- "chemo"
  
  rownames(phen_final) <- paste0("S_", phen_final$SeqID)
  
  # reorder columns and drop columns
  res <- phen_final %>%
    dplyr::select(all_of(columns), everything())
  
  return(res)
}

###### total accounts of the BJ GI patients with fecal samples of WGS and 16s ##########
get_total <- function(x = metadata,
                      y = profile,
                      round_name = c("all", "RoundCR",
                                     "RoundB", "RoundC", "RoundD", 
                                     "RoundE", "RoundF", "RoundG"),
                      profile_type = c("metaphlan", "humann"),
                      taxa_level = c("Phylum", "Class", "Order",
                                     "Family", "Genus", "Species")) {
  
  # x = metadata
  # y = profile
  # round_name = round
  # profile_type = type
  # taxa_level = taxalevel
  
  colnames(x)[which(colnames(x) == "SeqID")] <- "Seq.ID"
  # phenptype: Round and filter 16s V4 and Diagnosis
  phenotype <- x %>% 
    # dplyr::select(-Name) %>%
    # dplyr::select(columns, everything()) %>%
    dplyr::mutate(Batch = dplyr::case_when(
      Analysis.Project == "北肿_B轮" ~ "RoundB",
      Analysis.Project == "北肿_C轮" ~ "RoundC",
      Analysis.Project == "北肿_D轮" ~ "RoundD",
      Analysis.Project == "北肿_E轮" ~ "RoundE",
      Analysis.Project == "北肿_E轮_加测" ~ "RoundE",
      Analysis.Project == "北肿_F轮" ~ "RoundF",
      Analysis.Project == "北肿F轮自建库重测" ~ "RoundF",
      Analysis.Project == "北肿_G轮" ~ "RoundG",
      Analysis.Project == "北肿_CR项目" ~ "RoundCR",
      Analysis.Project == "古菌研究" ~ "Others",
      Analysis.Project == "北肿-口腔肠道" ~ "Others",
      Analysis.Project == "北肿_唾液" ~ "Others",
      Analysis.Project == "北肿联合治疗" ~ "Others",
      Analysis.Project == "北肿_患者" ~ "Others",
      Analysis.Project == "Research of 16S full length sequencing" ~ "Others",
      Analysis.Project == "" ~ "Others"
    )) %>%
    dplyr::rename(SeqID = Seq.ID) %>%
    # dplyr::filter(Content == "metagenomic") %>%
    dplyr::mutate(Diagnosis_new = 
                    ifelse(Diagnosis == "gastric cancer", "Gastric",
                           ifelse(Diagnosis == "esophageal carcinoma", "Esophageal",
                                  ifelse(Diagnosis == "colon cancer", "Colon", "Others")
                           )
                    )
    )
  
  # # select samples from B to G round
  # pheno <- phenotype %>%
  #   dplyr::filter(Batch %in% c("RoundB", "RoundC", "RoundD", 
  #                              "RoundE", "RoundF", "RoundG", 
  #                              "RoundCR"))
  
  # select samples by round 
  if (round_name == "all") {
    phen <- phenotype
  } else {
    phen <- phenotype %>%
      dplyr::filter(Batch %in% round_name)
  }
  
  # Sub_SubjectID with P
  phen$Sub_SubjectID <- paste0("P", phen$Sub_SubjectID)
  
  # BMI and category: https://www.code4example.com/r-lang/r-code-to-calculate-bmi/
  phen$BMI <- round(as.numeric(phen$Weight) / 
                    (as.numeric(phen$Height)/100)^2, 2)
  phen$BMI_category <- ifelse(phen$BMI < 18.5, "lean", 
                                   ifelse(phen$BMI < 25, "normal",
                                         ifelse(phen$BMI < 30, "overweight", 
                                                ifelse(is.na(phen$BMI), NA, "obese"))))
  
  # baseline and treatment: baseline; after
  phen$Day <- as.numeric(as.Date(phen$Sample.Collection.Date) - 
    as.Date(phen$Therapy_Start_Date))
  # 3 weeks after 1st therapy regards as baseline
  phen$Stage <- ifelse(phen$Day > 21, "after",
                      ifelse(phen$Day < -30 , "before_baseline", "baseline"))
  
  # WGS sample
  phen_WGS <- phen %>%
    dplyr::filter(Content == "metagenomic") %>%
    dplyr::filter(!Pipeline.Path %in% c("未测序/取消测序", "文件遗失"))
  
  # 16s sample
  phen_16 <- phen %>%
    dplyr::filter(Content %in% c("16S V4", "16S")) %>%
    dplyr::filter(!Pipeline.Path %in% c("未测序/取消测序", "文件遗失"))
  
  # WGS
  phen_WGS_baseline <- get_baseline(dat = phen_WGS)
  
  # 16s
  if (nrow(phen_16) == 0) {
    phen_16s_baseline <- NULL   
  } else {
    phen_16s_baseline <- get_baseline(dat = phen_16)    
  }
  
  # create WGS phyloseq 
  ## overlap of samples
  sid <- intersect(colnames(y), phen_WGS_baseline$SeqID)  
  
  ## profile
  prof <- y %>% 
    dplyr::select(all_of(c(colnames(y)[1], sid))) 
  colnames(prof)[1] <- "ID"
  
  if (profile_type == "metaphlan") {
    prof_cln <- import_metaphlan_taxa(data_metaphlan2 = prof,
                                      taxa_level = taxa_level)
    # otu_table
    otu_tab <- prof_cln$abu_tab
    colnames(otu_tab) <- paste0("S_", colnames(otu_tab))
    # tax_table
    tax_tab <- prof_cln$tax_tab
  } else {
    colnames(prof) <- paste0("S_", colnames(prof)) 
    otu_tab_temp <- prof %>%
      tibble::column_to_rownames("S_ID")
    otu_tab <- apply(otu_tab_temp, 2, function(x){as.numeric(x)}) %>%
      data.frame()
    rownames(otu_tab) <- rownames(otu_tab_temp)
  }
  
  ## phenotype subset
  phen_WGS_baseline_cln <- phen_WGS_baseline %>%
    dplyr::filter(SeqID %in% sid)
  
  ## phyloseq object
  otu_tab_res <- otu_tab %>%
    dplyr::select(rownames(phen_WGS_baseline_cln))
  
  if (!all(colnames(otu_tab_res)  == rownames(phen_WGS_baseline_cln))) {
    stop("Order of Samples between otu table and sample data is wrong")
  }
  
  # phyloseq object
  if (profile_type == "metaphlan") {
    ps_object <- get_metaphlan_phyloseq(
      sam_tab = phen_WGS_baseline_cln,
      otu_tab = otu_tab_res,
      tax_tab = tax_tab)
  } else {
    ps_object <- phyloseq::phyloseq(
      phyloseq::sample_data(data.frame(phen_WGS_baseline_cln)),
      phyloseq::otu_table(otu_tab_res, taxa_are_rows = TRUE))
  }
  
  ## trim occurrence less than 0
  phen_WGS_baseline_ps <- run_trim(ps_object, cutoff = 0, trim = "feature")
  
  res <- list(wgs=phen_WGS_baseline_ps,
              wgs_meta=phen_WGS_baseline,
              s16=phen_16s_baseline)
  
  return(res)
}

total_res <- get_total(x = metadata,
                       y = profile,
                       profile_type = type,
                       round_name = round,
                       taxa_level = taxalevel)

###### accounts of the BJ GI patients with fecal samples of WGS and 16s after filtering ##########
## Filtering samples according to the following criterion
##    a. Treatment_class are not "excluded" or "other";
##    b. Response_6 are "R", "CR" or "NR";
##    c. Batch are in Round B to G;
##    d. Diagnosis are "gastric cancer", "esophageal carcinoma" and "colon cancer" (no execute)
##    e. Features' occurrence are more than 0.
get_filtered_data <- function(dat,
                              data_type = c("wgs", "16s"),
                              round_name = c("all", "RoundCR",
                                             "RoundB", "RoundC", "RoundD", 
                                             "RoundE", "RoundF", "RoundG")) {
  
  # dat = total_res$wgs
  # data_type = "wgs"
  # profile_type = "metaphlan"
  
  if (data_type == '16s') {
    phenotype <- dat
  } else if (data_type == 'wgs') {
    phenotype <- dat@sam_data %>%
      data.frame()
  }
  
  phen <- phenotype
  # reorder  Treatment_class (20221012)
  phen$Treatment_class[which(phen$Treatment_class == "ICI|chemo|target")] <- "ICI|target|chemo"
  phen$Treatment_class[which(phen$Treatment_class == "ICI|targe")] <- "ICI|target"
  phen$Treatment_class[which(phen$Treatment_class == "XELOX")] <- "chemo"
  
  # remove Treatment_class is excluded and other (20221012)
  phen_cln_temp1 <- phen[which(phen$Treatment_class != "excluded"), , F]
  phen_cln_temp2 <- phen_cln_temp1[which(phen_cln_temp1$Treatment_class != "other"), , F]
  phen_cln <- phen_cln_temp2[which(phen_cln_temp2$Treatment_class != ""), , F]
  
  # remove patients without Response_6 of R or CR or NR (20221012)
  phen_cln2 <- phen_cln %>%
    # "R","失访","NR","R辅助","失访NR","术后辅助治疗","NR（2021-03后失访）","CR",""
    dplyr::filter(Response_6 %in% c("R", "CR", "NR"))
  
  # remove samples unmet Round B to round G
  phen_cln3 <- phen_cln2 %>%
    dplyr::filter(Batch %in% c("RoundB", "RoundC", "RoundD", 
                               "RoundE", "RoundF", "RoundG", 
                               "RoundCR")) 
  
  # remove Diagnosis
  # Diagnosis are "gastric cancer", "esophageal carcinoma" and "colon cancer"
  # phen_cln3 <- phen_cln3 %>%
  #   dplyr::filter(Diagnosis_new %in% c("Gastric", "Colon", "Esophageal"))   
  
  # recalculate the baseline sample
  phen_cln3_base <- get_baseline(dat = phen_cln3)
  
  # drop columns
  drop_columns <- c("Subject.ID",
                    "Sample.Name",
                    "Data.Source",
                    "Sample.Arrival.Date", 
                    "Sample.Remark",
                    "Config.Path")
  phen_res <- phen_cln3_base %>%
    dplyr::select(-all_of(drop_columns))  
  
  if (data_type == '16s') {
    res <- phen_res
  } else if (data_type == 'wgs') {
    ps_object <- dat
    phyloseq::sample_data(ps_object) <- phyloseq::sample_data(phen_res)
    res <- run_trim(ps_object, cutoff = 0, trim = "feature")
  }
  
  return(res)
}

# WGS
filter_wgs <- get_filtered_data(dat = total_res$wgs, 
                                data_type ="wgs")

# 16s
if (!is.null(total_res$s16)) {
  filter_16s <- get_filtered_data(dat = total_res$s16, 
                                  data_type ="16s")  
} else {
  filter_16s <- NULL
}


###### accounts of the BJ GI patients with unique baseline/after fecal samples and Response_6 (CR, R, NR) ##########
# Another dataset according to the unique Baseline and Response_6 and corresponding after samples
##    a. the unique baseline samples;
##    b. the corresponding after samples based on the unique baseline sample; 
##    c. Response_6 are "R" , "CR", and NR".
get_baseR6 <- function(dat, 
                       data_type = c("wgs", "16s")) {
  
  # dat = filter_wgs
  # data_type = "wgs"
  
  if (data_type == '16s') {
    phenotype <- dat
  } else if (data_type == 'wgs') {
    phenotype <- dat@sam_data %>%
      data.frame()
  }
  
  # Response 6 with R, R and NR
  phen <- phenotype %>%
    tibble::rownames_to_column("TempRowNames") %>%
    dplyr::filter(Response_6 %in% c("R", "CR", "NR")) %>%
    dplyr::mutate(Response_6_new = ifelse(Response_6 == "NR", "NR", "R")) %>%
    tibble::column_to_rownames("TempRowNames")
  
  # baseline
  phen_base <- get_baseline(dat = phen)
  phen_base_temp <- phen_base %>%
    tibble::rownames_to_column("TempRowNames") %>%
    dplyr::filter(Stage_v2 %in% c("baseline")) %>%
    tibble::column_to_rownames("TempRowNames")
  
  # the corresponding after samples based on the unique baseline sample
  phen_after <- phen_base %>%
    tibble::rownames_to_column("TempRowNames") %>%
    dplyr::filter(Stage_v2 %in% c("after")) %>%
    dplyr::filter(Sub_SubjectID %in% unique(phen_base_temp$Sub_SubjectID)) %>%
    tibble::column_to_rownames("TempRowNames")
  
  phen_res <- rbind(phen_base_temp, phen_after) %>%
    dplyr::arrange(Sub_SubjectID, Stage_v2)
  
  # reorder columns
  columns_order <- c("SubjectID", "Sub_SubjectID", "SeqID", 
               "Gender", "Age", "Height", "Weight", "BMI", "BMI_category",
               "Diagnosis", "Diagnosis_new", 
               "Treatment_Origin", "Treatment", "Treatment_class", 
               "Therapy_Start_Date", "Sample.Collection.Date", "Day", 
               "Response_3", 
               "Response_6", "Response_6_new",
               "MSI_Status", 
               "Analysis.Project", "Batch",
               "Stage", "Stage_v2",
               grep("irAE", colnames(dat), value = T))
  phen_final <- phen_res %>%
    dplyr::select(c(all_of(columns_order), everything()))
  
  if (data_type == '16s') {
    res <- phen_final
  } else if (data_type == 'wgs') {
    ps_object <- dat
    phyloseq::sample_data(ps_object) <- phyloseq::sample_data(phen_final)
    res <- run_trim(ps_object, cutoff = 0, trim = "feature")
  }
  
return(res)
}

# WGS
baseR6_wgs <- get_baseR6(dat = filter_wgs, 
                         data_type ="wgs")
# 16s
if (!is.null(filter_16s)) {
  baseR6_16s <- get_baseR6(dat = filter_16s, 
                           data_type ="16s") 
} else {
  baseR6_16s <- NULL
}


######### output ################################
if (!dir.exists(out)) {
  dir.create(out, recursive = T)
}
if (length(round) == 1) {
  round_name <- round
} else {
  round_name <- paste0("Round", 
                       paste(gsub("Round", "", round[order(round)]), 
                             collapse = "-"))
}

# filenames
## wgs
wgs_total_name <- paste0(out, "/", name, "_BJ_", round_name, "_ps_total.RDS")
wgs_total_metadata_name <- paste0(out, "/", "metaphlan2_BJ_", round_name, "_metadata_total.csv")
wgs_filter_name <- paste0(out, "/", name, "_BJ_", round_name, "_ps_filter.RDS")
wgs_baseR6_name <- paste0(out, "/", name, "_BJ_", round_name, "_ps_baseR6.RDS")

## 16s
if (!is.null(total_res$s16)) {
  s16_total_name <- paste0(out, "/", "16s", "_BJ_", round_name, "_metadata_total.csv")  
}
if (!is.null(filter_16s)) {
  s16_filter_name <- paste0(out, "/", "16s", "_BJ_", round_name, "_metadata_filter.csv") 
}
if (!is.null(baseR6_16s)) {
  s16_baseR6_name <- paste0(out, "/", "16s", "_BJ_", round_name, "_metadata_baseR6.csv")
}

# output
## wgs
saveRDS(total_res$wgs, wgs_total_name, compress = TRUE)
write.csv(total_res$wgs_meta, wgs_total_metadata_name, row.names = T)
saveRDS(filter_wgs, wgs_filter_name, compress = TRUE)
saveRDS(baseR6_wgs, wgs_baseR6_name, compress = TRUE)

## 16s
if (!is.null(total_res$s16)) {
  write.csv(total_res$s16, s16_total_name, row.names = T)
}
if (!is.null(filter_16s)) {
  write.csv(filter_16s, s16_filter_name, row.names = T) 
}
if (!is.null(baseR6_16s)) {
  write.csv(baseR6_16s, s16_baseR6_name, row.names = T)
}

message('Congratulations, Program Ended Without Problem')
