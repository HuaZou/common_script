suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

option_list <- list(
  make_option(c("-a", "--wgs_total"), 
              type = "character",
              help = "phyloseq of wgs total", 
              metavar = "character"),
  make_option(c("-b", "--wgs_filter"), 
              type = "character",
              help = "phyloseq of wgs filter", 
              metavar = "character"),
  make_option(c("-c", "--wgs_baseR6"), 
              type = "character",
              help = "phyloseq of wgs baseR6", 
              metavar = "character"),
  make_option(c("-d", "--s16_total"), 
              type = "character",
              help = "phyloseq of 16s total", 
              metavar = "character",
              default = NULL),
  make_option(c("-e", "--s16_filter"), 
              type = "character",
              help = "phyloseq of 16s filter", 
              metavar = "character",
              default = NULL),
  make_option(c("-f", "--s16_baseR6"), 
              type = "character",
              help = "phyloseq of 16s baseR6", 
              metavar = "character",
              default = NULL),
  make_option(c("-w", "--width"), 
              type = "integer",
              help = "width of plot", 
              default = 18), 
  make_option(c("-g", "--height"), 
              type = "integer",
              help = "height of plot", 
              default = 8),  
  make_option(c("-n", "--name"), 
              type = "character",
              help = "name of file", 
              metavar = "character"),  
  make_option(c("-o", "--out"), 
              type = "character",
              help = "output file path", 
              metavar = "character",
              default = "./")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# input parameters
## WGS
#wgs_total <- readRDS(opt$wgs_total)
wgs_total <- read.csv(opt$wgs_total, row.names = 1)
wgs_filter <- readRDS(opt$wgs_filter)
wgs_baseR6 <- readRDS(opt$wgs_baseR6)
## 16s
if (all(!is.null(opt$s16_total), 
        !is.null(opt$s16_filter),
        !is.null(opt$s16_baseR6))) {
  s16_total <- read.csv(opt$s16_total, row.names = 1)
  s16_filter <- read.csv(opt$s16_filter, row.names = 1)
  s16_baseR6 <- read.csv(opt$s16_baseR6, row.names = 1)  
} else {
  s16_total <- NULL
  s16_filter <- NULL
  s16_baseR6 <- NULL  
}

## output
Width <- opt$width
Height <- opt$height
name <- opt$name
out <- opt$out

# # WGS
# # wgs_total <- readRDS("dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_total.RDS")
# wgs_total <- read.csv("dataset/phyloseq/metaphlan2_BJ_all_metadata_total.csv", row.names = 1)
# wgs_filter <- readRDS("dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_filter.RDS")
# wgs_baseR6 <- readRDS("dataset/phyloseq/metaphlan2_phylum_BJ_all_ps_baseR6.RDS")
# # 16s
# s16_total <- read.csv("dataset/phyloseq/16s_BJ_all_metadata_total.csv", row.names = 1)
# s16_filter <- read.csv("dataset/phyloseq/16s_BJ_all_metadata_filter.csv", row.names = 1)
# s16_baseR6 <- read.csv("dataset/phyloseq/16s_BJ_all_metadata_baseR6.csv", row.names = 1)

# Width <- 18
# Height <- 8
# name <- "Summarize_Patients_Samples"
# out <- "./"

all_round <- c("RoundB", "RoundC", "RoundD", 
               "RoundE", "RoundF", "RoundG",
               "RoundCR", "Others")

# get samples' and patients' number from phyloseq object
get_ps_account <- function(round_type,
                           ps,
                           type = NULL) {
  
  # ps = wgs_filter
  # round_type = "RoundB"
  
  metadata <- ps@sam_data %>%
    data.frame()
  
  if (is.null(type)) {
    meta <- metadata
  } else {
    meta <- metadata %>% 
      dplyr::filter(Stage_v2 %in% "baseline")
  }
  
  if (!round_type %in% unique(meta$Batch)) {
    res <- data.frame(Round = round_type,
                      Samples = 0,
                      Patients = 0)
    return(res)
  }
  
  phen <- meta %>%
    dplyr::filter(Batch == round_type)
  
  data_sample <- nrow(phen)
  data_patient <- length(unique(phen$SubjectID))
  
  res <- data.frame(Round = round_type,
                    Samples = data_sample,
                    Patients = data_patient)
  
  return(res)
}

# get samples' and patients' number from csv object
get_csv_account <- function(round_type,
                            dat,
                            type = NULL) {
  
  # dat = wgs_total
  # round_type = "Others"
  
  if (is.null(type)) {
    meta <- dat
  } else {
    meta <- dat %>% 
      dplyr::filter(Stage_v2 %in% "baseline")
  }  
  
  
  if (!round_type %in% unique(meta$Batch)) {
    res <- data.frame(Round = round_type,
                      Samples = 0,
                      Patients = 0)
    return(res)
  }
  
  phen <- meta %>%
    dplyr::filter(Batch == round_type)
  
  data_sample <- nrow(phen)
  data_patient <- length(unique(phen$SubjectID))
  
  res <- data.frame(Round = round_type,
                    Samples = data_sample,
                    Patients = data_patient)
  
  return(res)
}

# Calculate
## WGS
wgs_total_NO <- do.call(rbind, lapply(all_round, get_csv_account, wgs_total))
wgs_filter_NO <- do.call(rbind, lapply(all_round, get_ps_account, wgs_filter))
wgs_baseR6_NO <- do.call(rbind, lapply(all_round, get_ps_account, wgs_baseR6, "no"))
## 16S
if (all(!is.null(s16_total), 
        !is.null(s16_filter),
        !is.null(s16_baseR6))) {
  s16_total_NO <- do.call(rbind, lapply(all_round, get_csv_account, s16_total))
  s16_filter_NO <- do.call(rbind, lapply(all_round, get_csv_account, s16_filter))
  s16_baseR6_NO <- do.call(rbind, lapply(all_round, get_csv_account, s16_baseR6, "no"))
}


# plotting 
get_plot <- function(dat1,
                     dat2,
                     dat3,
                     dat_total1,
                     dat_total2,
                     datype = c("wgs", "16s"),
                     All_round = all_round) {

  # dat1 = wgs_total_NO
  # dat2 = wgs_filter_NO
  # dat3 = wgs_baseR6_NO
  # dat_total1 = wgs_total
  # dat_total2 = wgs_filter@sam_data %>% data.frame()
  # datype = "wgs"
  # All_round = all_round
  
  if (datype == "wgs") {
    columns <- c("MGS_Samples_Total", "MGS_Patients_Total",
                 "MGS_Samples_Filter", "MGS_Patients_Filter",
                 "MGS_Samples_baseR6", "MGS_Patients_baseR6")
    columns_level <- c(columns[-c(5:6)], "MGS_Samples_Patients_baseR6")
    columns_new <- c("MGS (Total)\nSamples", 
                     "MGS (Total)\nPatients", 
                     "MGS (Filtering)\nSamples",
                     "MGS (Filtering)\nPatients", 
                     "MGS (Baseline & Response_6)\nSamples & Patients")
  } else {
    columns <- c("S16_Samples_Total", "S16_Patients_Total",
                 "S16_Samples_Filter", "S16_Patients_Filter",
                 "S16_Samples_baseR6", "S16_Patients_baseR6")
    columns_level <- c(columns[-c(5:6)], "S16_Samples_Patients_baseR6")
    columns_new <- c("16s (Total)\nSamples", 
                     "16s (Total)\nPatients", 
                     "16s (Filtering)\nSamples",
                     "16s (Filtering)\nPatients", 
                     "16s (Baseline & Response_6)\nSamples & Patients")    
  }
  
  # rename
  colnames(dat1) <- c("Round", columns[1], columns[2])
  colnames(dat2) <- c("Round", columns[3], columns[4])
  colnames(dat3) <- c("Round", columns[5], columns[6])
  
  # merge
  merge_res_temp <- dat1 %>%
    dplyr::inner_join(dat2,
                      by = "Round") %>%
    dplyr::inner_join(dat3,
                      by = "Round")
  if (datype == "wgs") {
    merge_res <- merge_res_temp %>%
      dplyr::select(-MGS_Patients_baseR6) %>%
      dplyr::rename(MGS_Samples_Patients_baseR6=MGS_Samples_baseR6)
  } else {
    merge_res <- merge_res_temp %>%
      dplyr::select(-S16_Patients_baseR6) %>%
      dplyr::rename(S16_Samples_Patients_baseR6=S16_Samples_baseR6)   
  }
  
  # plotdata for barplot
  plotdata <- merge_res %>%
    tidyr::gather(key = "kind", value = "count", -Round) %>%
    dplyr::mutate(Round = factor(Round, levels = All_round),
                  kind = factor(kind, 
                                levels = columns_level,
                                labels = columns_new))
  
  # barplot
  pl_box <- ggplot(data = plotdata, aes(x = kind, y = count, fill = Round)) +
    geom_bar(stat = "identity", position="dodge", color = "black") +
    geom_text(aes(label = count), 
              position = position_dodge(0.9), 
              vjust = 1.6, size = 4)+
    labs(x = "", y = "Number") +
    scale_y_continuous(expand = c(0, 2)) +
    guides(fill = guide_legend(bycol = TRUE, override.aes = list(size = 2))) +
    theme_bw()+
    theme(axis.title.y = element_text(face = "bold", color = "black", size = 12),
          axis.title.x = element_text(face = "bold", color = "black", size = 12, vjust = -1.2),
          axis.text.y = element_text(face = "bold", color = "black", size = 10),
          axis.text.x = element_text(face = "bold", color = "black", size = 10,
                                     angle = 0, vjust = 0.5),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key.height = unit(0.6, "cm"),
          legend.text = element_text(face = "bold", color = "black", size = 10))
  
  # count table
  merge_res_sum_table <- data.frame(Round="Sum",
                          Samples_Total = sum(merge_res[, 2]),
                          Patients_Total = sum(merge_res[, 3]),
                          Samples_Filter = sum(merge_res[, 4]),
                          Patients_Filter = sum(merge_res[, 5]),
                          Samples_Patients_baseR6 = sum(merge_res[, 6]))
  
  colnames(merge_res_sum_table) <- c("Round", columns_level)
  
  # real count table
  merge_res_sum_table_real <- data.frame(Round="Sum (Real)",
                                    Samples_Total = nrow(dat_total1),
                                    Patients_Total = length(unique(dat_total1$Sub_SubjectID)),
                                    Samples_Filter = nrow(dat_total2),
                                    Patients_Filter = length(unique(dat_total2$Sub_SubjectID)),
                                    Samples_Patients_baseR6 = sum(merge_res[, 6]))
  colnames(merge_res_sum_table_real) <- c("Round", columns_level)
  
  
  merge_table <- rbind(merge_res, merge_res_sum_table, merge_res_sum_table_real) 
  colnames(merge_table) <- c("Round", columns_new)
  count_tab <- ggtexttable(merge_table, 
                           rows = NULL,
                           theme = ttheme("mOrange",
                                          base_size = 10,
                                          padding = unit(c(3, 3), "mm")))
  
  pl <- ggarrange(pl_box, count_tab,
                  ncol = 1, nrow = 2,
                  heights = c(1, 0.5))
  
  res <- list(pl=pl,
              tab=merge_table)
  
  return(res)
}

# WGS 
wgs_plot <- get_plot(
      dat1 = wgs_total_NO,
      dat2 = wgs_filter_NO,
      dat3 = wgs_baseR6_NO,
      dat_total1 = wgs_total,
      dat_total2 = wgs_filter@sam_data %>% data.frame(),      
      datype = "wgs")

# 16s
if (all(!is.null(s16_total), 
        !is.null(s16_filter),
        !is.null(s16_baseR6))) {
  s16_plot <- get_plot(
    dat1 = s16_total_NO,
    dat2 = s16_filter_NO,
    dat3 = s16_baseR6_NO,
    dat_total1 = s16_total,
    dat_total2 = s16_filter,    
    datype = "16s")
  # plot
  pl <- ggarrange(wgs_plot$pl, s16_plot$pl, ncol = 2)
  
  # table
  tab <- wgs_plot$tab %>%
    dplyr::inner_join(s16_plot$tab, by = "Round")  
} else {
  # plot
  pl <- wgs_plot$pl
  
  # table
  tab <- wgs_plot$tab 
}


png_name <- paste0(out, "/", name, ".png")

ggsave(png_name,
       pl, width = Width, height = Height, dpi = 600,
       units = "in")

pdf_name <- paste0(out, "/", name, ".pdf")
ggsave(pdf_name,
       pl, width = Width, height = Height, dpi = 600,
       units = "in")

csv_name <- paste0(out, "/", name, ".csv")
write.csv(tab, csv_name, row.names = F)

message('Congratulations, Program Ended Without Problem')
