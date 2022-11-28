suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(ggpubr))  
suppressPackageStartupMessages(library(ggtext))


# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

option_list <- list(
  make_option(c("-a", "--phyloseq"), 
              type = "character",
              help = "phyloseq from all rounds", 
              metavar = "character",
              default = NULL),
  make_option(c("-p", "--phenotype"), 
              type = "character",
              help = "phenotype with all patients from data platform", 
              metavar = "character",
              default = NULL), 
  make_option(c("-c", "--cancer"), 
              type = "character",
              help = "type of cancer", 
              metavar = "character",
              default = NULL),   
  make_option(c("-n", "--name"), 
              type = "character",
              help = "name of swimming plot", 
              metavar = "character",
              default = "Swimming_Plot"),
  make_option(c("-w", "--width"), 
              type = "integer",
              help = "width of swimming plot", 
              default = 10), 
  make_option(c("-g", "--height"), 
              type = "integer",
              help = "height of swimming plot", 
              default = 20),   
  make_option(c("-o", "--out"), 
              type = "character",
              help = "output file path", 
              metavar = "character",
              default = "./")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# input parameters
if (!is.null(opt$phyloseq)) {
  phy <- readRDS(opt$phyloseq)  
} else {
  phy <- NULL
}
if (!is.null(opt$phenotype)) {
  phe <- read.csv(opt$phenotype)  
} else {
  phe <- NULL
}
# "Others" contains all the cancer except "gastric cancer", "esophageal carcinoma" and "colon cancer"
cancer <- opt$cancer # "gastric cancer" "esophageal carcinoma" "colon cancer" "Others"
name <- opt$name
Width <- opt$width
Height <- opt$height
out <- opt$out

# phy <- readRDS("dataset/phyloseq/metaphlan2_BJ_all_ps.RDS")
# phe <- read.csv("dataset/phenotype/bjch_sample_sequence_metadata_20220923.csv")
# cancer <- "colon_cancer"
# name <- "Swimming_Plot"
# Width <- 10
# Height <- 20
# out <- "./"


get_plotdata <- function(x1 = phy,
                         x2 = phe,
                         cancer_type = c("all", 
                                         "gastric_cancer", 
                                         "esophageal_carcinoma", 
                                         "colon_cancer", 
                                         "Others"),
                         round_name = c("RoundB", "RoundC", "RoundD", 
                                        "RoundE", "RoundF", "RoundG")) { 
  
  # x1 = phy
  # x2 = phe
  # cancer_type = cancer
  # round_name = c("RoundB", "RoundC", "RoundD",
  #                "RoundE", "RoundF", "RoundG")
  
  columns <- c("SubjectID", "Diagnosis", "Treatment_class", "Response_6",
               paste0("L", c(1:14)), "Therapy_Start_Date", "Sample.Collection.Date",
               "Seq.ID", "Analysis.Project")
  
  if (!is.null(x1)) {
    metadata <- x1@sam_data %>%
      data.frame() %>%
      dplyr::filter(Content == "metagenomic") 
    colnames(metadata)[which(colnames(metadata) == "SeqID")] <- "Seq.ID"
    
    phenotype <- metadata %>% 
       dplyr::select(all_of(columns))
  } else if (!is.null(x2)) {
    phenotype <- x2 %>%
      dplyr::filter(Content == "metagenomic") %>%
      dplyr::select(all_of(columns)) 
  }
  
  # Round & Diagnosis
  pheno <- phenotype %>% 
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
    dplyr::mutate(Diagnosis_new = 
      ifelse(Diagnosis == "gastric cancer", "gastric_cancer",
             ifelse(Diagnosis == "esophageal carcinoma", "esophageal_carcinoma",
                    ifelse(Diagnosis == "colon cancer", "colon_cancer", "Others")
                    )
             )
      )
  
  # filter by round and cancer 
  phen <- pheno %>%
    dplyr::filter(Batch %in% round_name) %>%
    dplyr::filter(Diagnosis_new %in% cancer_type)
  
  # filter by Response_6 and treatment_class
  phen$Treatment_class[which(phen$Treatment_class == "ICI|chemo|target")] <- "ICI|target|chemo"
  phen$Treatment_class[which(phen$Treatment_class == "ICI|targe")] <- "ICI|target"
  treat_type <- c("ICI/安慰剂|target|chemo",
                  "ICI/安慰剂|chemo",
                  "ICI/安慰剂",
                  "ICI|target|chemo", 
                  "ICI|target", 
                  "ICI|chemo",
                  "ICI|FMT",
                  "ICI",
                  "target|chemo",
                  "chemo")
  phen_cln <- phen %>%
    dplyr::filter(Response_6 %in% c("R", "CR", "NR")) %>%
    dplyr::filter(Treatment_class %in% treat_type)
  
  # remove columns with all NULL
  column_drop <- c()
  for (i in 1:ncol(phen_cln)) {
   temp_col <- colnames(phen_cln)[i]
   char_name <- unique(phen_cln[, i])
   if (all(char_name == "")) {
     column_drop <- c(column_drop, temp_col)
   }
  }
  phen_res <- phen_cln %>%
    dplyr::select(-all_of(column_drop))
  
  # Days and Weeks
  phen_res$Day <- as.numeric(as.Date(phen_res$Sample.Collection.Date) - 
                               as.Date(phen_res$Therapy_Start_Date))
  phen_res$Week <- round(phen_res$Day / 7, 2)
  
  # Weeks for evaluation of disease
  phen_eval <- phen_res %>%
    dplyr::select(starts_with("L"))
  phen_eval_res <- data.frame()
  for(i in 1:ncol(phen_eval)) {
    temp_df <- phen_res %>%
      dplyr::select(SubjectID, paste0("L", i)) %>% 
      tidyr::separate(col = paste0("L", i), into = c("Week", "Type")) %>% 
      data.frame()
    
    phen_eval_res <- rbind(phen_eval_res, temp_df)    
  }
  
  df_eval <- phen_eval_res %>%
    dplyr::filter(Week != "",
                  !is.na(Week)) %>%
    dplyr::filter(Type != "",
                  !is.na(Type)) %>%
    dplyr::mutate(Week = gsub("W", "", Week),
                  Week = as.numeric(Week))
  
  columns2 <- c("SubjectID",
                "Diagnosis_new",
                "Treatment_class",
                "Therapy_Start_Date",
                "Seq.ID",
                "Sample.Collection.Date",
                "Response_6")
  merge_week <- phen_res %>%
    dplyr::select(all_of(columns2)) %>%
    dplyr::inner_join(df_eval, by = "SubjectID")
  
  # merge metadata
  mdat <- phen_res %>%
    dplyr::select(all_of(c(columns2, "Week"))) %>%
    dplyr::mutate(Type = "Fecal_Sample") %>%
    rbind(merge_week) %>%
    dplyr::mutate(Size = ifelse(Type == "Fecal_Sample", 4, 5))

  
  return(mdat)
}


get_swimming_plot <- function(x, cancer_type = cancer) {
  
  # colors 
  R6Col <- data.frame(Response_6 = c("NR", "R", "CR"), 
                      ResponseCol = c("#48D1CC", "#32CD32", "#9400D3"))
  type_shape <- c("PR" = 2,
                  "CR" = 6, 
                  "pCR" = 11,
                  "SD" = 5, 
                  "PD" = 0,
                  "Fecal_Sample" = 8, 
                  "ICI/安慰剂|target|chemo" = NA,
                  "ICI/安慰剂|chemo" = NA,
                  "ICI/安慰剂" = NA,
                  "ICI|target|chemo" = NA, 
                  "ICI|target" = NA, 
                  "ICI|chemo" = NA,
                  "ICI|FMT" = NA,
                  "ICI" = NA,
                  "target|chemo" = NA,
                  "chemo" = NA)
  
  type_line <- c("PR" = "blank",
                 "CR" = "blank",
                 "pCR" = "blank",
                 "SD" = "blank",
                 "PD" = "blank",
                 "Fecal_Sample" = "blank", 
                 "ICI/安慰剂|target|chemo" = "solid",
                 "ICI/安慰剂|chemo" = "solid",
                 "ICI/安慰剂" = "solid",
                 "ICI|target|chemo" = "solid", 
                 "ICI|target" = "solid", 
                 "ICI|chemo" = "solid",
                 "ICI|FMT" = "solid",
                 "ICI" = "solid",
                 "target|chemo" = "solid",
                 "chemo" = "solid")
  
  type_col <- c("PR" = "#32CD32",
                "CR" = "#9400D3",
                "pCR" = "#9400D3",
                "SD" = "#48D1CC", 
                "PD" = "#48D1CC",
                "Fecal_Sample" = "#FF0000", 
                "ICI/安慰剂|target|chemo" = "#D51F26",
                "ICI/安慰剂|chemo" = "#272E6A",
                "ICI/安慰剂" = "#208A42",
                "ICI|target|chemo" = "#89288F", 
                "ICI|target" = "#F47D2B", 
                "ICI|chemo" = "#FEE500",
                "ICI|FMT" = "#8A9FD1",
                "ICI" = "#C06CAB",
                "target|chemo" = "#E6C2DC",
                "chemo" = "#90D5E4") 
  
  plotdata <- merge(x, R6Col)
  type_name <- c(unique(plotdata$Type), unique(plotdata$Treatment_class))
  index <- pmatch(type_name, names(type_shape))
  ordered_index <- index[order(index)]
  type_shape_cln <- type_shape[ordered_index]
  type_line_cln <- type_line[ordered_index]
  type_col_cln <- type_col[ordered_index]
  
  ## order plotdata by type*
  plotdata$Type <- factor(plotdata$Type, levels = names(type_col_cln))
  plotdata$SubjectID <- factor(plotdata$SubjectID, levels = unique(plotdata$SubjectID))
  
  # breaks
  week_range <- as.integer(range(plotdata$Week))
  MinWeek <- week_range[1] - 2
  MaxWeek <- week_range[2] + 2
  breaks_seq <- c(as.integer(rev(c(seq(0, MinWeek, by=-10), MinWeek))), 
                  3, 6, 8, 26, 
                  as.integer(seq(50, MaxWeek, by = 20)), MaxWeek)
  
  title_name <- paste0("Swimming Plot [", cancer_type, 
                       " (Patients=", length(unique(plotdata$SubjectID)),") ]")
  
  # y axis label
  y_axis <- plotdata %>%
    dplyr::select(SubjectID, Response_6, ResponseCol) %>%
    dplyr::distinct() %>%
    dplyr::mutate(SubjectID = factor(SubjectID, levels = unique(plotdata$SubjectID)))
  
  pl <- ggplot(plotdata, aes(x = as.integer(Week), y = SubjectID)) +
    geom_point(aes(shape = Type, color = Type), size = plotdata$Size) +
    geom_line(aes(color = Treatment_class)) +
    geom_vline(xintercept = c(0, 3, 6, 8, 26), color = "black", size = 0.6, linetype = "dashed") +
    scale_shape_manual(values = type_shape_cln) +
    scale_color_manual(values = type_col_cln) +
    labs(title = title_name, x = "Week", y = "PatientID") + 
    scale_x_continuous(breaks = breaks_seq,
                       limits = week_range) +
    guides(shape = "none", 
           color = guide_legend(title = "", 
                                override.aes = list(shape = type_shape_cln, 
                                                    linetype = type_line_cln, 
                                                    size = 3))) +  
    theme_bw() +
    theme(plot.title = element_text(face = "bold", color = "black", size = 14, hjust = 0.5),
          axis.title = element_text(face = "bold", color = "black", size = 12),
          axis.text.x = element_text(face = "bold", color = "black", size = 10),
          axis.text.y = ggtext::element_markdown(size = 12, color = y_axis$ResponseCol),
          #axis.text.y = element_text(size = 12, color = y_axis$ResponseCol),
          legend.key.size = unit(2, "line"),
          legend.text = element_text(face = "bold", color = "black", size = 10))
  
  return(pl)
}


# plotdata for swimming plot
plotdata_swim <- get_plotdata(cancer_type = cancer)

# plotting swimming plot
swim_pl <- get_swimming_plot(x = plotdata_swim)

if (!dir.exists(out)) {
  dir.create(out, recursive = TRUE)
}

file_name <- paste0(out, "/", name, "_", cancer, ".pdf")
ggsave(file_name,
       swim_pl,
       width = Width,
       height = Height,
       dpi = 600,
       units = "in")

message('Congratulations, Program Ended Without Problem')
