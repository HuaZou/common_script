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
  make_option(c("-a", "--phyloseq1"), 
              type = "character",
              help = "phyloseq from all rounds", 
              metavar = "character",
              default = "Combined_round/phyloseq/metaphlan2_BJ_all_ps.RDS"),
  make_option(c("-p", "--phyloseq2"), 
              type = "character",
              help = "phyloseq from all rounds with R6(R and NR) and baseline", 
              metavar = "character",
              default = "Combined_round/phyloseq/metaphlan2_BJ_all_ps_baseline_R6.RDS"),  
  make_option(c("-n", "--name"), 
              type = "character",
              help = "name of file", 
              metavar = "character",
              default = "Summarize_Patients_Samples"),  
  make_option(c("-o", "--out"), 
              type = "character",
              help = "output file path", 
              metavar = "character",
              default = "./")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# input parameters
phy1 <- readRDS(opt$phyloseq1)
phy2 <- readRDS(opt$phyloseq2)
name <- opt$name
out <- opt$out

# phy1 <- readRDS("Combined_round/phyloseq/metaphlan2_BJ_all_ps.RDS")
# phy2 <- readRDS("Combined_round/phyloseq/metaphlan2_BJ_all_ps_baseline_R6.RDS")
# name <- "Summarize_Patients_Samples"
# out <- "./"

all_round <- c("RoundCR",
               "RoundB", "RoundC", "RoundD", 
               "RoundE", "RoundF", "RoundG")

get_account <- function(round_type,
                        ps) {
  
  # ps = phy1
  # round_type = "RoundCR"
  
  metadata <- ps@sam_data %>%
    data.frame()
  
  phen <- metadata %>%
    dplyr::filter(Batch == round_type)
  
  data_sample <- nrow(phen)
  data_patient <- length(unique(phen$SubjectID))
  
  res <- data.frame(Round = round_type,
                    Samples = data_sample,
                    Patients = data_patient)
  
  return(res)
}

# origin
phy1_origin <- do.call(rbind, lapply(all_round, get_account, phy1))
# baseline with R6
phy2_baseR6 <- do.call(rbind, lapply(all_round, get_account, phy2))

# merge
colnames(phy1_origin) <- c("Round", "Samples_Origin", "Patients_Origin")
colnames(phy2_baseR6) <- c("Round", "Samples_baseR6", "Patients_baseR6")
phy_res <- dplyr::inner_join(phy1_origin, 
                             phy2_baseR6,
                             by = "Round") %>%
  dplyr::select(-Patients_baseR6)


plotdata <- phy_res %>%
  tidyr::gather(key = "kind", value = "count", -Round) %>%
  dplyr::mutate(Round = factor(Round, levels = all_round),
                kind = factor(kind, levels = c("Samples_Origin", 
                                               "Patients_Origin", 
                                               "Samples_baseR6")))

# barplot
pl_box <- ggplot(data = plotdata, aes(x = kind, y = count, fill = Round)) +
        geom_bar(stat = "identity", position="dodge", color = "black") +
        geom_text(aes(label = count), 
                  position = position_dodge(0.9), 
                  vjust = 1.6, size = 4)+
        labs(x = "", y = "Number") +
        #scale_y_continuous(expand = c(0, 0)) +
        guides(fill = guide_legend(bycol = TRUE, override.aes = list(size = 2))) +
        theme_bw()+
        theme(axis.title.y = element_text(face = "bold", color = "black", size = 14),
              axis.title.x = element_text(face = "bold", color = "black", size = 14, vjust = -1.2),
              axis.text.y = element_text(face = "bold", color = "black", size = 10),
              axis.text.x = element_text(face = "bold", color = "black", size = 12,
                                         angle = 45, vjust = 0.5),
              panel.grid = element_blank(),
              legend.position = "right",
              legend.key.height = unit(0.6, "cm"),
              legend.text = element_text(face = "bold", color = "black", size = 10))

# count table
phy_res_sum <- data.frame(Round="Sum",
                          Samples_Origin = sum(phy_res$Samples_Origin),
                          Patients_Origin = sum(phy_res$Patients_Origin),
                          Samples_baseR6 = sum(phy_res$Samples_baseR6))
phy_final <- rbind(phy_res, phy_res_sum)
count_tab <- ggtexttable(phy_final, 
                         rows = NULL, 
                         theme = ttheme("mOrange"))

# text
text_str <- paste("Patients and Samples from B to G round",
              "are in columns of Samples_origin and Patients_origin.",
              "Choosing patients whose baseline samples are closest to 1st therapy date", 
              " and have the Response_6 with R or NR.", sep = " ")
pl_text <- ggparagraph(text = text_str, face = "italic", size = 11, color = "black")

# merge boxplot and table
pl <- ggarrange(pl_box, count_tab, pl_text,
                ncol = 1, nrow = 3,
                heights = c(1, 0.5, 0.1))

png_name <- paste0(out, "/", name, ".png")
ggsave(png_name,
       pl, width = 10, height = 8, dpi = 600,
       units = "in")

pdf_name <- paste0(out, "/", name, ".pdf")
ggsave(pdf_name,
       pl, width = 10, height = 8, dpi = 600,
       units = "in")

csv_name <- paste0(out, "/", name, ".csv")
write.csv(phy_final, csv_name, row.names = F)

message('Congratulations, Program Ended Without Problem')
