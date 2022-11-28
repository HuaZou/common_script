suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
  make_option(c("-f", "--folder"), 
              type = "character",
              help = "file path", 
              metavar = "character"),
  make_option(c("-m", "--md5sum"), 
              type = "character",
              help = "previous md5sum file", 
              metavar = "character"),  
  make_option(c("-p", "--pattern"), 
              type = "character", 
              default = "Gmetadata.*.csv",
              help = "files with specific pattern", 
              metavar = "character"),  
  make_option(c("-o", "--out"), 
              type = "character",
              help = "output file path", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# input parameters
fold_path <- opt$folder
prev_md5 <- opt$md5sum
file_pattern <- opt$pattern
dir <- opt$out

# fold_path <- "./"
# prev_md5 <- "Gmetadata_from_Doctor_csv_md5sum.csv" # NULL
# file_pattern <- "csv"
# dir <- "./"

# previous md5
if (prev_md5 == "NULL") {
  prev_md5 <- NULL 
}

if (!is.null(prev_md5)) {
  previous_md5 <- data.table::fread(prev_md5)  
} else {
  previous_md5 <- NULL
}

# files from doctor to md5 checking
files <- list.files(path = fold_path,
                    pattern = paste0("Gmetadata.*.", file_pattern))
if (!is.null(prev_md5)) {
  prev_md5_cln <- gsub("\\S+\\/", "", prev_md5)
  files_cln <- files[-which(files == prev_md5_cln)]  
} else {
  files_cln <- files
}

# calculate files' md5sum
md5sum_res <- tools::md5sum(files = files_cln)
current_md5 <- data.frame(
                  file = names(md5sum_res),
                  md5sum = as.character(md5sum_res),
                  calculate_date = as.character(Sys.Date()))

# whether the files have been changed according to compare the previous and current md5sum
if (is.null(previous_md5)) {
  res <- current_md5
} else {
  temp_res <- data.frame()
  for (i in 1:nrow(previous_md5)) {
    
    pre_file <- previous_md5$file[i]
    pre_md5 <- previous_md5$md5sum[i]
    pre_date <- previous_md5$calculate_date[i]
    
    for (j in 1:nrow(current_md5)) {
      
      cur_file <- current_md5$file[j]
      cur_md5 <- current_md5$md5sum[j]
      cur_date <- current_md5$calculate_date[j]
      
      if (pre_file == cur_file) {
        if (pre_md5 == cur_md5) {
          temp <- data.frame(file = pre_file,
                             md5sum = pre_md5,
                             calculate_date = pre_date)
        } else if (pre_md5 != cur_md5) {
          temp <- rbind(data.frame(file = pre_file,
                                   md5sum = pre_md5,
                                   calculate_date = pre_date),
                        data.frame(file = cur_file,
                                   md5sum = cur_md5,
                                   calculate_date = cur_date))        
        }
      } else {
        temp <- data.frame(file = cur_file,
                           md5sum = cur_md5,
                           calculate_date = cur_date)      
      }
      
      temp_res <- rbind(temp_res, temp)
    }
  }
  
  res <- temp_res %>% 
    dplyr::distinct() %>%
    dplyr::arrange(file, desc(calculate_date))
}

res_final <- res[pmatch(unique(res$md5sum), res$md5sum), , F]

if (!dir.exists(dir)) {
  dir.create(dir)
}
res_file_path <- paste0(dir, "/Gmetadata_from_Doctor_csv_md5sum.csv")
write.csv(res_final, res_file_path, row.names = F)

message('Congratulations, Program Ended Without Problem')
