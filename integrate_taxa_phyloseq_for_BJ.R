suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(phyloseq))

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

option_list <- list(
  make_option(c("-p", "--phylum"), 
              type = "character",
              help = "phylum phyloseq", 
              metavar = "character"),
  make_option(c("-c", "--class"), 
              type = "character",
              help = "class phyloseq", 
              metavar = "character"),
  make_option(c("-o", "--order"), 
              type = "character",
              help = "order phyloseq", 
              metavar = "character"),
  make_option(c("-f", "--family"), 
              type = "character",
              help = "family phyloseq", 
              metavar = "character"), 
  make_option(c("-g", "--genus"), 
              type = "character",
              help = "genus phyloseq", 
              metavar = "character"), 
  make_option(c("-s", "--species"), 
              type = "character",
              help = "species phyloseq", 
              metavar = "character"), 
  make_option(c("-n", "--name"), 
              type = "character",
              help = "name of list file", 
              metavar = "character"),
  make_option(c("-r", "--rounds"), 
              type = "character",
              help = "name of rounds", 
              metavar = "character",
              default = "all"),  
  make_option(c("-d", "--out"), 
              type = "character",
              help = "output file path", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list,
                           add_help_option = FALSE)
opt <- parse_args(opt_parser)

# input parameters
taxa_phylum <- readRDS(opt$phylum)
taxa_class <- readRDS(opt$class)
taxa_order <- readRDS(opt$order)
taxa_family <- readRDS(opt$family)
taxa_genus <- readRDS(opt$genus)
taxa_species <- readRDS(opt$species)
rounds <- opt$rounds
name <- opt$name
out <- opt$out

# rounds <- "all"
# name <- "metaphlan2"
# out <- "phyloseq_v2"

ps_list <- list(phylum = taxa_phylum,
                class = taxa_class,
                order = taxa_order,
                family = taxa_family,
                genus = taxa_genus,
                species = taxa_species)

# output 
if (!dir.exists(out)) {
  dir.create(out, recursive = T)
}

round_names <- paste0("metaphlan2_BJ_", rounds, "_ps")

if (name == "total") {
  ps_filename <- paste0(out, "/", round_names, "_total_list.RDS")  
} else if (name == "filter") {
  ps_filename <- paste0(out, "/", round_names, "_filter_list.RDS")  
} else if (name == "baseR6") {
  ps_filename <- paste0(out, "/", round_names, "_baseR6_list.RDS")  
} else if (name == "total_norm") {
  ps_filename <- paste0(out, "/", round_names, "_total_norm_list.RDS")  
} else if (name == "filter_norm") {
  ps_filename <- paste0(out, "/", round_names, "_filter_norm_list.RDS")  
} else if (name == "baseR6_norm") {
  ps_filename <- paste0(out, "/", round_names, "_baseR6_norm_list.RDS")  
}

saveRDS(ps_list, ps_filename, compress = TRUE)

message('Congratulations, Program Ended Without Problem')
