suppressPackageStartupMessages({ 
  library(dplyr)
  library(tibble)
  library(optparse)
  library(data.table)
  library(APAlyzer)
  library(TBX20BamSubset)
  library(Rsamtools)
})


option_list <- list(
  make_option(c("-b", "--bam"), 
              type = "character",
              help = "bam csv file (1st column: sampleID; 2nd: bam path)", 
              metavar = "character"),
  make_option(c("-r", "--reference"), 
              type = "character", # RData/gtf
              help = "genomic reference type", 
              metavar = "character"),    
  make_option(c("-g", "--genome"), 
              type = "character",
              help = "genomic reference file", 
              metavar = "character"), 
  make_option(c("-c", "--chromosome"), 
              type = "character",
              default = "all", # chr19
              help = "chromosome to be selected", 
              metavar = "character"),  
  make_option(c("-e", "--expression"), 
              type = "character", 
              default = "all", # 3UTR/IPA
              help = "APA expression: 3UTR and intronic APA", 
              metavar = "character"),  
  make_option(c("-o", "--out"), 
              type = "character",
              help = "output file path", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# input parameters
bam_path <- opt$bam
ref_type <- opt$reference
ref_path <- opt$genome
chrom <- opt$chromosome
expr_type <- opt$expression
dir <- opt$out


# bam_path <- "bam_file.tsv"
# ref_type <- "RData"
# ref_path <- "mm9_REF.RData"
# chrom <- "chr19"
# expr_type <- "3UTR"
# dir <- "result"


# step1: bam file
bam_vector <- read.table("bam_file.tsv", header = TRUE)
bam_file <- bam_vector$BamPath
names(bam_file) <- bam_vector$SampleID

# step2: genomic reference
if (ref_type == "RData") {
  # data from built reference
  require(repmis)
  URL <- "https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
  source_data(paste0(URL, ref_path, "?raw=True"))
  
  if (ref_path == "mm9_REF.RData") {
    refUTRraw_temp <- refUTRraw
    dfIPAraw_temp <- dfIPA
    dfLEraw_temp <- dfLE
  } else if (ref_path == "hg19_REF.RData") {
    refUTRraw_temp <- refUTRraw_hg19
    dfIPAraw_temp <- dfIPA_hg19
    dfLEraw_temp <- dfLE_hg19
  }
  
} else if (ref_type == "gtf") {
  # building reference from gtf file
  PASREFraw <- PAS2GEF(ref_path)  
  
  refUTRraw_temp <- PASREFraw$refUTRraw
  dfIPAraw_temp <- PASREFraw$dfIPA
  dfLEraw_temp <- PASREFraw$dfLE
}

# step3: whether to choose chromosome
if (chrom == "all") {
  UTRdbraw <- refUTRraw_temp
  dfIPAraw <- dfIPAraw_temp
  dfLEraw <- dfLEraw_temp   
} else {
  # multiple chromosome or not
  if (length(grep(":", chrom)) > 0) {
    chroms <- unlist(strsplit(chrom, ":"))
  } else {
    chroms <- chrom
  }
  UTRdbraw <- refUTRraw_temp[which(refUTRraw_temp$Chrom %in% chroms), ]
  dfIPAraw <- dfIPAraw_temp[which(dfIPAraw_temp$Chrom %in% chroms), ]
  dfLEraw <- dfLEraw_temp[which(dfLEraw_temp$Chrom %in% chroms), ]
}
## aUTR cUTR
PASREF_temp <- REF4PAS(UTRdbraw, dfIPAraw, dfLEraw)
UTRdb <- PASREF_temp$UTRdbraw
dfIPA <- PASREF_temp$dfIPA
dfLE <- PASREF_temp$dfLE  

# step4: APA expression (3UTR and IPA)
if (expr_type == "all") {
  # 3UTR
  UTR_APA_OUT <- PASEXP_3UTR(UTRdb, bam_file, Strandtype = "forward")
  # IPA
  IPA_OUT <- PASEXP_IPA(dfIPA, dfLE, bam_file, Strandtype = "invert", nts = 4)
  
  final_OUT <- list(UTR = UTR_APA_OUT,
                    IPA = IPA_OUT)
} else if (expr_type == "3UTR") { 
  # 3UTR
  final_OUT <- PASEXP_3UTR(UTRdb, bam_file, Strandtype = "forward")  
} else if (expr_type == "IPA") { 
  final_OUT <- PASEXP_IPA(dfIPA, dfLE, bam_file, Strandtype = "invert", nts = 4)
}

# step5: output
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

if (!is.data.frame(final_OUT)) {
  file_name <- paste0(dir, "/APA_Expr_", expr_type, ".RDS")
  saveRDS(final_OUT, file_name, compress = TRUE)
} else {
  file_name <- paste0(dir, "/APA_Expr_", expr_type, ".tsv")
  write.table(final_OUT, file_name, quote = F, row.names = F, sep = "\t")
}

print("Program Ended without Problems")
