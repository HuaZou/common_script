## script contain 3 steps
### 1. filter data and coverage 
### 2. test result 
### 3. reportscore

args <- commandArgs(T)
if(length(args) < 10){
  stop("
       Rscript main.R [phen] [pro] [list]  [info] [cfg] [cfg2] [cfg3] [n] [dir] [type] [prefix]
       phen : phentype
       pro  : profile
       list : All_map_ko.list/All_module_ko.list
       info : map_KO_info.list.cleanonly.v2/Module_KO_info.list.cleanKOonly.v2
       cfg  : column to test, must be two factors
       cfg2 : column to split data for timescale
       cfg3 : column to connect profile colnames to phenotype 
       n : number which to filter profile [one row must have more than n times in one factor]
       dir : output dir
       type: pathway or module
       prefix : output file prefix
       ")	
}
######## load data #########################################

## load data 
phe <- read.csv("data/n50.phenotype.csv", header=T)
pro <- read.table("data/KO.profile", row.names=1, header=T)
ampk <- read.table("data/All_map_ko.list", row.names=1)
amdk <- read.table("data/All_module_ko.list", row.names=1)
mpki <- read.table("data/map_KO_info.list.cleanonly.v2", row.names = 1)
mdki <- read.table("data/Module_KO_info.list.cleanKOonly.v2", row.names = 1)
info <- read.table("data/ModuleAndPathway.tab", sep = "\t")
mrm <- read.table("data/MicroRelatedModule.tab", sep = "\t")
mrp <- read.table("data/MicroRelatedPathway.tab", sep="\t")
cfg <- "Group"
cfg1 <- "Baseline"
cfg2 <- "SampleID"
dir <- "./"
type <- "module"
prefix <- "base"
paired <- "FALSE"
method <- "wilcox"

### step1 ########################
source("step1.R")

### step2 ########################
source("step2.R")

### step3 #######################
source("step3.R")
