#!/usr/bin/R

#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This program Generate origin python/R/perl script with basic funciton      #
# Args :  1 args                                                             #
#   Type:	which type of script to generate								 #
#                                                                            #
# Output: 1 script		                                                     #
#----------------------------------------------------------------------------#


if (!require(pacman)) {
	install.packages("pacman", dependencies=T)
}
pacman::p_load(dplyr,ggplot2)

args <- commandArgs(T)

if (length(args) < num) {
	stop("Usage:")
}

args1 <- read.csv(args[1])
args2 <- args[2]

fun <- function(args){
	statements}

