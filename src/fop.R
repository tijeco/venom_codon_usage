library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(optparse)

option_list = list(
    make_option(c("-b", "--body"),
                type = "character",
                default = NULL,
                help = "stata dataset file name",
                metavar = "character"),
    make_option(c("-v", "--venom"),
                type = "character",
                default = NULL,
                help = "stata dataset file name",
                metavar = "character"),
	  make_option(c("-o", "--out"),
                type = "character",
                default = NULL,
                help = "output file name",
                metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
}

venom_quant <- read.table(opt$venom,header=T)
venom_quant <- subset(venom_quant, TPM > 2)
venom_quant$TPM_log10 <- log10(venom_quant$TPM)
body_quant <- read.table(opt$body,header=T)
body_quant <- subset(body_quant, TPM > 2)
body_quant$TPM_log10 <- log10(body_quant$TPM)

merge_quant <- merge(venom_quant, body_quant, by = "Name",suffixes = c(".venom",".body"))
merge_quant$diff <- merge_quant$TPM_log10.venom - merge_quant$TPM_log10.body
write.csv(merge_quant, file = opt$out)
