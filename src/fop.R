library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(optparse)

option_list = list(
    make_option(c("-q", "--quant"),
                type = "character",
                default = NULL,
                help = "stata dataset file name",
                metavar = "character"),
    make_option(c("-f", "--fop"),
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

merge_quant <- read.csv(opt$quant)

venom_up_1000TPM <- merge_quant %>%
  select(Name,TPM.venom,TPM.body,diff) %>%
  filter(diff >1 & TPM.venom >1000)
venom_up_1000TPM$tissue <- "venom"

body_up_1000TPM <- merge_quant %>%
  select(Name,TPM.venom,TPM.body,diff) %>%
  filter(diff < -1 & TPM.body >1000)
body_up_1000TPM$tissue <- "body"

combined_1000TPM <- merge(venom_up_1000TPM,body_up_1000TPM, all = T)

fop <- read.csv(opt$fop)

combined_1000TPM_fop <- merge(fop, combined)

write.csv(combined_1000TPM_fop, file = opt$out, row.names = F, quote = F)
