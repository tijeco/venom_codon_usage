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
    make_option(c("-v", "--violin"),
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
fop <- read.csv(opt$fop)

# venom_up_1000TPM <- merge_quant %>%
#   select(Name,TPM.venom,TPM.body,diff) %>%
#   filter(diff >1 & TPM.venom >1000)
# venom_up_1000TPM$tissue <- "venom"
#
# body_up_1000TPM <- merge_quant %>%
#   select(Name,TPM.venom,TPM.body,diff) %>%
#   filter(diff < -1 & TPM.body >1000)
# body_up_1000TPM$tissue <- "body"
#
# combined_1000TPM <- merge(venom_up_1000TPM,body_up_1000TPM, all = T)



# combined_1000TPM_fop <- merge(fop, combined_1000TPM)

# write.csv(combined_1000TPM_fop, file = opt$out, row.names = F, quote = F)

venom_up_5percent <- merge_quant %>%
  select(Name,TPM.venom,TPM.body,diff) %>%
  filter(diff > 1 & TPM.venom > quantile(TPM.venom, 0.95))

body_up_5percent <- merge_quant %>%
  select(Name,TPM.venom,TPM.body,diff) %>%
  filter(diff < -1 & TPM.body > quantile(TPM.body, 0.95))

venom_up_5percent$tissue <- "venom"
body_up_5percent$tissue <- "body"

combined_5percent <- merge(venom_up_5percent,body_up_5percent, all = T)
combined_5percent_fop <- merge(fop, combined_5percent)
write.csv(combined_5percent_fop, file = opt$out, row.names = F, quote = F)

violin_plot <- ggplot(Cc_5percent_fop_combined, aes(x = tissue, y = fop,fill = tissue))
violin_plot + geom_violin() + geom_boxplot(width=0.1,fill = "white")

print(violin_plot)
ggsave(opt$violin, width = 9, height = 5, dpi = 240)
