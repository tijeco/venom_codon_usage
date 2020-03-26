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
    make_option(c("-f", "--sc"),
                type = "character",
                default = NULL,
                help = "stata dataset file name",
                metavar = "character"),
    make_option(c("-s", "--test"),
                type = "character",
                default = NULL,
                help = "stata dataset file name",
                metavar = "character"),
    make_option(c("-m", "--mean"),
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
sc <- read.csv(opt$sc)

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



# combined_1000TPM_sc <- merge(sc, combined_1000TPM)

# write.csv(combined_1000TPM_sc, file = opt$out, row.names = F, quote = F)

venom_up_5percent <- merge_quant %>%
  select(Name,TPM.venom,TPM.body,diff) %>%
  filter(diff > 1 & TPM.venom > quantile(TPM.venom, 0.95))

body_up_5percent <- merge_quant %>%
  select(Name,TPM.venom,TPM.body,diff) %>%
  filter(diff < -1 & TPM.body > quantile(TPM.body, 0.95))

venom_up_5percent$tissue <- "venom"
body_up_5percent$tissue <- "body"

combined_5percent <- merge(venom_up_5percent,body_up_5percent, all = T)
combined_5percent_sc <- merge(sc, combined_5percent)
write.csv(combined_5percent_sc, file = opt$out, row.names = F, quote = F)
#
# violin_plot <- ggplot(combined_5percent_sc, aes(x = tissue, y = total_sc,fill = tissue))
# violin_plot + geom_violin() + geom_boxplot(width=0.1,fill = "white")
#
# print(violin_plot)
# ggsave(opt$violin, width = 9, height = 5, dpi = 240)

total_t_test <- tidy(t.test(combined_5percent_sc$total_sc~combined_5percent_sc$tissue))
write.csv(total_t_test, file = opt$test, row.names = F, quote = F)

mean_t_test <- tidy(t.test(combined_5percent_sc$mean_sc~combined_5percent_sc$tissue))
write.csv(mean_t_test, file = opt$mean, row.names = F, quote = F)
