library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(optparse)

option_list = list(
    make_option(c("-d", "--data"),
                type = "character",
                default = NULL,
                help = "stata dataset file name",
                metavar = "character"),
	make_option(c("-o", "--out"),
                type = "character",
                default = NULL,
                help = "output file name",
                metavar = "character"),
  make_option(c("-f", "--figure"),
                type = "character",
                default = NULL,
                help = "figure file name",
                metavar = "character"),
  make_option(c("-r", "--rscu"),
                type = "character",
                default = NULL,
                help = "rscu file name",
                metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$data)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
}

rscu_file <- read.csv(opt$data)
rscu_file <- rscu_file  %>% filter(aa != "Stop")

aa_table <- unique(rscu_file[c("aa","codon")])
deltaRSCU <- rscu_file %>% group_by(codon) %>% do(tidy(t.test(rscu~class, data =.,alternative = "less")))
deltaRSCU <- inner_join(aa_table,deltaRSCU)
write.csv(deltaRSCU,opt$rscu, row.names = F, quote = F)

optimized_codons <- deltaRSCU %>% group_by(aa) %>% filter(p.value == min(p.value)) %>% filter(p.value < 0.05)
write.csv(optimized_codons,opt$out, row.names = F, quote = F)

rscu_file$optimal <- rscu_file$codon %in% optimized_codons$codon
rscu_plot <- ggplot(rscu_file, aes(x=codon, y =rscu,fill = class,col = optimal))
rscu_plot <- rscu_plot + geom_boxplot() + facet_wrap(~ aa,scales="free_x") +theme_bw()+theme(legend.position = "top") +
  scale_color_manual(values=c("black","red")) +
  scale_fill_manual(values = c("#646464","#969696"))

print(rscu_plot)
ggsave(opt$figure, width = 12, height = 9, dpi = 240)
