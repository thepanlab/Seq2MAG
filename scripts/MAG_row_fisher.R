#!/usr/bin/env Rscript

# fisher test on each row
row_fisher <- function(row, alt = 'greater', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           p_val = f$p.value,
           or = f$estimate[[1]]
           #or_ll = f$conf.int[1],
           #or_ul = f$conf.int[2]
))
}

# generate sample data
data_in <- read.table("MAG_sum_out_E-6_coeff_Ncluster.tsv",head=T,sep="\t",stringsAsFactors=F)
row.names(data_in) <- data_in$MAG_name
data_in$background_total <- 64142
data_in$background_up <- 2355
data_in$background_down <- 2504
analysis_data <- subset(data_in, select=c(Cup_core,background_up,core_Cnumbers,background_total))
#analysis_data <- subset(data_in, select=c(Cdown_core,background_down,core_Cnumbers,background_total))

#a: MAG_up; b: background_up; c: MAG_total; D: background_total

# run
p <- data.frame(t(apply(analysis_data, 1, row_fisher)))
p$qvalue <- p.adjust(p$p_val,method="BH")
final_p <- subset(p,p$qvalue < 0.01)
final_p$MAG_name <- row.names(final_p)

# add taxonomy infor
taxonomy <- read.table("High_MAGs_gtdbtk.ar_bac.summary.tsv",head=T,sep="\t",stringsAsFactors=F)
merge_p <- merge(final_p,taxonomy,by="MAG_name")

write.table(merge_p,"up_enrichment.tsv",sep="\t")
