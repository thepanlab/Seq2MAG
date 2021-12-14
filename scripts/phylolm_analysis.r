# Rscript  FUS  g__Fusicatenibacter

set.seed(123456)
library(phylolm)

args <- commandArgs(trailingOnly = TRUE)

args[1] <- 'BAC'
args[2] <- 'g__Bacteroides'
dir <- paste0("/home/lizhang12/phylolm_test/",args[1])
setwd(dir)

# load tree
tree_file <- paste0(args[1], "_faa.tree.nwk")
tre <- read.tree(tree_file)

tre$edge.length[which(tre$edge.length == 0)] <- 0.00000001
# IA trait01
ia_file <- paste0(args[2], "_genome_type.tsv")
ia <- read.table(ia_file,head=T,sep="\t",stringsAsFactors = F)
ia$module <- gsub("^","GMS_",ia$module)
ia$trait01 <- ifelse(ia$groups == 'NOT',0,1)
ia <- subset(ia,select = c(module,trait01))

# KEGG module
modules <- paste0(args[2], "_genome_t.tsv")
mkegg <- read.table(modules,head=T,sep="\t")
mkegg$module <- gsub("^","GMS_",mkegg$module)
row.names(mkegg) <- mkegg$module



#### module M00570 ####
kegg1 <-subset(mkegg,select = c(module,M00064))
# dat
dat <- merge(ia,kegg1,by="module")
row.names(dat) <- dat$module
max_v <- max(dat$M00064)
dat$M00064 <- dat$M00064/max_v
dat <- subset(dat,select = c(trait01,M00064)) #M00024

fit <- phyloglm(trait01~M00064,phy=tre,data=dat,boot=100,btol = 20) #,

result1 <- summary(fit)
con1 <- as.data.frame(result1$coefficients)

write.table(dat,"BAC_data.txt",sep="\t",quote = F)
row.names(con1) <- c(paste0("M00570","_intercept"),"M00570")
bac_result <- con1

# KEGG module ### ALL ###
#mkegg <- subset(mkegg,select = -c(module))
all_kegg <- colnames(mkegg)

# max
max_value  <- apply(mkegg,2,max)

#length(all_kegg)
for (i in 2:length(all_kegg)){
  mod <- all_kegg[i]
  kegg_tmp <-subset(mkegg,select = c("module",mod))
  if (max_value[i] == 0) {
    
    next;
  }else{
  kegg_tmp[,mod] <- kegg_tmp[,mod]/as.integer(max_value[i])
  # dat
  dat_tmp <- merge(ia,kegg_tmp,by="module")
  row.names(dat_tmp) <- dat_tmp$module
  dat_tmp <- subset(dat_tmp,select = c("trait01",mod))
  names(dat_tmp) <- c("trait01","module_name")
  
  # fit
  fit_tmp <- phyloglm(trait01~module_name,phy=tre,data=dat_tmp,boot=100,btol = 20)
  
  # summary
  result_tmp <- summary(fit_tmp)
  con_tmp <- as.data.frame(result_tmp$coefficients)
  
  row.names(con_tmp) <- c(paste0(mod,"_intercept"),mod)
  bac_result <- rbind(bac_result,con_tmp)
}}

outfile <- paste0(args[1], "_phyloglm_all_max_btol20.txt")
write.table(bac_result,outfile,sep="\t",quote=F)