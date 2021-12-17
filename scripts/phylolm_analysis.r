library(phylolm)
set.seed(123456)

# load tree
tree_file <- "BAC_faa.tree.nwk" # tree file in newick
tre <- read.tree(tree_file)

tre$edge.length[which(tre$edge.length == 0)] <- 0.00000001
# IA trait01
ia_file <- "g__Bacteroides_genome_type.tsv"
############
# two columns separted with Tab:
# 1. MAGs and 2.their information of Seroconversion coded as "NOT" and "YES"
# mags	groups
# mag1	NOT
############
ia <- read.table(ia_file,head=T,sep="\t",stringsAsFactors = F)
ia$trait01 <- ifelse(ia$groups == 'NOT',0,1)
ia <- subset(ia,select = c(mags,trait01))

# KEGG module
modules <- "g__Bacteroides_genome_t.tsv"
###################################
# the first column is mags, and followings are gene counts in each module
# example:
# mags M00064	M00652	M00797
# mag_name 1	1	0
###################################
mkegg <- read.table(modules,head=T,sep="\t")
row.names(mkegg) <- mkegg$mags

#### module M00064 ####
kegg1 <-subset(mkegg,select = c(mags,M00064))
# dat
dat <- merge(ia,kegg1,by="mags")
row.names(dat) <- dat$mags
max_v <- max(dat$M00064)
dat$M00064 <- dat$M00064/max_v
dat <- subset(dat,select = c(trait01,M00064)) #M00064

fit <- phyloglm(trait01~M00064,phy=tre,data=dat,boot=100,btol = 20)

result1 <- summary(fit)
con1 <- as.data.frame(result1$coefficients)

write.table(dat,"BAC_data.txt",sep="\t",quote = F)
row.names(con1) <- c(paste0("M00064","_intercept"),"M00064")
bac_result <- con1

# KEGG module ### ALL ###
#mkegg <- subset(mkegg,select = -c(mags))
all_kegg <- colnames(mkegg)

# max
max_value  <- apply(mkegg,2,max)

#length(all_kegg)
for (i in 2:length(all_kegg)){
  mod <- all_kegg[i]
  kegg_tmp <-subset(mkegg,select = c("mags",mod))
  if (max_value[i] == 0) {
    
    next;
  }else{
  kegg_tmp[,mod] <- kegg_tmp[,mod]/as.integer(max_value[i])
  # dat
  dat_tmp <- merge(ia,kegg_tmp,by="mags")
  row.names(dat_tmp) <- dat_tmp$mags
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

outfile <- "BAC_phyloglm_all_max_btol20.txt"
write.table(bac_result,outfile,sep="\t",quote=F)
