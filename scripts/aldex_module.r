library(ALDEx2)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

# load input
ko.term_tmp = read.delim(args[4], sep="\t", row.names=1, check.names=FALSE,stringsAsFactors = F) # KEGG module summary counts

remove_file = read.table("remove_list.txt",head=F,stringsAsFactors = F) # remove the modules from plants, animals and fungi
remove_list = remove_file$V1

ko.term = ko.term_tmp[!row.names(ko.term_tmp) %in% remove_list,]


row_name = row.names(ko.term)
function_name = ko.term[,1]
ko.names = data.frame(row_name,function_name)
names(ko.names) = c("modules","function_des")
row.names(ko.names) = ko.names$modules

ko.dset_tmp = ko.term[,2:ncol(ko.term)]
ko.dset = ko.dset_tmp[rowSums(ko.dset_tmp[])>0,]

all.genome = read.delim(args[1], row.names=1) 

# define genomes to analyze
sig.genomes = rownames(all.genome)[which(all.genome$groups == 'SIG')]
not.genomes = rownames(all.genome)[which(all.genome$groups == 'NOT')]

# prepare dataset
analy.dset = ko.dset[,c(sig.genomes, not.genomes)]
analy.dset[is.na(analy.dset)] = 0 # replace NAs with 0
analy.conds = c(rep("SIG",length(sig.genomes)),rep("NOT", length(not.genomes)))

# perform aldex analysis
aldex.analy = aldex.clr(reads=analy.dset, conds=analy.conds, mc.samples=128)
aldex.eff = aldex.effect(aldex.analy, useMC=TRUE)
aldex.res = aldex.ttest(aldex.analy)
res.all = data.frame(rownames(aldex.eff), aldex.eff,aldex.res)

# summarize results and get most significant
effect_thresh = 0.5
res.sig = res.all[which(res.all$wi.eBH<0.01),]
res.eff = res.sig[which(abs(res.sig$effect) > effect_thresh),]
res.plot = merge(res.eff, ko.names, by="row.names")
write.table(res.plot,args[2],sep="\t",quote=F)
