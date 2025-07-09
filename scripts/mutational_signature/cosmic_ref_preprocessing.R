signatures.tmp=read.delim("ref/COSMIC_v3.2_SBS_GRCh37.txt",sep="\t",stringsAsFactors = F)
signatures.cosmicv32=data.frame(t(signatures.tmp[,-1]))
colnames(signatures.cosmicv32)=signatures.tmp$Type
signatures.cosmicv32
signatures.cosmicv32 = signatures.cosmicv32[,colnames(signatures.cosmic)]
saveRDS(signatures.cosmicv32, file = "ref/COSMIC_v3.2_SBS_GRCh37.rds")


signatures.tmp=read.delim("ref/COSMIC_v3.2_DBS_GRCh37.txt",sep="\t",stringsAsFactors = F)
signatures.cosmicv32=data.frame(t(signatures.tmp[,-1]))
colnames(signatures.cosmicv32)=signatures.tmp$Type
saveRDS(signatures.cosmicv32, file = "ref/COSMIC_v3.2_DBS_GRCh37.rds")


signatures.tmp=read.delim("ref/COSMIC_v3.2_ID_GRCh37.txt",sep = "\t",stringsAsFactors = F)
signatures.cosmicv32=data.frame(t(signatures.tmp[,-1]))
colnames(signatures.cosmicv32)=signatures.tmp$Type
saveRDS(signatures.cosmicv32, file = "ref/COSMIC_v3.2_ID_GRCh37.rds")

