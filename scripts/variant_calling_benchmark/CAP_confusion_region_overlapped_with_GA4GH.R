library(GenomicRanges)
all_bed=list.files("GRCh37@all/",recursive=T)
all_bed=all_bed[!grepl(all_bed,pattern="GenomeSpecific")]
all_bed=all_bed[grepl(all_bed,pattern="bed.gz")]
A = read.table("CAP_confusion_region.bed", header = FALSE, stringsAsFactors = FALSE)
A_gr = GRanges(seqnames = A$V1, ranges = IRanges(start = A$V2, end = A$V3))
A_gr$overlapped=""
for(i in all_bed){
	B = read.table(paste0("GRCh37@all/",i), header = FALSE, stringsAsFactors = FALSE)
	B_gr = GRanges(seqnames = B$V1, ranges = IRanges(start = B$V2, end = B$V3))
	overlaps = findOverlaps(A_gr, B_gr)
	if(length(overlaps)>0){
		A_gr$overlapped[overlaps@from]=paste0(A_gr$overlapped[overlaps@from],",",i)
	}
}

aa=t(sapply(1:length(A_gr),function(x){
	t(sapply(all_bed,function(y){
		grepl(A_gr$overlapped[x],pattern=y)
	}))
}))
colnames(aa)=all_bed
result=cbind(as.data.frame(A_gr),aa)
write.table(result,"CAP_confusion_region_overlapped_with_GA4GH.txt",sep="\t",row.names=F)

write.table(result[,c(1:6,which(colSums(result[,-c(1:6)])!=0)+6)],"CAP_confusion_region_filtered.txt",sep="\t",row.names=F)
