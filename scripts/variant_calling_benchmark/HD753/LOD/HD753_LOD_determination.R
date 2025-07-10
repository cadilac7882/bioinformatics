library(data.table)
TP=fread("../TP_hotspot/HD753_filter_by_intersect_bed_processed.vcf.avinput",sep="\t",data.table=F)
result=TP[,c(1:5,8)]
names(result)=c("chr","start","end","ref","alt","TP_DP")

## combine all testing samples
sampleList=read.delim("../testing/sampleList.txt",header=F)[,1]
for(i in sampleList){
	tmp=fread(paste0("../testing/",i,"/main.vcf_exclude_CNV.vcf.gz_filter_by_intersect_bed.vcf.avinput"),sep="\t",data.table=F)
	tmp=tmp[,c(1:5,15)]
	names(tmp)=c("chr","start","end","ref","alt","INFO")
	
	tmp[,c("DP","VAF")]=t(sapply(1:nrow(tmp),function(x){unlist(strsplit(tmp$INFO[x],split=":"))[c(4,5)]}))
	tmp$DP=as.numeric(tmp$DP)
	tmp$VAF=as.numeric(tmp$VAF)
	tmp=tmp[,-6]
	names(tmp)[-c(1:5)]=paste0(i,c("_DP","_VAF"))
	
	result=merge(result,tmp,by=c("chr","start","end","ref","alt"),all=T)
}

result$type=NA
result$type[which(result$ref=="-"|result$alt=="-")]="INDEL"
result$type[which(result$ref%in%c("A","T","C","G")&result$alt%in%c("A","T","C","G"))]="SNV"
#result$type[is.na(result$type)]="MNV"

## remove undetermined variant type
result=result[!is.na(result$type),]

## statistics for variants passing hard-filtering
result$PASS_Detection_rate=rowSums(!is.na(result[,names(result)[grep(names(result),pattern="^(VAL)-[0-9]+_(VAF)$")]]))/8
result$MEAN_PASS_DP=rowMeans(result[,names(result)[grep(names(result),pattern="^(VAL)-[0-9]+_(DP)$")]],na.rm=T)
result$MEAN_PASS_VAF=rowMeans(result[,names(result)[grep(names(result),pattern="^(VAL)-[0-9]+_(VAF)$")]],na.rm=T)

fwrite(result,"HD753_merge_TP_and_all_testing.tsv",sep="\t",row.names=F,quote=F)

## precision for different cutoff of depth (determine LOD)
DP_precision=data.frame(t(sapply(c(7, 9, 11, 13, 15, 17, 19, 21),function(x){

return(c(
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]<100&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]<100&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]<100&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]<100&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]<100))/length(which(!is.na(result[,x])&result[,x]<100)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=100&result[,x]<200&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>=100&result[,x]<200&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=100&result[,x]<200&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>=100&result[,x]<200&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=100&result[,x]<200))/length(which(!is.na(result[,x])&result[,x]>=100&result[,x]<200)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=200&result[,x]<300&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>=200&result[,x]<300&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=200&result[,x]<300&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>=200&result[,x]<300&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=200&result[,x]<300))/length(which(!is.na(result[,x])&result[,x]>=200&result[,x]<300)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=300&result[,x]<400&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>=300&result[,x]<400&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=300&result[,x]<400&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>=300&result[,x]<400&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=300&result[,x]<400))/length(which(!is.na(result[,x])&result[,x]>=300&result[,x]<400)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=400&result[,x]<500&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>=400&result[,x]<500&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=400&result[,x]<500&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>=400&result[,x]<500&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=400&result[,x]<500))/length(which(!is.na(result[,x])&result[,x]>=400&result[,x]<500)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>500&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>500&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>500&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>500&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>500))/length(which(!is.na(result[,x])&result[,x]>500))
))}
)),stringsAsFactors=F)

colnames(DP_precision)=c("SNV_DP<100","INDEL_DP<100","ALL_DP<100","SNV_100<=DP<200","INDEL_100<=DP<200","ALL_100<=DP<200","SNV_200<=DP<300","INDEL_200<=DP<300","ALL_200<=DP<300","SNV_300<=DP<400","INDEL_300<=DP<400","ALL_300<=DP<400","SNV_400<=DP<500","INDEL_400<=DP<500","ALL_400<=DP<500","SNV_DP>=500","INDEL_DP>=500","ALL_DP>=500")
DP_precision$ID=gsub(colnames(result)[c(7, 9, 11, 13, 15, 17, 19, 21)],pattern="_DP",replacement="")
fwrite(DP_precision,"HD753_DP_precision.tsv",sep="\t",quote=T,row.names=F)


VAF_precision=data.frame(t(sapply(c(8, 10, 12, 14, 16, 18, 20, 22),function(x){

return(c(
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]<0.03&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]<0.03&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]<100&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]<100&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]<100))/length(which(!is.na(result[,x])&result[,x]<100)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=100&result[,x]<200&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>=100&result[,x]<200&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=100&result[,x]<200&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>=100&result[,x]<200&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=100&result[,x]<200))/length(which(!is.na(result[,x])&result[,x]>=100&result[,x]<200)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=200&result[,x]<300&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>=200&result[,x]<300&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=200&result[,x]<300&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>=200&result[,x]<300&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=200&result[,x]<300))/length(which(!is.na(result[,x])&result[,x]>=200&result[,x]<300)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=300&result[,x]<400&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>=300&result[,x]<400&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=300&result[,x]<400&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>=300&result[,x]<400&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=300&result[,x]<400))/length(which(!is.na(result[,x])&result[,x]>=300&result[,x]<400)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=400&result[,x]<500&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>=400&result[,x]<500&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=400&result[,x]<500&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>=400&result[,x]<500&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>=400&result[,x]<500))/length(which(!is.na(result[,x])&result[,x]>=400&result[,x]<500)),

length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>500&result$type=="SNV"))/length(which(!is.na(result[,x])&result[,x]>500&result$type=="SNV")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>500&result$type=="INDEL"))/length(which(!is.na(result[,x])&result[,x]>500&result$type=="INDEL")),
length(which(!is.na(result$TP_DP)&!is.na(result[,x])&result[,x]>500))/length(which(!is.na(result[,x])&result[,x]>500))
))}

)),stringsAsFactors=F)
colnames(VAF_precision)=c("SNV_VAF<10","INDEL_VAF<10","ALL_VAF<10","SNV_10<=VAF<20","INDEL_10<=VAF<20","ALL_10<=VAF<20","SNV_20<=VAF<30","INDEL_20<=VAF<30","ALL_20<=VAF<30","SNV_30<=VAF<40","INDEL_30<=VAF<40","ALL_30<=VAF<40","SNV_40<=VAF<50","INDEL_40<=VAF<50","ALL_40<=VAF<50","SNV_VAF>=50","INDEL_VAF>=50","ALL_VAF>=50")
VAF_precision$ID=gsub(colnames(result)[c(11,16,21,26,31,36,41,46)],pattern="_GT",replacement="")
fwrite(VAF_precision,"NA12878_VAF_precision.tsv",sep="\t",quote=T,row.names=F)

##############################################################################
## recall for different cutoff of depth (determine LOD)
DP_recall=data.frame(t(sapply(c(11,16,21,26,31,36,41,46),function(x){

return(c(
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]<10&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+2]<10&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]<10&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+2]<10&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]<10))/length(which(!is.na(result$TP_GT)&result[,x+2]<10)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=10&result[x+2]<20&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=10&result[x+2]<20&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=10&result[x+2]<20&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=10&result[x+2]<20&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=10&result[x+2]<20))/length(which(!is.na(result$TP_GT)&result[,x+2]>=10&result[x+2]<20)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=20&result[x+2]<30&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=20&result[x+2]<30&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=20&result[x+2]<30&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=20&result[x+2]<30&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=20&result[x+2]<30))/length(which(!is.na(result$TP_GT)&result[,x+2]>=20&result[x+2]<30)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=30&result[x+2]<40&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=30&result[x+2]<40&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=30&result[x+2]<40&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=30&result[x+2]<40&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=30&result[x+2]<40))/length(which(!is.na(result$TP_GT)&result[,x+2]>=30&result[x+2]<40)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=40&result[x+2]<50&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=40&result[x+2]<50&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=40&result[x+2]<50&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=40&result[x+2]<50&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=40&result[x+2]<50))/length(which(!is.na(result$TP_GT)&result[,x+2]>=40&result[x+2]<50)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=50&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=50&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=50&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+2]>=50&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=50))/length(which(!is.na(result$TP_GT)&result[,x+2]>=50))
))}
)),stringsAsFactors=F)

colnames(DP_recall)=c("SNV_DP<10","INDEL_DP<10","ALL_DP<10","SNV_10<=DP<20","INDEL_10<=DP<20","ALL_10<=DP<20","SNV_20<=DP<30","INDEL_20<=DP<30","ALL_20<=DP<30","SNV_30<=DP<40","INDEL_30<=DP<40","ALL_30<=DP<40","SNV_40<=DP<50","INDEL_40<=DP<50","ALL_40<=DP<50","SNV_DP>=50","INDEL_DP>=50","ALL_DP>=50")
DP_recall$ID=gsub(colnames(result)[c(11,16,21,26,31,36,41,46)],pattern="_GT",replacement="")
fwrite(DP_recall,"NA12878_DP_recall.tsv",sep="\t",quote=T,row.names=F)


VAF_recall=data.frame(t(sapply(c(11,16,21,26,31,36,41,46),function(x){

return(c(
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]<0.1&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+4]<0.1&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]<0.1&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+4]<0.1&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]<0.1))/length(which(!is.na(result$TP_GT)&result[,x+4]<0.1)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.1&result[x+4]<0.2&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.1&result[x+4]<0.2&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.1&result[x+4]<0.2&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.1&result[x+4]<0.2&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.1&result[x+4]<0.2))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.1&result[x+4]<0.2)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.2&result[x+4]<0.3&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.2&result[x+4]<0.3&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.2&result[x+4]<0.3&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.2&result[x+4]<0.3&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.2&result[x+4]<0.3))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.2&result[x+4]<0.3)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.3&result[x+4]<0.4&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.3&result[x+4]<0.4&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.3&result[x+4]<0.4&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.3&result[x+4]<0.4&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.3&result[x+4]<0.4))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.3&result[x+4]<0.4)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.4&result[x+4]<0.5&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.4&result[x+4]<0.5&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.4&result[x+4]<0.5&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.4&result[x+4]<0.5&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.4&result[x+4]<0.5))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.4&result[x+4]<0.5)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.5&result$type=="SNV"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.5&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.5&result$type=="INDEL"))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.5&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.5))/length(which(!is.na(result$TP_GT)&result[,x+4]>=0.5))
))}
)),stringsAsFactors=F)
colnames(VAF_recall)=c("SNV_VAF<10","INDEL_VAF<10","ALL_VAF<10","SNV_10<=VAF<20","INDEL_10<=VAF<20","ALL_10<=VAF<20","SNV_20<=VAF<30","INDEL_20<=VAF<30","ALL_20<=VAF<30","SNV_30<=VAF<40","INDEL_30<=VAF<40","ALL_30<=VAF<40","SNV_40<=VAF<50","INDEL_40<=VAF<50","ALL_40<=VAF<50","SNV_VAF>=50","INDEL_VAF>=50","ALL_VAF>=50")
VAF_recall$ID=gsub(colnames(result)[c(11,16,21,26,31,36,41,46)],pattern="_GT",replacement="")
fwrite(VAF_recall,"NA12878_VAF_recall.tsv",sep="\t",quote=T,row.names=F)

