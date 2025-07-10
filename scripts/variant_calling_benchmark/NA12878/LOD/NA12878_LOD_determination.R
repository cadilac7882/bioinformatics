library(data.table)
TP=fread("../TP_hotspot/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz_filter_by_intersectBed_withVAF.vcf.avinput",sep="\t",data.table=F)
TP=TP[,c(1:8,15,18)]
names(TP)=c("chr","start","end","ref","alt","TP_GT","TP_QUAL","TP_DP","TP_FILTER","INFO")
TP$chr=paste0("chr",TP$chr)
TP$TP_VAF=sapply(1:nrow(TP),function(x){
		tmp_VAF=unlist(strsplit(TP$INFO[x],split=":"))
		tmp_VAF=tmp_VAF[length(tmp_VAF)]
		if(grepl(tmp_VAF,pattern=",")){
			if(x>1 & TP$INFO[x]==TP$INFO[x-1]){
				tmp_VAF=unlist(strsplit(tmp_VAF,split=","))[2]
			}else{
				tmp_VAF=unlist(strsplit(tmp_VAF,split=","))[1]
			}
		}
		return(as.numeric(tmp_VAF))
	})
		
result = TP[,-10]

## filtering variants with VAF==0 or DP<30 or FILTER != "PASS"
result = result[-which(result$TP_VAF==0|result$TP_DP<30|result$TP_FILTER!="PASS"),]
max(result$TP_VAF[result$TP_GT=="het"])
min(result$TP_VAF[result$TP_GT=="het"])

## combine all testing samples
sampleList=read.delim("../testing/sampleList.txt",header=F)[,1]
for(i in sampleList){

	tmp=fread(paste0("../testing/",i,"/",i,"_gpu_HF.vcf.gz_filter_by_intersectBed_withVAF.vcf.avinput"),sep="\t",data.table=F)
	tmp=tmp[,c(1:8,15,18)]
	names(tmp)=c("chr","start","end","ref","alt","GT","QUAL","DP","FILTER","INFO")
	
	tmp$VAF=sapply(1:nrow(tmp),function(x){
		tmp_VAF=unlist(strsplit(tmp$INFO[x],split=":"))
		tmp_VAF=tmp_VAF[length(tmp_VAF)]
		if(grepl(tmp_VAF,pattern=",")){
			if(x>1 & tmp$INFO[x]==tmp$INFO[x-1]){
				tmp_VAF=unlist(strsplit(tmp_VAF,split=","))[2]
			}else{
				tmp_VAF=unlist(strsplit(tmp_VAF,split=","))[1]
			}
		}
		return(as.numeric(tmp_VAF))
	})
	tmp=tmp[,-10]
	names(tmp)[-c(1:5)]=paste0(i,c("_GT","_QUAL","_DP","_FILTER","_VAF"))
	
	result=merge(result,tmp,by=c("chr","start","end","ref","alt"),all=T)
}

result$type=NA
result$type[which(result$ref=="-"|result$alt=="-")]="INDEL"
result$type[which(result$ref%in%c("A","T","C","G")&result$alt%in%c("A","T","C","G"))]="SNV"
table(result$type[which(!is.na(result$TP_GT)&result$TP_FILTER=="PASS")])

## remove undetermined variant type
result=result[!is.na(result$type),]

## statistics for variants passing hard-filtering
FILTER_ind=grep(names(result),pattern="^(VAL)-[0-9]+_S[0-9]+_(FILTER)$")
GT_ind=grep(names(result),pattern="^(VAL)-[0-9]+_S[0-9]+_(GT)$")
DP_ind=grep(names(result),pattern="^(VAL)-[0-9]+_S[0-9]+_(DP)$")
VAF_ind=grep(names(result),pattern="^(VAL)-[0-9]+_S[0-9]+_(VAF)$")

result$PASS_Detection_rate=rowSums(result[,FILTER_ind]=="PASS"&result[,GT_ind]==result$TP_GT,na.rm=T)/8
result$MEAN_PASS_DP=sapply(1:nrow(result),function(x){
	ifelse(length(which(result[x,FILTER_ind]=="PASS"&result[x,GT_ind]==result$TP_GT))==0,NA,mean(as.numeric(result[x,DP_ind[which(result[x,FILTER_ind]=="PASS"&result[x,GT_ind]==result$TP_GT)]]),na.rm=T))
})
result$MEAN_PASS_VAF=sapply(1:nrow(result),function(x){
        ifelse(length(which(result[x,FILTER_ind]=="PASS"&result[x,GT_ind]==result$TP_GT))==0,NA,mean(as.numeric(result[x,VAF_ind[which(result[x,FILTER_ind]=="PASS"&result[x,GT_ind]==result$TP_GT)]]),na.rm=T))
})
result$MIN_PASS_DP=sapply(1:nrow(result),function(x){
	ifelse(length(which(result[x,FILTER_ind]=="PASS"&result[x,GT_ind]==result$TP_GT))==0,NA,min(result[x,DP_ind[which(result[x,FILTER_ind]=="PASS"&result[x,GT_ind]==result$TP_GT)]],na.rm=T))
})

result$MIN_PASS_VAF=sapply(1:nrow(result),function(x){
	ifelse(length(which(result[x,FILTER_ind]=="PASS"&result[x,GT_ind]==result$TP_GT))==0,NA,min(result[x,VAF_ind[which(result[x,FILTER_ind]=="PASS"&result[x,GT_ind]==result$TP_GT)]],na.rm=T))
})
## statistics for variants ignoring hard-filtering
result$Detection_rate=rowSums(!is.na(result[,GT_ind]),na.rm=T)/8
result$MEAN_DP=rowMeans(result[,DP_ind],na.rm=T)
result$MEAN_VAF=rowMeans(result[,VAF_ind],na.rm=T)
result$MIN_DP=sapply(1:nrow(result),function(x){
	ifelse(length(which(!is.na(result[x,DP_ind])))==0,NA,min(result[x,DP_ind],na.rm=T))
})
result$MIN_VAF=sapply(1:nrow(result),function(x){
	ifelse(length(which(!is.na(result[x,VAF_ind])))==0,NA,min(result[x,VAF_ind],na.rm=T))
})
#result$MIN_PASS_VAF=sapply(1:nrow(result),function(x){
#	tmp=result[x,c(15,20,25,30,35,40,45,50)][which(result[x,c(14,19,24,29,34,39,44,49)]=="PASS")]
#	ifelse(ncol(tmp)==0,NA,tmp[,which.min(tmp)])
#})
#result$MIN_PASS_DP=sapply(1:nrow(result),function(x){
#	tmp=result[x,c(13,18,23,28,33,38,43,48)][which(result[x,c(14,19,24,29,34,39,44,49)]=="PASS")]
#	ifelse(ncol(tmp)==0,NA,tmp[,which.min(tmp)])
#})

#fwrite(result,"NA12878_merge_TP_and_all_testing.tsv",sep="\t",row.names=F,quote=F)

## precision for different cutoff of depth (determine LOD)
DP_precision=data.frame(t(sapply(c(11,16,21,26,31,36,41,46),function(x){

return(c(
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]<10&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+2]<10&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]<10&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+2]<10&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]<10))/length(which(result[,x+3]=="PASS"&result[,x+2]<10)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=10&result[x+2]<20&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=10&result[x+2]<20&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=10&result[x+2]<20&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=10&result[x+2]<20&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=10&result[x+2]<20))/length(which(result[,x+3]=="PASS"&result[,x+2]>=10&result[x+2]<20)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=20&result[x+2]<30&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=20&result[x+2]<30&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=20&result[x+2]<30&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=20&result[x+2]<30&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=20&result[x+2]<30))/length(which(result[,x+3]=="PASS"&result[,x+2]>=20&result[x+2]<30)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=30&result[x+2]<40&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=30&result[x+2]<40&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=30&result[x+2]<40&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=30&result[x+2]<40&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=30&result[x+2]<40))/length(which(result[,x+3]=="PASS"&result[,x+2]>=30&result[x+2]<40)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=40&result[x+2]<50&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=40&result[x+2]<50&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=40&result[x+2]<50&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=40&result[x+2]<50&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=40&result[x+2]<50))/length(which(result[,x+3]=="PASS"&result[,x+2]>=40&result[x+2]<50)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=50&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=50&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=50&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+2]>=50&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+2]>=50))/length(which(result[,x+3]=="PASS"&result[,x+2]>=50))
))}
)),stringsAsFactors=F)

colnames(DP_precision)=c("SNV_DP<10","INDEL_DP<10","ALL_DP<10","SNV_10<=DP<20","INDEL_10<=DP<20","ALL_10<=DP<20","SNV_20<=DP<30","INDEL_20<=DP<30","ALL_20<=DP<30","SNV_30<=DP<40","INDEL_30<=DP<40","ALL_30<=DP<40","SNV_40<=DP<50","INDEL_40<=DP<50","ALL_40<=DP<50","SNV_DP>=50","INDEL_DP>=50","ALL_DP>=50")
DP_precision$ID=gsub(colnames(result)[c(11,16,21,26,31,36,41,46)],pattern="_GT",replacement="")
fwrite(DP_precision,"NA12878_DP_precision.tsv",sep="\t",quote=T,row.names=F)


VAF_precision=data.frame(t(sapply(c(11,16,21,26,31,36,41,46),function(x){

return(c(
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]<0.1&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+4]<0.1&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]<0.1&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+4]<0.1&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]<0.1))/length(which(result[,x+3]=="PASS"&result[,x+4]<0.1)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.1&result[x+4]<0.2&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.1&result[x+4]<0.2&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.1&result[x+4]<0.2&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.1&result[x+4]<0.2&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.1&result[x+4]<0.2))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.1&result[x+4]<0.2)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.2&result[x+4]<0.3&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.2&result[x+4]<0.3&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.2&result[x+4]<0.3&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.2&result[x+4]<0.3&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.2&result[x+4]<0.3))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.2&result[x+4]<0.3)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.3&result[x+4]<0.4&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.3&result[x+4]<0.4&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.3&result[x+4]<0.4&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.3&result[x+4]<0.4&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.3&result[x+4]<0.4))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.3&result[x+4]<0.4)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.4&result[x+4]<0.5&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.4&result[x+4]<0.5&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.4&result[x+4]<0.5&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.4&result[x+4]<0.5&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.4&result[x+4]<0.5))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.4&result[x+4]<0.5)),

length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.5&result$type=="SNV"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.5&result$type=="SNV")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.5&result$type=="INDEL"))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.5&result$type=="INDEL")),
length(which(result[,x+3]=="PASS"&result[,x]==result$TP_GT&result[,x+4]>=0.5))/length(which(result[,x+3]=="PASS"&result[,x+4]>=0.5))
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

