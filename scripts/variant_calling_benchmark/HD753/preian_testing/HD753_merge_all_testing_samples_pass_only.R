library(data.table)
TP=fread("../TP_hotspot/HD753-normed-realigned-DP50-TSO500_AGILENT_intersect_notindifficult.vcf.avinput",sep="\t",data.table=F)
names(TP)[1:5]=c("chr","start","end","ref","alt")
TP$WES_DP=sapply(1:nrow(TP),function(x){
	as.numeric(gsub(unlist(strsplit(TP$V13[x],split=";"))[1],pattern="DP=",replacement=""))
})
result=TP[,c("chr","start","end","ref","alt","WES_DP")]

sampleList=read.delim("sampleList",header=F)[,1]
for(i in sampleList){
        tmp=fread(paste0(i,"/main.vcf_exclude_CNV-TSO500_AGILENT_intersect_notindifficult.vcf.avinput"),sep="\t",data.table=F)
        tmp=tmp[,c(1:5,14,15)]
        names(tmp)=c("chr","start","end","ref","alt","HEADER","INFO")

        tmp[,c("DP","VAF")]=t(sapply(1:nrow(tmp),function(x){
			tmp_header=unlist(strsplit(tmp$HEADER[x],split=":"))
			tmp_info=unlist(strsplit(tmp$INFO[x],split=":"))
			return(c(tmp_info[grep(tmp_header,pattern="DP")],tmp_info[grep(tmp_header,pattern="VF")]))
			
		}))
        tmp$DP=as.numeric(tmp$DP)
        tmp$VAF=as.numeric(tmp$VAF)
        tmp=tmp[,-c(6,7)]
        names(tmp)[-c(1:5)]=paste0(i,c("_DP","_VAF"))

        result=merge(result,tmp,by=c("chr","start","end","ref","alt"),all=T)
}

result$Variant_type=NA
result$Variant_type[which(result$ref=="-"|result$alt=="-")]="INDEL"
result$Variant_type[which(result$ref%in%c("A","T","C","G")&result$alt%in%c("A","T","C","G"))]="SNV"
#result$Variant_type[is.na(result$Variant_type)]="MNV"

## remove undetermined variant type
result=result[!is.na(result$Variant_type),]

#coa=read.delim("../TP_hotspot/HD753_CoA_list.txt",sep="\t")
#coa$Verified="v"
#result=merge(result,coa,by=c("chr","start","end","ref","alt"),all=T)


## statistics for variants passing hard-filtering
result$Detection_rate=rowSums(!is.na(result[,names(result)[grep(names(result),pattern="^(VAL-)[0-9]+_(VAF)$")]]))/8
result$MEAN_DP=rowMeans(result[,names(result)[grep(names(result),pattern="^(VAL-)[0-9]+_(DP)$")]],na.rm=T)
result$MEAN_VAF=rowMeans(result[,names(result)[grep(names(result),pattern="^(VAL-)[0-9]+_(VAF)$")]],na.rm=T)


result[,c("MIN_DP","MIN_VAF")]=data.frame(t(sapply(1:nrow(result),function(x){
	if(all(is.na(result[x,names(result)[grep(names(result),pattern="^(VAL-)[0-9]+_(DP)$")]]))){
		return(c(NA,NA))
	}else{
		return(c(min(result[x,names(result)[grep(names(result),pattern="^(VAL-)[0-9]+_(DP)$")]],na.rm=T),
			min(result[x,names(result)[grep(names(result),pattern="^(VAL-)[0-9]+_(VAF)$")]],na.rm=T)))
	}
})))

#result=result[order(result$Verified),]

fwrite(result,"HD753-merged-perian-testing-cases-pass-only-TSO500_AGILENT_intersect_notindifficult.tsv",sep="\t",row.names=F,quote=F)
