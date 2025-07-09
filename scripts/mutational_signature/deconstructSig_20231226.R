## Before conducting this program, make sure reference signautres had been pre-processed.
library(RColorBrewer)
library(deconstructSigs)
setwd("D:/Data/mutational signature/")
sampleID="24C00124"
result_Dir=paste0("Result/",sampleID)
dir.create(result_Dir)
message("sample ID: ",sampleID)


#step 1: prepare an input data frame containing sample ID, chr, position, ref, alt
## 1.1 extract DP and VAF from INFO for further filtering
ann_variant=read.delim(paste0("data/Lynch/",sampleID,"_ann.txt"),sep="\t",quote = "")
ann_variant[,c("DP","VAF")]=t(sapply(1:nrow(ann_variant),function(x){
  # head(ann_variant)
  info=unlist(strsplit(ann_variant$Otherinfo13[x],split=":"))
  names(info)=unlist(strsplit(ann_variant$Otherinfo12[x],split=":"))
  return(as.numeric(info[c("DP","VF")]))
}))
message("Before filtering: ",nrow(ann_variant)," variants")


## 1.2 filter by population allele frequency to remove common variant (if you only have tumor sample)
gnomad_filter = ((ann_variant$AF==".")|(as.numeric(ann_variant$AF)<0.01))
X1000G_filter = ((ann_variant$X1000G_ALL==".")|(as.numeric(ann_variant$X1000G_ALL)<0.01))
VAF_filter = as.numeric(ann_variant$VAF)>=0
DP_filter = ann_variant$DP>=0
# Exonic_filter = ann_variant$Func.refGeneWithVer%in%c("exonic","exonic;splicing","splicing")
# Nonsysnonymous_filter = !ann_variant$ExonicFunc.refGene%in%c("synonymous SNV","unknown") 

## 1.3 apply filtering
filtered_variant = ann_variant[gnomad_filter & X1000G_filter & VAF_filter & DP_filter,]
message("After filtering: ",nrow(filtered_variant)," variants")

## 1.4 adjust input format
filtered_variant$Sample=sampleID
maf = filtered_variant[,c("Sample","Chr","Start","Ref","Alt")]
colnames(maf) = c("Sample","chr","pos","ref","alt")
dim(maf)

# step2: Convert to deconstructSigs input
### extract and estimate trinucleotide context (single base subsitution, SBS) by mut.to.sigs.input function
sigs.input.SBS = mut.to.sigs.input(mut.ref = maf,sig.type = "SBS",
                                   sample.id = "Sample",
                                   chr = "chr",
                                   pos = "pos",
                                   ref = "ref",
                                   alt = "alt")

#sum(sigs.input.SBS)==sum(maf$ref%in%c("A","T","C","G")&maf$alt%in%c("A","T","C","G"))

### read reference signatures
signatures.cosmicv34.SBS=readRDS("ref/COSMIC_v3.4_SBS_GRCh37.rds")


### use whichSignatures which takes input context and refence signature to 
### determine weights to assign to each signature in order to best recontruct 
### the mutational profile of the input tumor sample.

## 2.1 Single Base Substitution
sample_SBS = whichSignatures(tumor.ref = sigs.input.SBS, 
                             signatures.ref = signatures.cosmicv34.SBS, 
                             sample.id = sampleID, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'exome2genome')


pdf(paste0(result_Dir,"/",sampleID,"_plotSignatures(SBS).pdf"),width = 10)
plotSignatures(sample_SBS, sub = 'normalized')
dev.off()

tmp_SBS=data.frame(names(sample_SBS$weights),t(sample_SBS$weights))
names(tmp_SBS)=c("signature","weight")
aetiology_map=read.delim("aetiology_map.tsv",sep="\t")
tmp_SBS=merge(tmp_SBS,aetiology_map,by="signature",all.x=T)
tmp_SBS=tmp_SBS[tmp_SBS$weight!=0,]
tmp_SBS$aetiology[is.na(tmp_SBS$aetiology)]="Germline variants contamination"

write.table(tmp_SBS,paste0(result_Dir,"/",sampleID,"_weights(SBS).tsv"),sep="\t",row.names = F,col.names = T)


display.brewer.all()
color_SBS=brewer.pal(sum(sample_SBS$weights!=0), name="Set3")
?brewer.pal
# makePie(sample_SBS, sub = 'example',add.color = color_SBS)

tmp_plot=as.data.frame(t(sample_SBS$weights))
tmp_plot$signature=rownames(tmp_plot)
names(tmp_plot)=c("weight","signature")
tmp_plot=tmp_plot[tmp_plot$weight!=0,]
tmp_plot=merge(tmp_plot,aetiology_map,by="signature",all.x=T)
tmp_plot$aetiology[is.na(tmp_plot$aetiology)]="Possible sequencing artefact"
tmp_plot=rbind(tmp_plot,data.frame(weight=1-sum(tmp_plot$weight),signature="unknown",aetiology=''))


rownames(tmp_plot)=paste0(tmp_plot$signature,"\n",tmp_plot$aetiology,"\n(",round(tmp_plot$weight,4)*100,"%)")
pdf(paste0(result_Dir,"/",sampleID,"_PieSignatures(SBS).pdf"),width = 10)
pie(tmp_plot$weight,labels = rownames(tmp_plot),
    col=c(color_SBS,"grey"),
    main = paste0("SBS signature for case ",sampleID))
dev.off()

### 2.2 Doublet Base Substitution
signatures.cosmicv34.DBS=readRDS("ref/COSMIC_v3.4_DBS_GRCh37.rds")

sigs.input.DBS = mut.to.sigs.input(mut.ref = maf,sig.type = "DBS",
                                   sample.id = "Sample",
                                   chr = "chr",
                                   pos = "pos",
                                   ref = "ref",
                                   alt = "alt")

sample_DBS = whichSignatures(tumor.ref = sigs.input.DBS, 
                             signatures.ref = signatures.cosmicv34.DBS, 
                             sample.id = sampleID, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'exome2genome')

pdf(paste0(result_Dir,"/",sampleID,"_plotSignatures(DBS).pdf"),width = 10)
plotSignatures(sample_DBS, sub = 'default')
dev.off()


tmp_DBS=data.frame(names(sample_DBS$weights),t(sample_DBS$weights))
names(tmp_DBS)=c("signature","weight")
tmp_DBS=merge(tmp_DBS,aetiology_map,by="signature",all.x=T)
tmp_DBS=tmp_DBS[tmp_DBS$weight!=0,]
tmp_DBS$aetiology[is.na(tmp_DBS$aetiology)]="Germline variants contamination"

write.table(tmp_DBS,paste0(result_Dir,"/",sampleID,"_weights(DBS).tsv"),sep="\t",row.names = F,col.names = T)


color_DBS=brewer.pal(sum(sample_DBS$weights!=0), name="Set3")
# makePie(sample_DBS, sub = 'example',add.color = color_DBS)


tmp_plot=as.data.frame(t(sample_DBS$weights))
tmp_plot$signature=rownames(tmp_plot)
names(tmp_plot)=c("weight","signature")
tmp_plot=tmp_plot[tmp_plot$weight!=0,]

tmp_plot=rbind(tmp_plot,data.frame(weight=1-sum(tmp_plot$weight),signature="unknown"))

rownames(tmp_plot)=paste0(tmp_plot$signature," (",round(tmp_plot$weight,4)*100,"%)")
pdf(paste0(result_Dir,"/",sampleID,"_PieSignatures(DBS).pdf"),width = 10)
pie(tmp_plot$weight,labels = rownames(tmp_plot),
    col=c(color_DBS,"grey"),
    main = paste0("DBS signature for case ",sampleID))
dev.off()


# ## 2.3 indel signature
# signatures.cosmicv34.ID=readRDS("ref/COSMIC_v3.4_ID_GRCh37.rds")
#


# sigs.input.ID = mut.to.sigs.input(mut.ref = maf,sig.type = "ID",
#                                     sample.id = "Sample",
#                                     chr = "chr",
#                                     pos = "pos",
#                                     ref = "ref",
#                                     alt = "alt")
#
# ?mut.to.sigs.input()
# sigs.input.ID=read.delim("data/testcc/output/ID/testcc.ID83.all",sep="\t")
# sigs.input.ID2=data.frame(t(sigs.input.ID$test))
# colnames(sigs.input.ID2)=sigs.input.ID$MutationType
# 
# row.names(sigs.input.ID2)=sampleID
# 
# 
# sample_ID = whichSignatures(tumor.ref = sigs.input.ID2,
#                              signatures.ref = signatures.cosmicv34.ID,
#                              sample.id = sampleID,
#                              contexts.needed = TRUE,
#                              tri.counts.method = 'exome2genome')
# 
# 
# #############################

SBS_union=data.frame()
DBS_union=data.frame()
for(i in c("AAS501_DP20","23C00169_DP20","23C00170_DP20")){
  # i="AAS501_DP20"
  tmp_SBS=read.delim(paste0("Result/",i,"/",i,"_weights(SBS).tsv"),sep="\t")
  tmp_DBS=read.delim(paste0("Result/",i,"/",i,"_weights(DBS).tsv"),sep="\t")
  
  if(nrow(SBS_union)==0){
    SBS_union=tmp_SBS
  }else{
    SBS_union=cbind(SBS_union,tmp_SBS)
  }
  
  if(nrow(DBS_union)==0){
    DBS_union=tmp_DBS
  }else{
    DBS_union=cbind(DBS_union,tmp_DBS)
  }
  
}
write.table(SBS_union[which(rowSums(SBS_union)!=0),],"test_SBS.tsv",sep="\t",row.names = T,col.names = T)



write.table(DBS_union[which(rowSums(DBS_union)!=0),],"test_DBS.tsv",sep="\t",row.names = T,col.names = T)
