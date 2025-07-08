library(ggplot2)
library(ggExtra)
library(data.table)
setwd("D:/CAPQC/")
#' --------------------------------
#' Define density plot function
#' --------------------------------
depth_gc_original=function(data,output_figure){
  p=ggplot(data=data,aes(x=pct_GC,y=mean_depth_original))+
    geom_point(alpha=0.1,size=0.5)+
    geom_density_2d_filled(aes(alpha = after_stat(level))) +
    scale_y_continuous(limits = c(0,300),name = "Mean depth (x)")+
    scale_x_continuous(limits=c(0,1),breaks = seq(0,1,0.1),labels = seq(0,100,10),
                       name="GC content (%)")+
    theme(axis.line = element_line(),
          legend.position = "none",
          # panel.background = element_rect(fill="white"),
          axis.title = element_text(size=15))
  
  q=ggMarginal(p, type = "densigram", size = 5,
             xparams = list(fill="red",alpha=0.2),
             yparams = list(fill="blue",alpha=0.2))
  pdf(output_figure)
  plot(q)
  dev.off()
}
depth_gc_filtered=function(data,output_figure){
  p=ggplot(data=data,aes(x=pct_GC,y=mean_depth_filtered))+
    geom_point(alpha=0.1,size=0.5)+
    geom_density_2d_filled(aes(alpha = after_stat(level))) +
    scale_y_continuous(limits = c(0,300),name = "Mean depth (x)")+
    scale_x_continuous(limits=c(0,1),breaks = seq(0,1,0.1),labels = seq(0,100,10),
                       name="GC content (%)")+
    theme(axis.line = element_line(),
          legend.position = "none",
          # panel.background = element_rect(fill="white"),
          axis.title = element_text(size=15))
  
  q=ggMarginal(p, type = "densigram", size = 5,
             xparams = list(fill="red",alpha=0.2),
             yparams = list(fill="blue",alpha=0.2))
  
  pdf(output_figure)
  plot(q)
  dev.off()
}

#' ------------------------------
#' evaluate by target regions
#' ------------------------------
samples=read.delim("sampleList.txt",header=F)$V1

gc_content=fread("TruSeq_Exome_TargetedRegions_v1.2_hg19.bed.gc_content.txt",
                 col.names = c("chr","start","end","contig_id","pct_AT","pct_GC","num_A","num_C","num_G","num_T","num_N","num_oth","seq_len"),
                 sep="\t",data.table=F)

for(sample in samples){
  # sample=samples[2]

  cov_origianl_path=paste0(sample,"/bed_coverage_original.txt")
  cov_filtered_path=paste0(sample,"/bed_coverage_filtered.txt")
  
  cov_original=fread(cov_origianl_path,
                     col.names = c("chr","start","end","contig_id","mean_depth_original"),
                     sep="\t",data.table=F)
  cov_original$contig_id=gsub(cov_original$contig_id,pattern="\r",replacement="")
  
  
  # if(all(cov_original[,c(1:3)]==gc_content[,c(1:3)])){
  #   cov_original$contig_id=gc_content$contig_id
  #   write.table(cov_original,paste0(cov_origianl_path,"_name_adjust.txt"),sep="\t",row.names = F,col.names = F,quote = F)
  # }
  
  cov_filtered=fread(cov_filtered_path,
                     col.names = c("chr","start","end","contig_id","mean_depth_filtered"),
                     sep="\t",data.table=F)
  cov_filtered$contig_id=gsub(cov_filtered$contig_id,pattern="\r",replacement="")

  
  # if(all(cov_filtered[,c(1:3)]==gc_content[,c(1:3)])){
  #   cov_filtered$contig_id=gc_content$contig_id
  #   write.table(cov_filtered,paste0(cov_filtered_path,"_name_adjust.txt"),sep="\t",row.names = F,col.names = F,quote = F)
  # }

  merged_result=merge(gc_content,cov_original[,-c(1:3)],by="contig_id")
  merged_result=merge(merged_result,cov_filtered[,-c(1:3)],by="contig_id")

  depth_gc_original(merged_result,paste0(sample,"/",sample,"_target_original_depth_gc.pdf"))
  depth_gc_filtered(merged_result,paste0(sample,"/",sample,"_target_filtered_depth_gc.pdf"))
}


#' --------------------------
#' Evaluate by bins (size~150)
#' -------------------------
gc_content=fread("TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.bed_processed.bed.gc_content.txt",
                 col.names = c("chr","start","end","contig_id","pct_AT","pct_GC","num_A","num_C","num_G","num_T","num_N","num_oth","seq_len"),
                 sep="\t",data.table=F)

for(sample in samples){
  # sample=samples[2]
  cov_origianl_path=paste0(sample,"/bed_coverage_original_bin150.txt")
  cov_filtered_path=paste0(sample,"/bed_coverage_filtered_bin150.txt")
  
  cov_original=fread(cov_origianl_path,
                     col.names = c("chr","start","end","contig_id","mean_depth_original"),
                     sep="\t",data.table=F)
  cov_original$contig_id=gsub(cov_original$contig_id,pattern="\r",replacement="")
  
  
  # if(all(cov_original[,c(1:3)]==gc_content[,c(1:3)])){
  #   cov_original$contig_id=gc_content$contig_id
  #   write.table(cov_original,paste0(cov_origianl_path,"_name_adjust.txt"),sep="\t",row.names = F,col.names = F,quote = F)
  # }
  
  cov_filtered=fread(cov_filtered_path,
                     col.names = c("chr","start","end","contig_id","mean_depth_filtered"),
                     sep="\t",data.table=F)
  cov_filtered$contig_id=gsub(cov_filtered$contig_id,pattern="\r",replacement="")
  
  
  # if(all(cov_filtered[,c(1:3)]==gc_content[,c(1:3)])){
  #   cov_filtered$contig_id=gc_content$contig_id
  #   write.table(cov_filtered,paste0(cov_filtered_path,"_name_adjust.txt"),sep="\t",row.names = F,col.names = F,quote = F)
  # }
  
  merged_result=merge(gc_content,cov_original[,-c(1:3)],by="contig_id")
  merged_result=merge(merged_result,cov_filtered[,-c(1:3)],by="contig_id")
  
  depth_gc_original(merged_result,paste0(sample,"/",sample,"_bin150_original_depth_gc.pdf"))
  depth_gc_filtered(merged_result,paste0(sample,"/",sample,"_bin150_filtered_depth_gc.pdf"))
}




merged_result[grep(merged_result$contig_id,pattern="CEX-chr15-89876276-89877036"),]


median(merged_result$mean_depth_original)



aa=fread("D:/ForUpload/mybam.depth.gz",sep="\t")
head(aa)
round(mean(aa$V3[aa$V1=="chr9" & aa$V2>=101867550 & aa$V2<=101867570]),digits = 1)

2/13
