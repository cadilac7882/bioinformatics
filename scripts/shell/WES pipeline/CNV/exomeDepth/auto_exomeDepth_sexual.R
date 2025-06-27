library(GenomicRanges)
library(ExomeDepth)
setwd("/workdir")
rm(list=ls())

## load regions for evaluation
data(exons.hg19.X)
exons.hg19.X.GRanges = GenomicRanges::GRanges(seqnames = exons.hg19.X$chromosome,
                                            IRanges::IRanges(start=exons.hg19.X$start,end=exons.hg19.X$end),
                                            names = exons.hg19.X$name)

if(!file.exists("mergedSexualBamCount.txt")){

        ### select samples for evaluation
        sampleList = read.delim("sampleListWithGender.txt",header=FALSE,stringsAsFactors=FALSE)
		names(sampleList)=c("sampleID","gender")

        sampleList$bam = sapply(sampleList$sampleID,function(x){
          paste0(x,"/nv_result/",x,"_gpu.bam")
        })

        ### calculate read counts in evaluated regions
        my.counts = getBamCounts(bed.frame = exons.hg19.X,
                                 bam.files = sampleList$bam,
                                 include.chr = TRUE)


        myCountsDataFrame = as(my.counts,'data.frame')
		names(myCountsDataFrame)[-c(1:5)] = sampleList$sampleID
		
		write.table(myCountsDataFrame,"mergedSexualBamCount.txt", sep="\t", col.name=T, row.name=F, quote=F)
}else{
	sampleList = read.delim("sampleListWithGender.txt",header=FALSE,stringsAsFactors=FALSE)
        names(sampleList)=c("sampleID","gender")
	myCountsDataFrame = read.delim("mergedSexualBamCount.txt",sep="\t",check.names=FALSE)
	
	myCountsDataFrame$space = gsub(myCountsDataFrame$space,pattern="chr",replacement="")
	
	all.samples = colnames(myCountsDataFrame)
	all.samples = all.samples[!all.samples%in%c("space","start","end","width","names" )]
	message(paste0(length(all.samples)," samples are founded."))
	
	##input gender information
	## store the correlation matrix
	corMatrix = matrix(nrow = length(all.samples),ncol=2)
	row.names(corMatrix) = all.samples
	for( i in 1:length(all.samples)){
		message(paste0("Sample",i,":",all.samples[i]))
		
		tmpID = all.samples[i]
		my.test = myCountsDataFrame[,tmpID]

		## select samples with the same gender as the testing sample
		my.ref.samples = sampleList$sampleID[sampleList$gender==sampleList$gender[sampleList$sampleID==tmpID]]
		
		my.reference.set = as.matrix(myCountsDataFrame[, my.ref.samples])
		
		my.choice = select.reference.set(test.counts = my.test,
									  reference.counts = my.reference.set,
									  bin.length = (myCountsDataFrame$end - myCountsDataFrame$start)/1000,
									  n.bins.reduced = 10000)
									  
		print(my.choice[[1]])
									  
		my.matrix = as.matrix( myCountsDataFrame[, my.choice$reference.choice, drop = FALSE])
		
		my.reference.selected = apply(my.matrix,1,sum)
		
		print("select test and reference:")
		all.exons = new('ExomeDepth',
					 test = my.test,
					 reference = my.reference.selected,
					 formula = 'cbind(test, reference) ~ 1')
		print("call variant:")			 
		all.exons = CallCNVs(x = all.exons,
						  transition.probability = 10^-4,
						  chromosome = myCountsDataFrame$space,
						  start = myCountsDataFrame$start,
						  end = myCountsDataFrame$end,
						  name = myCountsDataFrame$names)

		if(nrow(all.exons@CNV.calls)==0){
			message("No CNV has been detected!")
		}else{				  
			print("annotation:")	
			all.exons = AnnotateExtra(x = all.exons,
									   reference.annotation = exons.hg19.X.GRanges,
									   min.overlap = 0.0001,
									   column.name = 'exons.hg19.X')
			print("sorting:")	
			all.exons@CNV.calls = all.exons@CNV.calls[order(all.exons@CNV.calls$BF, decreasing=TRUE),]
		}
		output.file = paste0(tmpID,'/', tmpID,'_sexual_exome_call.csv')
		print("write")
		write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)
		
		## keep correlations
		corMatrix[tmpID,1]=cor(all.exons@test,all.exons@reference)
		corMatrix[tmpID,2]=paste0(my.choice[[1]],collapse=",")
		message("Success!")
	}

	write.table(corMatrix,"mergedSexualCorrelationTable.txt",sep="\t",col.names=F,quote=F)
}
