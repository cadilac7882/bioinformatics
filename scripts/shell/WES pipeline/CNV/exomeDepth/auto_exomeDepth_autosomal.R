library(GenomicRanges)
library(ExomeDepth)
setwd("/workdir")
rm(list=ls())

## load regions for evaluation 
data(exons.hg19)
exons.hg19.GRanges = GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
					    IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
					    names = exons.hg19$name)

if(!file.exists("mergedAutosomalBamCount.txt")){

	### select samples for evaluation
	sampleList = read.delim("sampleList.txt",header=F,stringsAsFactor=F)$V1

	my.bam = sapply(sampleList,function(x){
	  paste0(x,"/nv_result/",x,"_gpu.bam")
	})

	### calculate read counts in evaluated regions
	my.counts = getBamCounts(bed.frame = exons.hg19,
							  bam.files = my.bam,
							  include.chr = TRUE)


	myCountsDataFrame = as(my.counts,'data.frame')
	names(myCountsDataFrame)[-c(1:5)] = sampleList
	
	write.table(myCountsDataFrame,"mergedAutosomalBamCount.txt", sep="\t", col.name=T, row.name=F, quote=F)
}else{
	## estimate CNV for specific sample
	myCountsDataFrame = read.delim("mergedAutosomalBamCount.txt",check.names=FALSE)

	
	myCountsDataFrame$space = gsub(myCountsDataFrame$space,pattern="chr",replacement="")
	
	all.samples = colnames(myCountsDataFrame)
	all.samples = all.samples[!all.samples%in%c("space","start","end","width","names" )]

	message(paste0(length(all.samples)," samples are founded."))

	## store the correlation matrix
	corMatrix = matrix(nrow = length(all.samples),ncol=2)
	row.names(corMatrix) = all.samples
	for( i in 1:length(all.samples)){
		message(paste0("Sample",i,":",all.samples[i]))
		
		tmpID = all.samples[i]
		my.test = myCountsDataFrame[,tmpID]
		
		my.ref.samples = all.samples[!all.samples %in% tmpID ]
		
		my.reference.set = as.matrix(myCountsDataFrame[, my.ref.samples])
		
		my.choice = select.reference.set(test.counts = my.test,
						  reference.counts = my.reference.set,
						  bin.length = (myCountsDataFrame$end - myCountsDataFrame$start)/1000,
						  n.bins.reduced = 10000)
									  
		print(my.choice[[1]])
									  
		my.matrix = as.matrix( myCountsDataFrame[, my.choice$reference.choice, drop = FALSE])
		
		my.reference.selected = apply(my.matrix,1,sum)
		
		all.exons = new('ExomeDepth',
					 test = my.test,
					 reference = my.reference.selected,
					 formula = 'cbind(test, reference) ~ 1')
					 
		all.exons = CallCNVs(x = all.exons,
						  transition.probability = 10^-4,
						  chromosome = myCountsDataFrame$space,
						  start = myCountsDataFrame$start,
						  end = myCountsDataFrame$end,
						  name = myCountsDataFrame$names)
						  
		
		all.exons = AnnotateExtra(x = all.exons,
									   reference.annotation = exons.hg19.GRanges,
									   min.overlap = 0.0001,
									   column.name = 'exons.hg19')
		
		all.exons@CNV.calls = all.exons@CNV.calls[order(all.exons@CNV.calls$BF, decreasing=TRUE),]

		output.file = paste0(tmpID,'/', tmpID ,'_autosomal_exome_call.csv')
		
		write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)
		
		## keep correlations
		corMatrix[tmpID,1]=cor(all.exons@test,all.exons@reference)
		corMatrix[tmpID,2]=paste0(my.choice[[1]],collapse=",")
		message("Success!")
	}

	write.table(corMatrix,"mergedAutosomalCorrelationTable.txt",sep="\t",col.names=F,quote=F)
}




