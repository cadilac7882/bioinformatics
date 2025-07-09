args <- commandArgs(trailingOnly = TRUE)

show_usage <- function() {
  cat("Usage: select fragments from text file\n")
  cat("arg1: number of fragments\n")
  cat("arg2: read file\n")
  cat("arg3: output folder\n")
  q(status=1)
}

# Check if we have exactly 3 arguments
if (length(args) != 3) {
  show_usage()
}

n_fragments <- as.numeric(args[1])
read_file <- args[2]
output_folder <- args[3]

if(n_fragments<2){
  cat("n_fragments must > 1\n")
  q(status=1)
}

reads_in_region=read.delim(read_file,header=F)
reads_in_region=reads_in_region[,c(1:11)]
names(reads_in_region)=c("read_name","flag","chr","start","mapq","cigar","mate_chr","mate_start","mate_dist","seq","qual")
reads_in_region=reads_in_region[rowSums(is.na(reads_in_region))==0,]
reads_in_region$flag=as.numeric(reads_in_region$flag)

family_table=data.frame(read_name=unique(reads_in_region$read_name))
family_table[,c("family_group","alt_pair","alt_split","max_softclip","sup_aligned")]=data.frame(t(sapply(unique(reads_in_region$read_name),function(x){
  # x="NB552499:205:HVWHYBGXW:1:11110:24464:16471"
  #message(x)
  ## check if any supplementary aligned read exists
  sup_aligned=any(bitwAnd(reads_in_region$flag[reads_in_region$read_name==x],0x800)!=0)

  criteria_met_reads=reads_in_region[reads_in_region$read_name==x & bitwAnd(reads_in_region$flag,0x800)<0x800 & reads_in_region$mapq>=10,]
  fragment_group=NA
  alt_pair=NA
  alt_split=NA
  
  if(nrow(criteria_met_reads)<2){
    return(c(fragment_group,alt_pair,alt_split,NA,sup_aligned))
  }else{
    ## deal with the "chr" of mates aligned on the same chromosome (=)
    criteria_met_reads$mate_chr=ifelse(criteria_met_reads$mate_chr=="=",criteria_met_reads$chr,criteria_met_reads$mate_chr)
    
    ## min and max combinition of chr:position:CIGAR indicate the start and end points
    positions=c(paste0(criteria_met_reads$chr,":",criteria_met_reads$start,":",criteria_met_reads$cigar),
                paste0(criteria_met_reads$mate_chr,":",criteria_met_reads$mate_start,":",criteria_met_reads$cigar))
    
    fragment_group=paste0(min(positions),"-",max(positions))
    
    ## alt_pair = TRUE if chr != mate_chr => two reads are aligned in different chromosome
    alt_pair=TRUE
    if(max(criteria_met_reads$chr,criteria_met_reads$mate_chr)==min(criteria_met_reads$chr,criteria_met_reads$mate_chr)){
      alt_pair=FALSE 
    }
    
    ## alt_split = TRUE if CIGAR contains more than "5S" => reads containing soft-clipped bases
    alt_split=FALSE
    mismatch=sapply(criteria_met_reads$cigar,function(y){
      match=gregexpr(y,pattern="[0-9]+S$|^[0-9]+S")
      max_softclip=max(sapply(1:length(match[[1]]),function(z){
        if(match[[1]][z]<0){
          return(0)
        }else{
          return(as.numeric(substr(y,match[[1]][z],match[[1]][z]+attr(match[[1]],"match.length")[z]-2)))
        }  
      }))
      return(max_softclip)
    })
    if(max(mismatch)>5 && sup_aligned){
      alt_split=TRUE
    }
    # if(any(grepl(criteria_met_reads$cigar,pattern="[0-9]+S$")) | any(grepl(criteria_met_reads$cigar,pattern="^[0-9]+S")))
    #   alt_split=TRUE
    # }
    return(c(fragment_group,alt_pair,alt_split,max(mismatch),sup_aligned))

  }
})))

family_table=family_table[!is.na(family_table$family_group),]
family_table$alt_pair=ifelse(family_table$alt_pair=="TRUE",TRUE,FALSE)
family_table$alt_split=ifelse(family_table$alt_split=="TRUE",TRUE,FALSE)
family_table$max_softclip=as.numeric(family_table$max_softclip)
family_table$sup_aligned=ifelse(family_table$sup_aligned=="TRUE",TRUE,FALSE)
message(paste0(length(unique(family_table$family_group))," groups of fragments are detected."))

family_table=merge(family_table,
      data.frame(family_group=unique(family_table$family_group),family_id=1:length(unique(family_table$family_group))),
      by="family_group")

family_table=family_table[order(family_table$family_id),]
write.table(family_table,
			paste0(output_folder,"/fragment_family_table.tsv"),
			sep="\t",row.names=F,quote=F)

### select fragments
ratio=0.5

split_junction=sample(unique(family_table$family_id[family_table$alt_split&family_table$alt_pair]), ceiling(n_fragments*ratio), replace=FALSE)
pair_junction=sample(unique(family_table$family_id[(!family_table$alt_split)&family_table$alt_pair]), n_fragments-ceiling(n_fragments*ratio), replace=FALSE)

selected_family=c(split_junction,pair_junction)

selected_reads=family_table[family_table$family_id%in%selected_family,]
#selected_reads=selected_reads[order(selected_reads$family_id),]

excluded_reads=family_table[!family_table$family_id%in%selected_family,]
#excluded_reads=excluded_reads[order(excluded_reads$family_id),]

write.table(unique(selected_reads$read_name),paste0(output_folder,"/selected_reads.txt"),sep="\t",row.names=F,col.names = F,quote = F)
write.table(unique(excluded_reads$read_name),paste0(output_folder,"/excluded_reads.txt"),sep="\t",row.names=F,col.names = F,quote = F)
