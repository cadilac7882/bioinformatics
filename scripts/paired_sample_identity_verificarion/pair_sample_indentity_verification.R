#########################################################
#' Date:2025-05-09
#'  Sample identity verification for paired tumor-normal 
#'  sequencing data, using VAF concordance as the primary metric
#########################################################

library(data.table)
library(ggplot2)
library(reshape2)
library(dplyr)

# --- PART 1: Analysis of a single tumor-blood pair ---
setwd("D:/Data/WES_interpretation")

# Load variant data for a blood sample and a tumor sample
blood=fread("22000418_ann.txt",sep="\t")
tumor=fread("25C00126_ann.txt",sep="\t")

# Load a BED file defining targeted regions (e.g., Nextera TSO500 panel regions)
intersect_bed=fread("nextera_TSO500_intersect.bed",sep="\t")

# Extract VAF (Variant Allele Frequency) for the blood sample
blood$blood_VAF=sapply(1:nrow(blood),function(x){
  # Otherinfo13 contains colon-separated values, Otherinfo12 contains corresponding colon-separated keys
  tmp_info=unlist(strsplit(blood$Otherinfo13[x],split=":"))
  names(tmp_info)=unlist(strsplit(blood$Otherinfo12[x],split=":"))
  # Extract the value associated with the 'VAF' key
  as.numeric(tmp_info['VAF'])
})
# Filter blood variants to keep only those that passed quality filters
# blood=blood[blood$Otherinfo10=="PASS",]

# Extract VAF for the tumor sample
tumor$tumor_VAF=
  sapply(1:nrow(tumor),function(x){
    # Similar to blood, but uses 'VF' (Variant Frequency) key, often synonymous with VAF in tumor
    tmp_info=unlist(strsplit(tumor$Otherinfo13[x],split=":"))
    names(tmp_info)=unlist(strsplit(tumor$Otherinfo12[x],split=":"))
    as.numeric(tmp_info['VF'])
  })
# Filter tumor variants to keep only those that passed quality filters
# tumor=tumor[tumor$Otherinfo10=="PASS",]

# Merge tumor and blood data based on variant coordinates (Chr, Start, End, Ref, Alt)
# Keep all tumor variants, and add blood_VAF if a match is found in the blood sample
merged_table=merge(tumor,blood[,c("Chr","Start","End","Ref","Alt","blood_VAF")],all.x=TRUE)

# Determine if each variant in the merged table falls within a targeted region
merged_table$targeted=sapply(1:nrow(merged_table),function(x){
  # Check if the variant's chromosome and start position overlap with any interval in intersect_bed
  return(any(intersect_bed$V1==merged_table$Chr[x] & 
               intersect_bed$V2<=merged_table$Start[x] & 
               intersect_bed$V3>=merged_table$Start[x]))
})

# Convert AF (Allele Frequency from tumor data) to numeric
merged_table$AF=as.numeric(merged_table$AF)
# Replace NA VAFs with 0 (assuming no call means 0 VAF) and ensure numeric type
merged_table$tumor_VAF=ifelse(is.na(merged_table$tumor_VAF),0,as.numeric(merged_table$tumor_VAF))
merged_table$blood_VAF=ifelse(is.na(merged_table$blood_VAF),0,as.numeric(merged_table$blood_VAF))

# Filter for targeted variants
filter=which(merged_table$targeted)

# Create a table with only blood_VAF and tumor_VAF for targeted variants
filterd_table=merged_table[filter,c("blood_VAF","tumor_VAF")]
dim(filterd_table)

# Calculate a metric: sum of absolute VAF differences divided by the number of variants present in blood
# This gives an average VAF difference for variants found in blood.
sum(abs(filterd_table$blood_VAF-filterd_table$tumor_VAF))/sum(filterd_table$blood_VAF!=0)

# Commented out: Alternative denominator - number of variants present in BOTH blood and tumor
#sum(abs(filterd_table$blood_VAF-filterd_table$tumor_VAF))/sum((filterd_table$tumor_VAF!=0)&(filterd_table$blood_VAF!=0))

# Preliminary germline classification based on blood VAF in targeted regions
merged_table$is_Germline=NA
# If targeted and blood VAF is 0, classify as NOT germline (likely somatic or artifact),
merged_table$is_Germline[merged_table$targeted & merged_table$blood_VAF==0]=FALSE
# If targeted and blood VAF is non-zero, classify as germline
merged_table$is_Germline[merged_table$targeted & merged_table$blood_VAF!=0]=TRUE
# Remove variants that were not classified (i.e., not in targeted regions)
merged_table=merged_table[!is.na(merged_table$is_Germline),]

# Select relevant columns for the output table
merged_table=merged_table[,c("Chr","Otherinfo5","Otherinfo7","Otherinfo8","blood_VAF","tumor_VAF","is_Germline")]
# Rename columns for clarity
names(merged_table)=c("Chr","Pos","Ref","Alt","VAF_blood","VAF_tumor","is_Germline")
# Write the table to a CSV file
write.table(merged_table,"25C00126_testing_data_for_germline_predict.csv",sep=",",quote = F,row.names=F)

############
# --- PART 2: Cohort-wide Paired Sample Comparison ---

# Load a pre-collected large CSV file containing variant annotations for a cohort of paired samples
# data.table = F ensures it's read as a standard data.frame, not a data.table
paired_cohort=fread("D:/Data/TSO500/paired_variant_ann_v7_3.csv",sep=",",data.table = F)
dim(paired_cohort) # Display dimensions of the cohort data

# Clean AF column: replace "." with 0 and convert to numeric
paired_cohort$AF=ifelse(paired_cohort$AF==".",0,as.numeric(paired_cohort$AF))

# Filter the cohort data:
# Keep only SNVs (Single Nucleotide Variants)
# Keep only variants that are "overlapped" (presumably meaning they fall within targeted regions)
# Commented out: other potential filters like AF <= 0.1 or exonic variants
paired_cohort=paired_cohort[which(paired_cohort$type=="SNV" &
                                    # paired_cohort$AF<=0.1 &
                                    # paired_cohort$`Func.refGeneWithVer` == "exonic" &
                                    paired_cohort$overlapped ),]
dim(paired_cohort)

# Example of columns structure (commented out):
# sss=paired_cohort[,c("overlapped","AF","Gene.refGeneWithVer","AAChange.refGeneWithVer",paste0("WES_",x,"_GT"),paste0("CGP_",y,"_VAF"))]

# Extract sample IDs from column names (e.g., from "is_Germline_SAMPLEID")
samples=gsub(colnames(paired_cohort)[grep(colnames(paired_cohort),pattern="is_Germline")],pattern="is_Germline_",replacement="")

# Calculate the difference matrix:
# Outer sapply iterates through each sample 'x' treating it as a BLOOD sample
diff_matrix=sapply(samples,function(x){ # x represents the blood sample ID
  # Print progress every 100 samples
  if(grep(x,samples)%%100==0){print(grep(x,samples))}
  
  # Inner sapply iterates through each sample 'y' treating it as a TUMOR sample
  # t() transposes the result so rows are tumors and columns are bloods
  t(sapply(samples,function(y){ # y represents the tumor sample ID
    # Select WES Genotype (GT) for blood sample 'x' and CGP VAF for tumor sample 'y'
    tmp=paired_cohort[,c(paste0("WES_",x,"_GT"),paste0("CGP_",y,"_VAF"))]
    # Convert blood genotype (WES_GT) to numeric VAF: hom=1, het=0.5, ref=0
    tmp[,1]=ifelse(tmp[,1]=="hom",1,ifelse(tmp[,1]=="het",0.5,0))
    # Replace NA tumor VAFs (CGP_VAF) with 0
    tmp[,2]=ifelse(is.na(tmp[,2]),0,tmp[,2])
    
    # CRITICAL FILTER: Consider only variants present in the blood sample (WES_GT != 0)
    tmp=tmp[tmp[,1]!=0,]
    
    # Calculate the mean absolute difference in VAFs between this blood-tumor pair
    # for variants present in the blood sample.
    # If nrow(tmp) is 0 (no variants in blood), this will result in NaN.
    sum(abs(tmp[,1]-tmp[,2]))/nrow(tmp)
  }))
})
# Set column names of the difference matrix (blood samples)
colnames(diff_matrix)=paste0("blood_",colnames(diff_matrix))
# Set row names of the difference matrix (tumor samples)
rownames(diff_matrix)=gsub(colnames(diff_matrix),pattern="blood",replacement="tumor")


# Plotting and analyzing the difference matrix
# Melt the matrix into a long format for ggplot
plotdata=melt(diff_matrix)
# Create columns for tumor and blood IDs by removing prefixes
plotdata$tumor=gsub(plotdata$Var1,pattern="tumor_",replacement="")
plotdata$blood=gsub(plotdata$Var2,pattern="blood_",replacement="")
# Create a boolean column indicating if it's a same-subject pair
plotdata$Same_subject=plotdata$blood==plotdata$tumor

# Order factor levels for plotting based on the VAF difference of same-subject pairs
# This helps in visualizing the diagonal of the heatmap more clearly
tmp=plotdata[plotdata$Same_subject,] # Subset for true pairs
plotdata$tumor=factor(plotdata$tumor,levels = tmp$tumor[order(tmp$value)])
plotdata$blood=factor(plotdata$blood,levels = tmp$blood[order(tmp$value)])

# Create a heatmap of VAF differences
# Low values (red) indicate good concordance, high values (white) indicate poor concordance
# The diagonal should ideally be red.
ggplot(data=plotdata,aes(x=tumor,y=blood,fill=value))+
  geom_tile()+
  scale_fill_gradient(low = "red", high = "white") 

# Calculate variance and summary statistics for VAF differences
# for same-subject pairs vs. different-subject pairs
var(plotdata$value[plotdata$Same_subject],na.rm=TRUE) # Variance for same-subject (true pairs)
summary(plotdata$value[plotdata$Same_subject])      # Summary for same-subject

var(plotdata$value[!plotdata$Same_subject],na.rm=TRUE) # Variance for different-subject (unrelated pairs)
summary(plotdata$value[!plotdata$Same_subject])       # Summary for different-subject

# Identify blood/tumor samples from same-subject pairs that have a high VAF difference (e.g., >= 0.1)
# These could be problematic pairs (sample swaps, poor quality, high tumor heterogeneity/evolution)
tmp$blood[tmp$value>=0.1]
tmp$tumor[tmp$value>=0.1]


# Example of filtering out specific problematic samples from plotdata for further analysis
# These sample IDs were likely identified from previous inspection as outliers or confusing
confusing_samples_list_1 = c("11950458","21558929","14037653","17678442","3175073","20375660","12815172","21769432","15082960","21431894")
plotdata2=plotdata[!plotdata$tumor%in%confusing_samples_list_1 &
                     !plotdata$blood%in%c("11950458","21558929","14037653","17678442"),] # Note: blood list is a subset

# Identify different-subject pairs that have a suspiciously low VAF difference (<0.1)
# These could indicate sample mislabeling or highly related individuals.
plotdata2[!plotdata2$Same_subject & plotdata2$value<0.1,]
head(plotdata2) # Display head of this filtered data

# Example: Using dplyr to filter and arrange data for a specific tumor sample
plotdata2%>%
  filter(tumor=="3175073")%>% # Filter for tumor "3175073"
  arrange(value)%>%           # Arrange by VAF difference (value)
  head()                      # Show the top few matches

# Create a violin plot to compare the distribution of VAF differences
# between same-subject and different-subject pairs (using plotdata2)
a=ggplot(data=plotdata2,aes(x=Same_subject,y=value,fill=Same_subject))+
  geom_violin()
# Add boxplots inside the violins for better summary visualization
# Add a horizontal line at the median VAF difference for same-subject pairs
a+geom_boxplot(width=0.1,fill="white",outlier.size = 0)+
  # geom_jitter(alpha=0.5,width=0.2)+ # Jitter plot (optional)
  geom_hline(yintercept = median(plotdata$value[plotdata$Same_subject],na.rm=TRUE),linetype = "dashed") # Median from original plotdata
# Annotation for the median line (optional)
# annotate("text", x = -Inf, 
#          y = median(plotdata$value[plotdata$Var1==plotdata$Var2],na.rm=TRUE), 
#          label = round(median(plotdata$value[plotdata$Var1==plotdata$Var2],na.rm=TRUE),4),
#          vjust = -1.5,  hjust=-0.5,color = "red")

# Investigate "bad" same-subject pairs (high VAF difference)
same_subject=plotdata[plotdata$Same_subject,] # Get all same-subject pairs from original plotdata
high_diff_blood_samples = same_subject$blood[same_subject$value>0.1 & !is.na(same_subject$value)] # Blood IDs of high-diff pairs

# For each "bad" blood sample (from a high-diff same-subject pair),
# find which tumor sample it actually matches best (lowest VAF difference)
sapply(as.character(high_diff_blood_samples),function(x){ # x is a blood sample ID
  tmp_blood_subset=plotdata[plotdata$blood==x,] # Get all pairs involving this blood sample
  tmp_blood_subset=tmp_blood_subset[order(tmp_blood_subset$value),] # Order by VAF difference
  closest_tumor_sample=as.character(tmp_blood_subset$tumor[1]) # Tumor with the lowest difference
  closest_value=tmp_blood_subset$value[1] # The lowest difference value
  
  # If the closest tumor is indeed its paired sample, print TRUE (it's correctly paired but has high diff)
  # Otherwise, print the blood ID, the ID of the tumor it *actually* matches best, and that VAF difference
  if(closest_tumor_sample==x){ # x is used here because blood and tumor IDs are the same for true pairs
    print(paste0(x, ": Original pair is closest, despite high diff: ", closest_value))
  }else{
    print(paste0("Potential Mismatch for Blood ", x,": Best match is Tumor ",closest_tumor_sample,": Diff=",closest_value))
  }
})

# Example: Inspecting a specific blood sample's matches
aaa=plotdata[plotdata$blood=="17373878",]
aaa[order(aaa$value),] # Show all tumor matches for blood "17373878", ordered by VAF difference

# Investigate "suspiciously good" different-subject pairs (low VAF difference)
diff_subject=plotdata[!plotdata$Same_subject,] # Get all different-subject pairs
# Identify unique sample IDs that are part of a different-subject pair with VAF difference < 0.1
candidates_for_mislabelling=as.character(unique(c(
  diff_subject$tumor[diff_subject$value<0.1 & !is.na(diff_subject$value)],
  diff_subject$blood[diff_subject$value<0.1 & !is.na(diff_subject$value)]
)))

# For each of these candidate samples, check if their *true* pair is actually the best match,
# or if another (incorrect) sample appears to be a better match.
sapply(candidates_for_mislabelling,function(x){ # x is a sample ID (could be blood or tumor context)
  # Check assuming x is a blood sample
  tmp_blood_context=plotdata[plotdata$blood==x,]
  tmp_blood_context=tmp_blood_context[order(tmp_blood_context$value),]
  closest_tumor_to_blood_x=as.character(tmp_blood_context$tumor[1])
  closest_value_blood_context=tmp_blood_context$value[1]
  
  message = paste0("Sample ", x, ": ")
  # If its true pair (where blood_id == tumor_id) is the best match
  if(closest_tumor_to_blood_x==x){
    message = paste0(message, "When considered as BLOOD, its true TUMOR pair is the closest (Diff=", closest_value_blood_context, ").")
  }else{
    message = paste0(message, "When considered as BLOOD, a DIFFERENT Tumor (", closest_tumor_to_blood_x, ") is closer (Diff=", closest_value_blood_context, "). True pair diff: ", plotdata$value[plotdata$blood==x & plotdata$tumor==x], ".")
  }
  
  # Check assuming x is a tumor sample
  tmp_tumor_context=plotdata[plotdata$tumor==x,]
  tmp_tumor_context=tmp_tumor_context[order(tmp_tumor_context$value),]
  closest_blood_to_tumor_x=as.character(tmp_tumor_context$blood[1])
  closest_value_tumor_context=tmp_tumor_context$value[1]
  
  if(closest_blood_to_tumor_x==x){
    message = paste0(message, " When considered as TUMOR, its true BLOOD pair is the closest (Diff=", closest_value_tumor_context, ").")
  }else{
    message = paste0(message, " When considered as TUMOR, a DIFFERENT Blood (", closest_blood_to_tumor_x, ") is closer (Diff=", closest_value_tumor_context, "). True pair diff: ", plotdata$value[plotdata$blood==x & plotdata$tumor==x], ".")
  }
  print(message)
})


# Define a list of "confusing samples" (likely identified from previous runs as problematic)
confusing_samples=c("11950458","21558929","14037653","17678442","3175073","20375660","12815172","21769432","15082960","21431894")
# Create a new plotdata subset excluding these confusing samples from both tumor and blood columns
plotdata2_cleaned=plotdata[!plotdata$tumor%in%confusing_samples &
                             !plotdata$blood%in%confusing_samples,]

# Re-plot violin plot, but this time using the original `plotdata` (not the filtered plotdata2_cleaned)
# This shows the overall distribution including all samples initially processed.
a=ggplot(data=plotdata,aes(x=Same_subject,y=value,fill=Same_subject))+ # Using original `plotdata`
  geom_violin()
a+geom_boxplot(width=0.1,fill="white",outlier.size = 0) # Add boxplots

# Using dplyr with the `plotdata2_cleaned` (after removing confusing samples):
# Filter for pairs that are NOT the same subject AND have a VAF difference < 0.1
# This shows potential remaining mislabels/related pairs after initial filtering of known problematic ones.
plotdata2_cleaned%>%
  filter(!Same_subject)%>%
  filter(value<0.1)

# --- END OF FILE ---