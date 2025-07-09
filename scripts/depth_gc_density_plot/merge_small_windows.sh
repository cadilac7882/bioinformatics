#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_bed_file> <output_bed_file>"
  exit 1
fi

input_bed="$1"
output_bed="$2"

echo "Input BED file: $input_bed"
echo "Output BED file: $output_bed"

awk -v input="$input_bed" -v output="$output_bed" '
BEGIN{OFS="\t"; prev_chr=""; prev_start=-1; prev_end=-1; prev_line=""; prev_name=""}
{
  chr = $1;
  start = $2;
  end = $3;
  name = $4;
  line = $0;
  current_size = end - start;

  if (current_size >= 100) {
    if (prev_chr != "") {
      print prev_line > output;
    }
    prev_chr = chr;
    prev_start = start;
    prev_end = end;
    prev_line = line;
    prev_name = name;
  } else {
    if (prev_chr == chr && prev_start != -1) {
      merged_start = prev_start;
      merged_end = end;
      merged_line = prev_chr OFS merged_start OFS merged_end OFS prev_name;
      for (i=5; i<=NF; i++) {
        merged_line = merged_line OFS $i;
      }
      prev_start = merged_start;
      prev_end = merged_end;
      prev_line = merged_line;
      # prev_name remains the same as the former window
    } else {
      # If it is the first small window of a chromosome, buffer it
      prev_chr = chr;
      prev_start = start;
      prev_end = end;
      prev_line = line;
      prev_name = name;
    }
  }
}

END {
  if (prev_chr != "") {
    print prev_line > output;
  }
}' "$input_bed"

echo "Merged small windows saved to: $output_bed"