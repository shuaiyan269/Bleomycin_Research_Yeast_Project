#!/bin/bash
#calculate the copy number, defined as (chromosome median depth / median depth) × 2.

awk -v OFS='\t' '
BEGIN {
    header = "sample";
}

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "median_depth") {
            median_col = i;  
        }
        if (i % 2 == 1 && i >= 11 && $i ~ /_median_depth$/ && $i !~ /MT/) {
            chr_cols[++n] = i;      
            chr_names[n] = $i;     
            sub(/_median_depth/, "", chr_names[n]);  
            header = header OFS chr_names[n] "_ratio";
        }
    }
    if (NR == 1) print header;
}

NR > 1 {
    printf "%s", $1;
    
    for (i = 1; i <= n; i++) {
        ratio = ($(chr_cols[i]) / $median_col)*2;  # (chr_median / total_median)*2
        printf "%s%.3f", OFS, ratio;
    }
    printf "\n";
}' ./../../01_variant_calling/output/01_mapping/all_samples.coverage_summary.txt > ./../output/copy_number.txt

echo "finish!"
