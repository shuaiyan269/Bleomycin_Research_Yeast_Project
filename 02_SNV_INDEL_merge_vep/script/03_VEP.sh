#!/bin/bash
##annotate genes by VEP
source /public/shuaiyan/miniconda3/etc/profile.d/conda.sh
conda activate vep_env

id="NA" #BY, WE, WA, NA, SA, MA.
input_vcf="./../output/merge_sample/${id}_snv_indel_merged.vcf"
output="./../output/vep/${id}/all_${id}_snv_indel_vepanno.vcf"
gff="./../input/gff/${id}.gff.gz"
fa="./../input/Reference_Genomes/${id}.tidy.lite.fa"
bed1="./../input/bed1/${id}.repeats_tagged.bed.gz"
###bed1=repeats information
bed2="./../input/bed2/${id}.subtelomere_tagged.bed.gz"
###bed2=subtelomere and telomere
mkdir -p "$(dirname "$output")"

vep -i "$input_vcf" \
 -o "$output" \
 --gff "$gff" \
 --fasta "$fa" \
 --custom file=${bed1},short_name=repeats,format=bed \
 --custom file=${bed2},short_name=subtelomere,format=bed \
 --distance 500 \
 --pick \
 --force_overwrite \
 --vcf 

if [[ $? -eq 0 ]]; then
    echo "VEP annotation completed: $output"
else
    echo "Error: VEP annotation failed!"
    exit 1
fi



