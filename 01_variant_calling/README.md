# 01_variant_calling
## Server Environment (Bash)

This project uses the Varathon pipeline (https://github.com/yjx1217/Varathon) to perform mapping and variant calling on sequencing data. The workflow involves the following modules: `00.Short_Reads`, `01.Short_Read_Mapping`, `02.Short_Read_SNV_INDEL_Calling`, `03.Short_Read_SV_Calling`, and `04.Short_Read_CNV_Calling`.

Specifically:

First, place the sequencing files in the `00.Short_Reads` folder. The reference genomes for the BY, WE, WA, NA, SA, and MA genetic backgrounds (obtained from https://yjx1217.github.io/Yeast_PacBio_2016/data/) are prepared using the `Varathon.00.Prepare_Reference_Genome.sh` script in the `00.Reference_Genomes` folder.

Second, perform mapping in the `01.Short_Read_Mapping` folder. Edit the `Master_Sample_Table.Batch_yeast_sim_illumina.txt` file to add the corresponding file names and their associated fastq file names. Because mapping requires the corresponding genetic background, one `Master_Sample_Table.Batch_yeast_sim_illumina.txt` file is needed per genetic background; rename "Batch_yeast_sim_illumina" in the file name to "BY", etc., as appropriate. Modify the `Varathon.01.Short_Read_Mapping.sh` script:
```bash
batch_id="Batch_yeast_sim_illumina"   # Replace "Batch_yeast_sim_illumina" with the actual value.
ref_genome="./../00.Reference_Genomes/XX"   # "XX" is the reference genome for the corresponding genetic background.
```

Third, perform SNV and INDEL calling in the `02.Short_Read_SNV_INDEL_Calling` folder. Edit the `Master_Sample_Table.Batch_yeast_sim_illumina.txt` file in the same way, and modify the `Varathon.02.Short_Read_SNV_INDEL_Calling.sh` script:
```bash
batch_id="Batch_yeast_sim_illumina"   # Replace "Batch_yeast_sim_illumina" with the actual value.
ref_genome="./../00.Reference_Genomes/XX"   # "XX" is the reference genome for the corresponding genetic background.
ploidy=2 #The organism in this study is diploid.
snv_indel_caller="gatk4"   # Run both gatk4 and freebayes.
```

Fourth, perform SV calling in the `03.Short_Read_SV_Calling` folder. Edit the `Master_Sample_Table.Batch_yeast_sim_illumina.txt` file, and modify the `Varathon.03.Short_Read_SV_Calling.sh` script:
```bash
batch_id="Batch_yeast_sim_illumina"   # Replace "Batch_yeast_sim_illumina" with the actual value.
ref_genome="./../00.Reference_Genomes/XX"   # "XX" is the reference genome for the corresponding genetic background.
ploidy=2 #The organism in this study is diploid.
sv_caller="manta"   # Run both manta and delly.
```

Fifth, perform CNV calling in the `04.Short_Read_CNV_Calling` folder. Edit the `Master_Sample_Table.Batch_yeast_sim_illumina.txt` file, and modify the `Varathon.04.Short_Read_CNV_Calling.sh` script:
```bash
batch_id="Batch_yeast_sim_illumina"   # Replace "Batch_yeast_sim_illumina" with the actual value.
ref_genome="./../00.Reference_Genomes/XX"   # "XX" is the reference genome for the corresponding genetic background.
ploidy=2 #The organism in this study is diploid.
```
