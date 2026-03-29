#!/bin/bash
#SBATCH --job-name=vcf_filter_assoc
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=long
#SBATCH --output=vcf_filter_assoc.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackgrayner@gmail.com

#FILTER BATCH 1 (high depth)
source activate vcftools
vcftools --gzvcf ../Sw_g3_all.vcf.gz --not-chr scaffold_1 --keep keep_batch1 --minDP 5 --maxDP 100 --minGQ 20 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > Sw_g3_b1_filtered_auto.vcf.gz
vcftools --gzvcf ../Sw_g3_all.vcf.gz --chr scaffold_1 --keep keep_batch1 --minDP 5 --maxDP 70 --minGQ 20 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > Sw_g3_b1_filtered_x.vcf.gz
source activate bcftools
bcftools sort Sw_g3_b1_filtered_auto.vcf.gz -o Sw_g3_b1_filtered_auto_sorted.vcf.gz
bcftools sort Sw_g3_b1_filtered_x.vcf.gz -o Sw_g3_b1_filtered_x_sorted.vcf.gz
bcftools index Sw_g3_b1_filtered_auto_sorted.vcf.gz
bcftools index Sw_g3_b1_filtered_x_sorted.vcf.gz
bcftools concat Sw_g3_b1_filtered_x_sorted.vcf Sw_g3_b1_filtered_auto_sorted.vcf -Oz -o Sw_g3_b1_filtered_concat.vcf.gz


#FILTER BATCH 2 (med. depth)
source activate vcftools
vcftools --gzvcf ../Sw_g3_all.vcf.gz --not-chr scaffold_1 --remove keep_batch1 --minDP 5 --maxDP 60 --minGQ 20 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > Sw_g3_b2_filtered_auto.vcf.gz
vcftools --gzvcf ../Sw_g3_all.vcf.gz --chr scaffold_1 --keep keep_batch1 --minDP 5 --maxDP 42 --minGQ 20 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > Sw_g3_b2_filtered_x.vcf.gz
source activate bcftools
bcftools sort Sw_g3_b2_filtered_auto.vcf.gz -o Sw_g3_b2_filtered_auto_sorted.vcf.gz
bcftools sort Sw_g3_b2_filtered_x.vcf.gz -o Sw_g3_b2_filtered_x_sorted.vcf.gz
bcftools index Sw_g3_b2_filtered_auto_sorted.vcf.gz
bcftools index Sw_g3_b2_filtered_x_sorted.vcf.gz
bcftools concat Sw_g3_b2_filtered_x_sorted.vcf Sw_g3_b2_filtered_auto_sorted.vcf -Oz -o Sw_g3_b2_filtered_concat.vcf.gz

#MERGE and FILTER
source activate bcftools
bcftools merge Sw_g3_b1_filtered_concat.vcf.gz Sw_g3_b2_filtered_concat.vcf.gz -o Sw_merged_b1_b2.vcf.gz
source activate vcftools
vcftools --gzvcf Sw_merged_b1_b2.vcf.gz  --maf 0.1 ---min-alleles 2 --max-alleles 2 --max-missing 0.75 --recode --recode-INFO-all --stdout | bgzip -c > Sw_merged_b1_b2_filtered.vcf.gz

#create bed files
source activate plink2

plink2 --vcf Sw_merged_b1_b2_filtered.vcf.gz  --memory 286368 --maf 0.1 --allow-extra-chr --make-bed --out Sw_merged_b1_b2_filtered
cp ./Sw_g3_all.fam ./Sw_merged_b1_b2_filtered.fam

plink2 --allow-extra-chr --bfile Sw_merged_b1_b2_filtered \
 --chr scaffold_1 --double-id --make-bed --snps-only --pca  \
 --out Sw_merged_b1_b2_filtered_chr1

mkdir gemma_b1_b2
cd ./gemma_b1_b2

source activate GEMMA
gemma -bfile ../Sw_merged_b1_b2_filtered -gk 1 -o Sw_merged_b1_b2_filtered
gemma -bfile ../Sw_merged_b1_b2_filtered -k ./output/Sw_merged_b1_b2_filtered.cXX.txt -lmm 2 -o Sw_merged_b1_b2_filtered
