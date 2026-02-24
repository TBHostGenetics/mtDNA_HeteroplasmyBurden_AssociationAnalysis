# Heteroplasmy burden association analyses 

This repository includes all files and scripts used in the work performed by Croock _et al._ (2025). 
In this study, we aimed to identify the primary factors influencing heteroplasmic burden in southern African cohorts. 

Scripts used in this study are organised into the following directories: 
## 1. Process mitochondrial sequence data
Includes scripts and files used to process VCF files and remove haplogroup-defining variants. Removal of haplogroup-defining variants requires the output from Haplogrep v3.
Following removal of haplogroup-defining vatiants, remaining variants are classified into low and intermediate level heteroplasmic and homoplasmic variants based on variant level or variant allele frequency (VAF) information. These scripts also count the number of variantas in each category and calculate the total heteroplasmic burden per individual. 
This step also generates basic descriptive statistics. 

## 2. Association models 
Includes scripts and files used to run negative binomial regression models. Using these models, we aimed to identify the variables influencing heteroplasmic and homoplasmic variant burden. 
These scripts also include multiple testing correction using FDR-correction methods. 

## 3. Regional analysis
Includes scripts and files used to annotate variants with gene region information (based on rCRS coordinates) and count the number of heteroplasmic variants in each gene region. 
These scripts also include the statistical tests used to determine significant differences in regional heteroplasmic burden between different ancestry categories. 

## Publication
For more information, please see our publication. Please cite our publication when making use of these scripts. 
