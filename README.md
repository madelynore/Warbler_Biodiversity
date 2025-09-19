**Signals of corresponding genetic diversity loss in four warbler species exhibiting regional or rangewide declines**

Here, we leverage DNA sequencing from museum specimens to examine genetic variation between historic and contemporary populations of four species of warblers that vary with respect to degree of population changes over this period: 

Yellow-rumped Warbler (Setophaga coronata), Black-throated Blue Warbler (Setophaga caerulescens),Prairie Warbler (Setophaga discolor), and Blackpoll Warbler (Setophaga striata). 

To explore the genetic impacts of varying population declines, we assembled SNPs across 157 PCR-amplified loci in 341 individuals sampled in two time periods—historic (1789-1955) vs contemporary (2001-2020). 

[Calculating diversity statistics](./Calculate_diversity_stats.Rmd)


Description of the data and file structure

Dataset Overview

These Variant Call Format (VCF) file contains Single Nucleotide Polymorphism (SNP) data from 341 individuals across four species (Black-throated Blue, Yellow-rumped, Blackpoll, and Prairie warblers) and two sampling periods (historic samples from 1789 – 1955: n= 161, range n=25-48 per species; contemporary period from 2000 – 2020: n= 230, range n=42-76 per species. warbler_341_unlinked_variants_highcov.recode.vcf contains 113 SNPs and warbler_341_all_sites_highcov.recode.vcf.gz contains 13,556 SNPs.

Files and variables
File: data/warbler_341_unlinked_variants_highcov.recode.vcf

Description: 

    File Name: warbler_341_unlinked_variants_highcov.recode.vcf
    File Format: VCF (Variant Call Format)
    Date: September 19, 2025
    
File: data/warbler_341_all_sites_highcov.recode.vcf.gz
File Details

    File Name: warbler_341_all_sites_highcov.recode.vcf.gz
    File Format: VCF (Variant Call Format)
    Date: September 18, 2025

Data Description

The VCF files includes the following columns standard to the format:

    #CHROM: Chromosome number
    POS: Position of the SNP on the chromosome
    ID: Identifier of the SNP
    REF: Reference base
    ALT: Alternate base(s)
    QUAL: Quality score of the SNP
    FILTER: Filter status
    INFO: Additional information (e.g., allele frequency, number of samples)
    FORMAT: Data format
    Sample columns: One per individual, containing genotype information

Data was derived from museum specimens housed at the following sources:

    New York State Museum
    Cornell University Museum of Vertebrates
    North Carolina Museum of Natural History
    Yale Peabody Museum
    Delaware Museum of Natural History
    UC Berkeley Museum of Vertebrate Zoology
    Los Angeles County Museum of Natural History
    Chicago Academy of Sciences
    Harvard Museum of Comparative Zoology
    The Ohio State University Museum
    Carnegie Museum

