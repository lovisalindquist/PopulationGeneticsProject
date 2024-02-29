
# Project Title: Circos Plot Application based on IBD data of ancient individuals

Author: Lovisa Lindquist
Email: lo3303li-s@student.lu.se
Date: 2024-02-29

# Project Description: 

The input file is a tsv-table, with the following columns:

    Start: start position of segment
    End: end position of segment
    StartM: Morgans at start of segment
    EndM: Morgans at end of segment
    length: End - Start, length of segment
    lengthM: total length in Morgans
    ch: Chromosome
    iid1: ID of individual 1
    iid2: ID of individual 2
    SNPdens: Density of SNPs, calculated by the total length / centiMorgans

    1. Transform data using Python
    2. Generate circos plot in R
    3. Develop application using R Shiny

Input file:
The data used to produce the application can be downloaded at:
https://zenodo.org/records/8417049 

Group reference file downloaded at:
https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data

# File tree
IBD_Application/
    REAMDE.md
    Data/  
        ibd220.ibd.v54.1.pub.tsv.zip
    Src/
        IBD_Transform.py
        IBD_Circos_App.R
    Results/


# Required software:

python v3.11.5
Installation:
conda install python=3.11.5

R v4.2.0
Installation:

Packages:


Installation:

Installation:

###########################################################################################


# 1. Transform dataset using Python

In order to facilitate the downstream analysis in R, the input data file is transformed using the python script IBD_transform.py. 
The script takes one input file, and requires a specified output file. 

It specifically looks for the header names:
    'iid1' - represnting the first individual in the comparison
    'iid2' - representing teh second individual in the comparison
    'ch' - representing the chromosome
    'start' - representing start posiiton of the segment on the given chromosome.
    'end' - representing start end of the segment on the given chromosome.
    'SNP_Dens' - represnting the SNP density, calculated by the Morgans per unit segment length. 

If these column names are not present in the input file, the script terminates with an error message. 
To solve this, either change the column names accordingly, modify the python scirpt to search for the column names specified in the file. 


# 2. Generate circos plot in R

