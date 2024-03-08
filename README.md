
# Project Title: Circos Plot Application based on IBD data of ancient individuals

Author: Lovisa Lindquist
Email: lo3303li-s@student.lu.se
Date: 2024-03-07

# Project Description: 
In this project, we develop an application that takes an input file with IBD distances for pairwise comparisons of ancient individuals, and plots relationships between individiduals using circos plots, both at an individual-individual level, population-population level, and individual-population level. 

The input file is a tab separated table, with the following columns:

    Start: start position of segment (bp)
    End: end position of segment (bp)
    StartM: Morgans at start of segment (M)
    EndM: Morgans at end of segment (M)
    length: End - Start, length of segment (bp)
    lengthM: total length in Morgans (bp)
    ch: Chromosome
    iid1: ID of individual 1
    iid2: ID of individual 2
    SNPdens: Density of SNPs, calculated by the total length / centiMorgans
An additional file containing the IDs and their associated geographic origin will be used for the parsing:

The procedure involves the following steps:

    1. Data transformation using Python
        1.1 Read input files and identify column indices of headers
        1.2 Generate dictionaries connecting IDs to information and origin, respectively.
        1.3 Write to output using dictionaries.

    2. Generation of circos plot function in R
        2.1 Import packages and data
        2.2 Generate subdataframes depending on individuals/groups of interest
        2.3 Create matrix for plotting from subdataframe 
        2.4 Generate plot using circlize::chordDiagram

    3. Development of application using R Shiny
        3.1 Define user interface function
            3.1.1 Generate side panel
            3.1.2 Generate main panel
        3.2 Define server function
            3.2.1 Identify and update interactive options
            3.2.2 Plot chordDiagram using function described above
            3.2.3 Generate summary statistics table
            3.2.4 Generate detailed statistics table as drop down from summary table

# Input files
The data used to produce the application can be downloaded at:
https://zenodo.org/records/8417049 
or
https://github.com/lovisalindquist/PopulationGeneticsProject.git

Group reference file downloaded at:
https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
or 
https://github.com/lovisalindquist/PopulationGeneticsProject.git

# File tree

    IBD_Project/
        README.md
        Data/  
            ibd220.ibd.v54.1.pub.tsv
            AADR_Annotation.xlsx
        Scripts/
            Transform_IBD_File.py
            CircosPlotter.R

# Required software:

python v3.11.5
Installation:
conda install python=3.11.5

R v4.2.0
Installation:
conda install -c r=4.2.0
Packages:
- circlize v0.4.16
- shiny v1.8.0
- tibble v3.1.7
- reactable v0.4.4


# 1. Transform dataset using Python

In order to facilitate the downstream analysis in R, and to include data on demographic origin, the input data file is transformed using the python script IBD_transform.py. The script takes two input files, and requires a specified output file. 

It specifically looks for the header names in teh first input file:
    'iid1' - represnting the first individual in the comparison
    'iid2' - representing teh second individual in the comparison
    'ch' - representing the chromosome
    'lengthM' - representing the Morgans for the given segment and comparison
    'Start' - representing start position of segment (bp)
    'End' - representing end position of segment (bp)

If these column names are not present in thefile, the script terminates with an error message. 
To solve this, either change the column names accordingly, modify the python scirpt to search for the column names specified in the file. 

The data is parsed with data on the demographic origin, named AADR_Annotation.xlsx
The file contains cryptic columns that interfere with data extraction by column names or indices, thus,
created a new file with only the two columns of relevance: "Genetic ID" and "Political entity", saved as a .csv file named "AADR_Short.csv"

In IBD_Project/Data/:
`python ../Scripts/Transform_IBD_File.py ibd220.ibd.v54.1.pub.tsv AADR_Short.csv Parsed_Ancient_data.tsv`

The output is stored in Parsed_Ancient_data.tsv

# 2. Generate circos plot in R
In IBD_Project/Scripts/:
`Rscript CircosPlotter.R` 

Copy and Paste the resulting link into a web browser. 
