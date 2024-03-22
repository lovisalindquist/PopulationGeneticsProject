#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Description:
    This script takes two input files, one with all pairwise comparisons and their associated information,
    and one with the origin and age of each ID. The aim is to generate a table with each pairwise comparison, all relevant data 
    and the origins+age of the individuals. Using dictionaries, the origin, age and other information are separately attached to their IDs, 
    and the dictionaries are later used to connect this information and print to teh output file. Note that the input
    IBD_Data file must have the following headers: "iid1", "iid2", "ch" and "snp_d". If not, either change in the input file,
    or adapt the script to the correct headers. 

    
Procedure:
    1. Open files
    2. Identify columns of relevance in IBD_data file
    3. Identify all unique individuals in the file
    4. Create two dictionaries: one for IDs-info, one for IDs-origin+age
        4.1 In both dictionaries, name keys after each individual
    5. Go through file with ID origin and append the origin+age to each ID in the dictionary
    6. Go through IBD_data file and append all necessary info to each comparison to a list
    7. Print results from list to output file

Author: Lovisa Lindquist
Date: 2024-03-16


Usage: python Transform_IBD_File.py IBD_data_file.tsv[IBD data file] AADR_Short.csv[annotation file] Parsed_Ancient_Data.tsv[output]

"""
# import modules
import sys
import os


try:       
    infile = open(sys.argv[1], 'r') #open input file in read mode
    ref_file = open(sys.argv[2], 'r') # open origin file in read mode
    if os.path.exists(sys.argv[3]):
        print("A file with the name of the specified output file already exists. Provide a new name to avoid overwriting important files.")
        sys.exit()
    else: 
        outfile = open(sys.argv[3], 'w') # open output file in write mode
            
    lines = infile.readlines() # read all lines in input file and save in variable 
    firstline = lines[0].split("\t") #extract the first line (header) and save in variable as a list
    ref_lines = ref_file.readlines()
    # iteratate through firstline, identify which columns contain information of interest, save indices to separate variables

    for index, header in enumerate(firstline): 
        if header == "iid1":
            id1 = index # save column index of ID1
        if header == "iid2":
            id2 = index # save column index of ID2
        if header == "ch":
            ch = index
        if header == "lengthM":
            distM = index
        if header == "Start":
            start = index
        if header == "End":
            end = index
                        
    # terminate program if columns could not be extracted
    if id1 is not None and id2 is not None and distM is not None and ch is not None and start is not None and end is not None:
        pass
    else: 
        print("Columns could not be extracted from input file 1. Make sure that the following header names exist: 'iid1', 'iid2', 'SNP_Dens', 'ch', 'Start', 'End'")
        sys.exit()
    
           
    # create list of all individuals present in file
    individuals = [] # create empty list to hold all unique individuals 
    for i in range(1, len(lines)): #iterate through input file
        split_lines = lines[i].split("\t") #create lists of each row
        individuals.append(split_lines[id1]) #append all individuals in ID1 column
        individuals.append(split_lines[id2]) #append all individuals in ID2 column
                
                
    individuals = list(set(individuals)) # only keep one copy of each individual
    dct = {} # create a dictionary that will hold a list for each individual 
    dct_group = {} # create a dictionary that will hold the group belonging
    for i in individuals: # iterate thorugh all individuals
        dct['%s' % i] = [] #create an empty list named as each individual in dictionary
        dct_group['%s' % i] = [] #create an empty list named as each individual in dictionary        
        
    for ind, ind_list in dct_group.items(): #iterate through the dictionary
        for i in range(1, len(ref_lines)): #iterate through all lines in reference file
            if ref_lines[i].split(";")[0] == ind: #if ID matches
                ind_list.append(ref_lines[i].split(";")[1].strip()) #append the group to the ID
                ind_list.append(ref_lines[i].split(";")[2].strip())
        if len(ind_list) < 2: # if not match found, append Unknown origin and age
            ind_list.append("Unknown")
            ind_list.append("Unknown")
               
    for ind, ind_list in dct.items(): # iterate through all IDs
        for i in range(1, len(lines)): # iterate through all lines in input file
            split_lines = lines[i].split("\t") # split files by tab
            if split_lines[id1] == ind: # if the first id matches the key
                ind_list.extend([split_lines[id1], dct_group.get(split_lines[id1])[0], dct_group.get(split_lines[id1])[1], split_lines[id2], dct_group.get(split_lines[id2])[0], dct_group.get(split_lines[id2])[1], split_lines[ch], split_lines[start], split_lines[end], split_lines[distM], "\n"]) #add all necessary information to the key      
            elif split_lines[id2] == ind: # if the second id matches the key
                ind_list.extend([split_lines[id2], dct_group.get(split_lines[id2])[0], dct_group.get(split_lines[id2])[1], split_lines[id1], dct_group.get(split_lines[id1])[0], dct_group.get(split_lines[id1])[1], split_lines[ch], split_lines[start], split_lines[end], split_lines[distM], "\n"])  #add all necessary information to the key    
            else: pass  
                                    
                                    
    # write to output file from dictionary
    outfile.write("ID1\tID1_Group\tID1_Age\tID2\tID2_Group\tID2_Age\tChr\tStart\tEnd\tSegmentLengthM\n") # create header line
    for ind, ind_list in dct.items(): #iterate through dictionary
        i=0 # start at first individual
        while i < len(ind_list): # as long as additional information remains
            outfile.write("\t".join(ind_list[i:i+11])) #write lines to output file
            
            i=i+11 # continue to the next line
           
    print("Output stored in ", sys.argv[3])
 
    infile.close()
    ref_file.close()
    outfile.close()                       
    
except IndexError:
    print("Incorrect file input. Make sure both input and output files are specified and structured according to instructions.")
    sys.exit()
except FileNotFoundError:
    print("Input file does not exist. Make sure path to input file is correct")
    sys.exit()
except NameError:
    print("Incorrect file input. Make sure both input files are specified and structured according to instructions.")
    sys.exit()
except UnicodeDecodeError:
    print("Did you use the correct annotation file? Tip: use AADR_Short.csv")
    sys.exit()

