# smpcr

Download smpcr.py. That's the only file you will need to run the PCR simualation.
###NOTE: 
You need python version 3 to run the script.

### USAGE:
python smpcr.py Nmut Nrounds output_file.txt

Nmut - Number of mutations
Nrounds  - Number of rounds of PCR
outout_file.txt - Mutation Proportions

SMPCR library method by Dong hee Chung and Michael D Toney (UC Davis)
This program calculates the distribution of mutations for each round of SMPCR
Code written by: Krishna Ravikumar (07/27/2016)
For errors/comments/suggestions contact: mdtoney@ucdavis.edu

Output to outfile:    Round   %1mutation   %2mutations  . . .  %Nmutations   Unique_mutations

Unique_combinations DOES NOT include the wild type gene, only mutations. 
