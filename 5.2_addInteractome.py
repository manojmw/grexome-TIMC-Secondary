#!/usr/bin/python

# manojmw
# 15 Mar, 2022

import argparse, sys
import pandas as pd
import scipy.stats as stats
import re
import gzip
import logging

###########################################################

# Parses tab-seperated canonical transcripts file
# Required columns are: 'ENSG' and 'GENE' (can be in any order,
# but they MUST exist)
#
# Returns a dictionary:
# - Key -> ENSG
# - Value -> Gene
def ENSG_Gene(inCanonicalFile):

    logging.info("starting to run\n")

    # Dictionary to store ENSG & Gene data
    ENSG_Gene_dict = {}

    # Opening canonical transcript file (gzip or non-gzip)
    try:
        if inCanonicalFile.endswith('.gz'):
            Canonical_File = gzip.open(inCanonicalFile, 'rt')
        else:
            Canonical_File = open(inCanonicalFile)
    except IOError as e:
        sys.exit("E: At Step 5.2_addInteractome - Failed to read the Canonical transcript file: %s" % inCanonicalFile)

    Canonical_header_line = Canonical_File.readline() # Grabbing the header line

    Canonical_header_fields = Canonical_header_line.split('\t')

    # Check the column headers and grab indexes of our columns of interest
    (ENSG_index, Gene_index) = (-1,-1)

    for i in range(len(Canonical_header_fields)):
        if Canonical_header_fields[i] == 'ENSG':
            ENSG_index = i
        elif Canonical_header_fields[i] == 'GENE':
            Gene_index = i

    if not ENSG_index >= 0:
        sys.exit("E: At Step 5.2_addInteractome - Missing required column title 'ENSG' in the file: %s \n" % inCanonicalFile)
    elif not Gene_index >= 0:
        sys.exit("E: At Step 5.2_addInteractome - Missing required column title 'GENE' in the file: %s \n" % inCanonicalFile)
    # else grabbed the required column indexes -> PROCEED

    # Data lines
    for line in Canonical_File:
        line = line.rstrip('\n')
        CanonicalTranscripts_fields = line.split('\t')

        # Key -> ENSG
        # Value -> Gene

        ENSG_Gene_dict[CanonicalTranscripts_fields[ENSG_index]] = CanonicalTranscripts_fields[Gene_index]

    # Closing the file
    Canonical_File.close()

    return ENSG_Gene_dict

###########################################################

# Parses the candidateGenes file in .xlsx format
# Required columns are: 'Gene', 'pathologyID' and
# 'Confidence score' (can be in any order,
# but they MUST exist)
#
# Returns a list with sublist(s)
# One sublist per candidate gene
# Each sublist contains:
# - Gene
# - pathologyID
# - Confidence score)
def CandidateGeneParser(inCandidateFile):

    # Input - List of candidate gene file(s)
    candidate_files = inCandidateFile

    # Initializing an empty data frame to store the data from all files
    meta_data = pd.DataFrame()

    # Data lines
    # Iterating over the list of files and appending to the DataFrame (meta_data)
    for file in candidate_files:
        data = pd.read_excel(file, engine = 'openpyxl')
        meta_data = pd.concat([meta_data, data])

    # Extract Gene, pathologyID and Confidence score and drop rows with missing values(na)
    CandidateGene_data = pd.DataFrame(meta_data, columns=['Gene', 'pathologyID', 'Confidence score']).dropna()

    return CandidateGene_data.values.tolist()

###########################################################

# Parses the dictionary {ENSG_Gene_dict} returned
# by the function: ENSG_Gene
# Also parses the list [CandidateGene_data]
# returned by the function: CandidateGeneParser
#
# Maps Candidate Gene names to ENSG
#
# Returns 2 lists:
# The first list contains sublists
# One sublist per candidate ENSG
# - ENSG identifier of candidate gene
# - pathologyID
# - Confidence score
#
# The second list contains all the pathologies/Phenotypes
def CandidateGene2ENSG(ENSG_Gene_dict, CandidateGene_data):

    # List of pathologies
    pathologies_list = []

    # list to store data associated with each candidate gene
    candidateENSG_out_list = []

    # Data
    for data in CandidateGene_data:
        # data[0] -> Candidate Gene
        # data[1] -> pathologyID
        # data[2] -> Confidence Score

        # ENSG_Gene_dict: Key -> ENSG; Value -> Gene
        for ENSG in ENSG_Gene_dict.keys():
            # Check if candidate gene is present in ENSG_Gene_dict
            if data[0] == ENSG_Gene_dict[ENSG]:
                # Store the ENSG along with other data
                candidateENSG_out = [ENSG,  data[1],  data[2]]
                candidateENSG_out_list.append(candidateENSG_out)

        if not data[1] in pathologies_list:
            pathologies_list.append(data[1])

    return candidateENSG_out_list, pathologies_list

###########################################################

# Parses the candidateENSG_out_list & pathologies_list
#
# Counts the total number of candidate genes
# associated with each pathology
#
# Returns a list with the total count candidate genes (for each pathology)
def CountCandidateGenes(candidateENSG_out_list, pathologies_list):

    # List for counting total candidate genes
    # associated with each pathology
    pathology_CandidateCount = [0] * len(pathologies_list)

    # Data
    for candidateGenedata in candidateENSG_out_list:
        for i in range(len(pathologies_list)):
            if candidateGenedata[1] == pathologies_list[i]:
                pathology_CandidateCount[i] += 1

    return pathology_CandidateCount

###########################################################

# Parses the High-quality Interactome produced by Interactome.py
#
# Required columns:
# First column - ENSG of Protein A
# Second column - ENSG of Protein B
#
# Returns 2 lists:
# First list contains sublists
# One sublist per interaction
# Each sublist contains:
# - ENSG of Protein A
# - ENSG of Protein B
#
# Second list contains all the interacting proteins from the interactome
def Interacting_Proteins(inInteractome):

    # Input - Interactome file
    Interactome_File = open(inInteractome)

    # Dictionary to Interacting proteins
    # from the Interactome
    Interactome_list = []

    # Keeping count of Self Interactions
    SelfInteracting_PPICount = 0

    # List of all interactors from the Interactome
    All_Interactors_list = []

    # Data lines
    for line in Interactome_File:
        line = line.rstrip('\n')

        Interactome_fields = line.split('\t')

        if Interactome_fields[0] != Interactome_fields[1]:
            Interacting_Proteins = [Interactome_fields[0], Interactome_fields[1]]
            Interactome_list.append(Interacting_Proteins)

            # Storing all the interactors in All_Interactors_list
            if not Interactome_fields[0] in All_Interactors_list:
                All_Interactors_list.append(Interactome_fields[0])
            elif not Interactome_fields[1] in All_Interactors_list:
                All_Interactors_list.append(Interactome_fields[1])
        else:
            SelfInteracting_PPICount += 1

    # Closing the file
    Interactome_File.close()

    return Interactome_list, All_Interactors_list

###########################################################

# Parses the UniProt Primary Accession file produced by Uniprot_parser.py
# Required columns are: 'Primary_AC' and 'ENSGs' (can be in any order,
# but they MUST exist)
#
# Also parses the dictionary ENSG_Gene_dict
# returned by the function ENSG_Gene
#
# Maps UniProt Primary Accession to ENSG
# Returns the count of UniProt accessions with unique ENSGs
#
# Count corresponds to total human genes
# This count is later used for calculating
# Benjamini-Hochberg adjusted P-values
def Uniprot_ENSG(inPrimAC, ENSG_Gene_dict):

    # Initializing the dictionary
    Uniprot_ENSG_dict = {}

    UniprotPrimAC_File = open(inPrimAC)

    # Grabbing the header line
    UniprotPrimAC_header = UniprotPrimAC_File.readline()

    UniprotPrimAC_header = UniprotPrimAC_header.rstrip('\n')

    UniprotPrimAC_header_fields = UniprotPrimAC_header.split('\t')

    # Check the column header and grab indexes of our columns of interest
    (UniProt_PrimAC_index, ENSG_index) = (-1, -1)

    for i in range(len(UniprotPrimAC_header_fields)):
        if UniprotPrimAC_header_fields[i] == 'Primary_AC':
            UniProt_PrimAC_index = i
        elif UniprotPrimAC_header_fields[i] == 'ENSGs':
            ENSG_index = i

    if not UniProt_PrimAC_index >= 0:
        sys.exit("E: At Step 5.2_addInteractome - Missing required column title 'Primary_AC' in the file: %s \n" % inPrimAC)
    elif not ENSG_index >= 0:
        sys.exit("E: At Step 5.2_addInteractome - Missing required column title 'ENSG' in the file: %s \n" % inPrimAC)
    # else grabbed the required column indexes -> PROCEED

    # Compiling regular expression

    # Eliminating Mouse ENSGs
    re_ENSMUST = re.compile('^ENSMUSG')

    # Counter for accessions with single canonical human ENSG
    Count_UniqueENSGs = 0

    # Data lines
    for line in UniprotPrimAC_File:
        line = line.rstrip('\n')
        UniprotPrimAC_fields = line.split('\t')

        # ENSG column  - This is a single string containing comma-seperated ENSGs
        # So we split it into a list that can be accessed later
        UniProt_ENSGs = UniprotPrimAC_fields[ENSG_index].split(',')

        # Initializing empty lists
        human_ENSGs = []
        canonical_human_ENSGs = []

        # Eliminating Mouse ENSGs
        for UniProt_ENSG in UniProt_ENSGs:
            if not re_ENSMUST.match(UniProt_ENSG):
                human_ENSGs.append(UniProt_ENSG)
        if not human_ENSGs:
            continue

        # If ENSG is in the canonical transcripts file
        # Append it to canonical_human_ENSGs
        for ENSG in human_ENSGs:
            if ENSG in ENSG_Gene_dict.keys():
                canonical_human_ENSGs.append(ENSG)

        # Keeping the count of protein with single ENSGs
        if len(canonical_human_ENSGs) == 1:
            Count_UniqueENSGs += 1

    return Count_UniqueENSGs

###########################################################

# Parses the Interactome_list & All_Interactors_list returned
# by the function: Interacting_Proteins
#
# Checks the number of interactors for each gene
# Checks the number of known interactors
# using the candidateGene_out_list returned by the function: CandidateGene2ENSG
#
# Returns a Dictionary
# Key -> Gene name; 
# Value -> A List of Interactome data associated with each pathology
# For each pathology, following Interactome data is added:
# - Known Interactors count
# - list of Known Interactors
# - P-value
def Interactors_PValue(Interactome_list, All_Interactors_list, candidateENSG_out_list, pathologies_list, pathology_CandidateCount, Count_UniqueENSGs, ENSG_Gene_dict):

    # Dictionary to store Gene and
    # Interactome data associated with each pathology
    # Key-> Gene; 
    # Value -> A List of Interactome data associated 
    # with each pathology
    Gene_IntAllpatho = {}

    # Checking the number of interactors for each gene
    for ENSG_index in range(len(All_Interactors_list)):

        # List to store Gene name and interactome
        # data associated with each pathology
        # For each pathology, following Interactome 
        # data is added:
        # - Known Interactors count
        # - list of Known Interactors
        # - P-value
        Gene_AllPatho = []

        Gene_AllPatho.append(All_Interactors_list[ENSG_index])

        # List of interactors
        Interactors = []

        for Proteins in Interactome_list:
            # If Protein is the first protein
            if (All_Interactors_list[ENSG_index] == Proteins[0]):
                # Get the interacting protein
                if not Proteins[1] in Interactors:
                    Interactors.append(Proteins[1])
            # If Protein is the Second protein
            elif (All_Interactors_list[ENSG_index] == Proteins[1]):
                if not Proteins[0] in Interactors:
                    # Get the interacting protein
                    Interactors.append(Proteins[0])

        for i in range(len(pathologies_list)):

            # List for known interactor(s)
            Known_Interactors = []

            # Initializing a list to store data for each pathology
            Output_eachPatho = []

            # Checking if the interactor is a known ENSG (candidate ENSG)
            for interactor in Interactors:
                for candidateENSG in candidateENSG_out_list:
                    if interactor in candidateENSG:
                        if candidateENSG[1] == pathologies_list[i]:
                            Known_Interactors.append(interactor)

            # Getting the Gene name for Known Interactors
            for Known_InteractorIndex in range(len(Known_Interactors)):
                Known_Interactors[Known_InteractorIndex] = ENSG_Gene_dict[Known_Interactors[Known_InteractorIndex]]

            if Known_Interactors:
                # Applying Fisher's exact test to calculate p-values
                ComputePvalue_data = [[len(Known_Interactors), len(Interactors)],[pathology_CandidateCount[i], Count_UniqueENSGs]]
                (odds_ratio, p_value) = stats.fisher_exact(ComputePvalue_data)
            # If there are no Known Interactors, 
            # there is no point is computing P-value,
            # So we assign P-value as 1
            else:
                p_value = 1

            if Known_Interactors:
                Output_eachPatho = [len(Known_Interactors), Known_Interactors, p_value]
            else:
                Output_eachPatho = [len(Known_Interactors), '', p_value]

            for data in Output_eachPatho:
                Gene_AllPatho.append(data)

        #  Storing it in the dictionary 
        Gene_IntAllpatho[Gene_AllPatho[0]] = Gene_AllPatho[1:len(Gene_AllPatho)]

    return Gene_IntAllpatho

###########################################################

# Function takes args
# - Canonical transcripts file
# - Candidate Gene(s) file
# - UniProt Primary Accession file
# - Interactome file
#
# Reads on STDIN a TSV file as produced by 5.1_addGTEX.pl
#
# Print to stdout a TSV file with added columns holding
# Interactome data (non-clustering approach)
#
# New columns containing Interactome data are added
# immediately after the 'SYMBOL' column
# This is checked using the Symbol_index
def addInteractome(args):

    # Calling the functions
    CandidateGene_data = CandidateGeneParser(args.inCandidateFile)
    ENSG_Gene_dict = ENSG_Gene(args.inCanonicalFile)
    (Interactome_list, All_Interactors_list) = Interacting_Proteins(args.inInteractome)
    (candidateENSG_out_list, pathologies_list) = CandidateGene2ENSG(ENSG_Gene_dict, CandidateGene_data)
    pathology_CandidateCount = CountCandidateGenes(candidateENSG_out_list, pathologies_list)
    Count_UniqueENSGs = Uniprot_ENSG(args.inPrimAC, ENSG_Gene_dict)
    Gene_IntAllpatho = Interactors_PValue(Interactome_list, All_Interactors_list, candidateENSG_out_list, pathologies_list, pathology_CandidateCount, Count_UniqueENSGs, ENSG_Gene_dict)

    # Take STDIN file as produced by 5.1_addGTEX.pl
    addGTEX_output = sys.stdin

    # Get header line
    addGTEX_out_headerLine = addGTEX_output.readline()

    addGTEX_out_headerLine = addGTEX_out_headerLine.rstrip('\n')

    # Stores headers in a list
    addGTEX_out_headers = addGTEX_out_headerLine.split('\t')

    # Check the column header and grab index of our column of interest
    # Here, we want the grab the indices of columns 'SYMBOL' and 'Gene'
    (Symbol_index, Gene_index) = (-1,-1)

    for i in range(len(addGTEX_out_headers)):
        if addGTEX_out_headers[i] == 'SYMBOL':
            Symbol_index = i
        elif addGTEX_out_headers[i] == 'Gene': 
            Gene_index = i    

    if not Symbol_index >= 0:
        sys.exit("E: At Step 5.2_addInteractome - Missing required column title 'SYMBOL' in STDIN")
    elif not Gene_index >= 0:
        sys.exit("E: At Step 5.2_addInteractome - Missing required column title 'Gene' in STDIN")    
    # else grabbed the required column index -> PROCEED

    # Patho_header_list is a list containing sublists
    # One sublist per pathology (COHORT)
    # Each Sublist contains following header names:
    # - patho_INTERACTORS_COUNT
    # - patho_INTERACTORS
    # - patho_INTERACTORS_PVALUE
    Patho_header_list = [[patho+'_INTERACTORS_COUNT', patho+'_INTERACTORS', patho+'_INTERACTORS_PVALUE'] for patho in pathologies_list]

    # Initializing a list to store all
    # all the elements of a sublist in a
    # single list
    Patho_headers = []

    for headers in Patho_header_list:
        for header in headers:
            Patho_headers.append(header)

    # Adding Interactome data headers to
    # the header of STDIN file as
    # immediately next to the 'SYMBOL' column
    for headerIndex in range(len(Patho_headers)):
        addGTEX_out_headers.insert(headerIndex + (Symbol_index+1), Patho_headers[headerIndex])

    # Printing header
    print('\t'.join(header for header in addGTEX_out_headers))

    # Data lines
    # Add Interactome data and print to STDOUT
    for line in addGTEX_output:
        line = line.rstrip('\n')
        line_fields = line.split('\t')

        # Checking if the Gene present in addGTEX_output exists 
        # as a key in the dictionary: Gene_IntAllpatho
        # If yes, then adding the interactome data associated
        # with this particular gene
        if line_fields[Gene_index] in Gene_IntAllpatho.keys():
            # Inserting Interactome data immediately after the 'SYMBOL' column
            line_fields[Symbol_index+1:Symbol_index+1] = Gene_IntAllpatho[line_fields[Gene_index]]
            print('\t'.join(str(data) for data in line_fields))
        else:
            # If the gene is not present
            # then we leave the fields empty to avoid
            # messing up other columns
            # Since there are 3 types of Interactome data
            # associated with each pathology, we multiply
            # empty string by 3 which is further multiplied by the no. of pathologies
            line_fields[Symbol_index+1:Symbol_index+1] = [''] * 3 * len(pathologies_list)
            print('\t'.join(str(data) for data in line_fields))  

    # Closing the file
    addGTEX_output.close()

    logging.info("5.2_addInteractome.py - ALL DONE, completed successfully!\n")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
------------------------------------------------------------------------------------------------------------------
Program: Reads on STDIN a TSV file as produced by 5.1_addGTEX.pl. Also parses the files provided to the arguments.
         Prints to stdout a TSV file with added columns holding Interactome data (non-clustering approach)
------------------------------------------------------------------------------------------------------------------

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inPrimAC', metavar = "Input File", dest = "inPrimAC", help = 'Uniprot Primary Accession File generated by the UniProt_parser.py', required = True)
    required.add_argument('--inCandidateFile', metavar = "Input File", dest = "inCandidateFile", nargs = '+', help = 'Candidate Genes Input File name(.xlsx)', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file (.gz or non .gz)', required = True)
    required.add_argument('--inInteractome', metavar = "Input File", dest = "inInteractome", help = 'Input File Name (High-quality Human Interactome(.tsv) produced by Build_Interactome.py)', required = True)

    args = file_parser.parse_args()
    addInteractome(args)

if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
