#!/usr/bin/python

# manojmw
# 15 Mar, 2022

import argparse, sys
import logging

###########################################################

# Parses the Interactome results file produced by
# InteractomeResults.py
#
# Returns a List and a dictionary:
# List contains the headers from Interactome results file
# excluding the Gene header
#
# Dictionary contains:
# - Key: Gene (ENSG)
# - Value: a list containing interactome data for each pathology
def InteractomeResults(ininteractomeresults):

    logging.info("starting to run")

    # Initializing the dictionary to store
    # Gene and assoicated interactome data
    Gene_Interactome = {}

    # Initalizing list for headers
    Interactome_headers = []

    # Excute below code only if Interactome file is provided
    if ininteractomeresults:
        IntResults_File = open(ininteractomeresults)

        # Grabbing the header line
        Interactome_headerL = IntResults_File.readline()

        Interactome_headerL = Interactome_headerL.rstrip('\n')

        # Stores headers in a list
        Interactome_headers = Interactome_headerL.split('\t')

        # The first item in the Interactome_headers list
        # is the Gene header, we want to exclude this
        Interactome_headers = Interactome_headers[1:]

        # Data lines
        for line in IntResults_File:
            line = line.rstrip('\n')
            line_fields = line.split('\t')

            # The first item in line_fields list is the Gene (ENSG)
            # The remaining items are Interactome data assoicated with each pathology
            ENSG = line_fields[0]
            Interactome_data = line_fields[1:]

            # Storing the data in dictionary
            # - Key: Gene (ENSG)
            # - Value: a list containing interactome data for each pathology
            if ENSG in Gene_Interactome.keys():
                logging.error("At Step 5.2_addInteractome - Gene found twice in the Interactome results file, impossible. Please fix the file")
                sys.exit()
            else:
                Gene_Interactome[ENSG] = Interactome_data
        
        # Closing the file
        IntResults_File.close()

    return Interactome_headers, Gene_Interactome

###########################################################

# Reads on STDIN a TSV file as produced by 5.1_addGTEX.pl
#
# Print to stdout a similar TSV file with additional columns holding
# Interactome data 
#
# New columns containing Interactome data are added
# immediately after the 'SYMBOL' column
# This is checked using the Symbol_index
def addInteractome(args):

    # Calling the function
    (Interactome_headers, Gene_Interactome) = InteractomeResults(args.ininteractomeresults)

    # Take STDIN file as produced by 5.1_addGTEX.pl
    step5_File = sys.stdin

    # Get header line
    step5_File_headerL = step5_File.readline()

    step5_File_headerL = step5_File_headerL.rstrip('\n')

    # Stores headers in a list
    step5_File_headers = step5_File_headerL.split('\t')

    # Check the column header and grab index of our column of interest
    # Here, we want the grab the indices of columns 'SYMBOL' and 'Gene'
    (Symbol_index, Gene_index) = (-1,-1)

    for i in range(len(step5_File_headers)):
        if step5_File_headers[i] == 'SYMBOL':
            Symbol_index = i
        elif step5_File_headers[i] == 'Gene': 
            Gene_index = i    

    # Sanity check
    if not Symbol_index >= 0:
        logging.error("At Step 5.2_addInteractome - Missing required column title 'SYMBOL' in STDIN")
        sys.exit()
    elif not Gene_index >= 0:
        logging.error("At Step 5.2_addInteractome - Missing required column title 'Gene' in STDIN")    
        sys.exit()
    # else grabbed the required column index -> PROCEED

    # Adding Interactome data headers to
    # the header of STDIN file 
    # immediately next to the 'SYMBOL' column
    for headerIndex in range(len(Interactome_headers)):
        step5_File_headers.insert(headerIndex + (Symbol_index+1), Interactome_headers[headerIndex])

    # Printing header
    print('\t'.join(header for header in step5_File_headers))

    # Data lines
    # Add Interactome data and print to STDOUT
    for line in step5_File:
        line = line.rstrip('\n')
        line_fields = line.split('\t')

        # Checking if the Gene present in addGTEX_output exists 
        # as a key in the dictionary: Gene_IntAllpatho
        # If yes, then adding the interactome data associated
        # with this particular gene
        if line_fields[Gene_index] in Gene_Interactome.keys():
            # Inserting Interactome data immediately after the 'SYMBOL' column
            line_fields[Symbol_index+1:Symbol_index+1] = Gene_Interactome[line_fields[Gene_index]]
            print('\t'.join(str(data) for data in line_fields))
        else:
            pass

    # Closing the file
    step5_File.close()

    logging.info("ALL DONE, completed successfully!")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """

Parse on STDIN a TSV file as produced by 5.1_addGTEX.pl. Print to STDOUT a TSV file (similar to the one produced by 5.1_addGTEX.pl) with additional columns holding Interactome data

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--ininteractomeresults', metavar = "Input File", dest = "ininteractomeresults", help = 'Input File Name (Interactome Results File (.tsv) produced by InteractomeResults.py)')

    args = file_parser.parse_args()
    addInteractome(args)

if __name__ == "__main__":
    # Logging to Standard Error
    logging.basicConfig(format = "%(levelname)s %(asctime)s: %(filename)s - %(message)s", datefmt='%Y-%m-%d %H:%M:%S', stream = sys.stderr, level = logging.DEBUG)
    logging.addLevelName(logging.INFO, 'I' )
    logging.addLevelName(logging.ERROR, 'E')
    logging.addLevelName(logging.WARNING, 'W')
    main()
