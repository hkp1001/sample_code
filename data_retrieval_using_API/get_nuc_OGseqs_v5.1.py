#! /usr/bin/env python3

__author__ = "Hannah K. Pare"
__email__ = "Hannah.Pare@unh.edu"
__version__ = "5.0"

"""
This script reads in OrthoGroup protein sequence fasta files from a directory and uses the accession numbers to query the Entrez
database for the nucleotide sequence which codes for the protein. The input sequences should be from the NCBI RefSeq database.
It takes a path (string) to the directory containing the fasta files and the user's email address for NCBI to contact in case
of problems with the Entrez server. For each input fasta the script will output a fasta file with the corresponding nucleotide
sequences and a json file with corresponding protein and nucleotide accession numbers. If the json file has previously been
created the script will read in the json file in order to avoid needing to re-create the dictionary.
The user should make sure the json file is in the working directory.
"""

import argparse
import Bio.Entrez
import re
import time
import Bio.SeqIO
import os
import json

def main():
    # Argument takes the path to a directory containing sequence fasta files for the OrthoGroups of interest
    parser = argparse.ArgumentParser()
    parser.add_argument("--OGseqs", type=str, required=True, help="Absolute path to the directory containing fasta files for the OrthoGroups of interest")
    parser.add_argument("--email", type=str, required=True, help="Email address to receive any problems from NCBI")
    args=parser.parse_args()

    # make a file for each OrthoGroup containing the nucleotide sequences encoding each protein sequence
    for OG_file in os.scandir(args.OGseqs):
        name = OG_file.name.split("_")[0]
        print("Starting file: {0}".format(OG_file.name))
        get_nuc_seqs(OG_file.path, args.email)
        print()

# function takes the path (a string) to one of the OrthoGroup sequence fasta files and the user's email
# function will write a fasta file containing the corresponding nucleotide sequences
def get_nuc_seqs(prot_file_path, email):
    # store protein accession numbers
    prot_acc_list = []
    try:
        for record in Bio.SeqIO.parse(prot_file_path, "fasta"):
            header_list = record.description.split("_")
            # edited to account for non-formatted headers
            # headers for the sus scrofa RefSeq genome on NCBI are formatted differently
            if not "scrofa" in header_list:
                letters = header_list[5]
                num = header_list[6]
            else:
                letters = header_list[2]
                num = header_list[3]
            accession = "{0}_{1}".format(letters, num)
            prot_acc_list.append(accession)
    except IOError as err:
        print("Error opening file {0}:{1}".format(prot_file_path, err))

    # get the name of the original file for naming purposes. The name is anything before the first "."
    name_tmp = prot_file_path.split("/")[len(prot_file_path.split("/"))-1]
    name = name_tmp.split(".")[0]

    # if the accession dictionary has already been exported as a json file then read it in, if not populate the dictionary and export is as a json file
    if os.path.isfile("{0}_dict.json".format(name)):
        try:
            with open ("{0}_dict.json".format(name), "r") as in_json:
                accession_dict = json.load(in_json)
        except IOError as errA:
            print("Error opening file {0}:{1}".format("{0}_dict.json".format(name), errA))
    else:
        # get a dictionary of corresponding protein and nucleotide accession numbers for the sequences in the file
        accession_dict = get_accessions(prot_acc_list, email)
        try:
            with open ("{0}_dict.json".format(name), "w") as out_json:
                json.dump(accession_dict, out_json)
        except IOError as errB:
            print("Error creating lookup dictionary: {0}".format(errB))

    print("{0}:{1} entries in dictionary".format(name, len(accession_dict)))

    # for each protein sequence in the file, search the NCBI Nucleotide database for the nucleotide sequence in fasta format and write that sequence to a new nucleotide fasta file
    Bio.Entrez.email = email
    try:
        with open("{0}_nucs.fa".format(name), "w") as nuc_file:
            for prot_acc in accession_dict:
                time.sleep(0.5)
                handle = Bio.Entrez.esearch(db="Nucleotide", term=accession_dict[prot_acc])
                results_2 = Bio.Entrez.read(handle)
                handle.close()

                # return nucleotide sequence in fasta format
                if int(results_2["Count"]) > 0:
                    handle = Bio.Entrez.efetch(db="Nucleotide", id=results_2["IdList"][0], rettype = "fasta", retmode = "text")
                    text_data_2 = handle.read()
                    edited_data = text_data_2.rstrip()
                    handle.close()
                    nuc_file.write(edited_data)
                    nuc_file.write("\n")
                else:
                    print("Could not grab nucleotide sequence for {0}".format(prot_acc))
    except IOError as err:
        print("Error writing nuc file: {0}".format(err))
        pass
    except RuntimeError as rterr:
        print("RuntimeError: {0} \n Seq: {1}".format(rterr, prot_acc))
        pass

# function takes a list of refseq protein accession numbers and the user's email
# returns a dictionary with an entry for each protein accession (key) and its corresponding nucleotide accession number (value)
def get_accessions(prot_list, email):
    # regular expression matches the line in the protein entry with the refseq nucleotide accession
    regex = "^DBSOURCE    REFSEQ: accession ([A-Z]{2}_[0-9]+\.[0-9]{1})$"
    accession_dict = {}

    Bio.Entrez.email = email

    # query the NCBI protein database for each protein accession number
    for prot in prot_list:
        # pause between each query in order to not overwhelm the server
        time.sleep(0.5)
        handle = Bio.Entrez.esearch(db="Protein", term=prot)
        results = Bio.Entrez.read(handle)
        handle.close()

        # check to see that the search resulted in at least one record
        if int(results["Count"]) > 0:
            handle = Bio.Entrez.efetch(db="Protein", id=results["IdList"][0], rettype = "gb", retmode = "text")
            text_data = handle.read()
            handle.close()
            lines = text_data.split("\n")
            for line in lines:
                match = re.search(regex, line)
                # identify the line with the nucleotide accession
                if match:
                    accession_dict[prot] = match.group(1)
        # check to see if a match was found for each protein accession. If no match was found, print a warning
        if not prot in accession_dict:
            print("Could not find a match for {0}, possibly a format issue with the NCBI entry".format(prot))

    return accession_dict

if __name__ == "__main__":
    main()
