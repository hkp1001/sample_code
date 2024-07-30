#! /usr/bin/env python3

import argparse
import json
import os
import csv

"""
This script parses the raw output json files from the RELAX analysis and condenses them into a single csv file. It uses information from 
the file names according to my specific format. 
"""

parser = argparse.ArgumentParser()
parser.add_argument("-i", type = str, required = True, help = "Path to directory containing raw results json files")
args=parser.parse_args()

# store output file entries
lines = []

# grab the json output files from the directory
for file in os.scandir(args.i):
    # grab info from file name
    # gene name defined as everything before the first "_" in the file name
    gene = file.name.split("_")[0]
    # fix for file name issue with RP1 genes
    if not gene == "RP1":
        char = file.name.split("_")[4][4]
    else:
        gene = file.name.split("_")[0] + "_" + file.name.split("_")[1]
        char = file.name.split("_")[5][4]
    # read in and grab info from json result file
    try:
        with open (file.path, "r") as in_json:
            result_json = json.load(in_json)
            LRT = result_json["test results"]["LRT"]
            pval = result_json["test results"]["p-value"]
            K = result_json["test results"]["relaxation or intensification parameter"]
    except IOError as err:
        print("Error reading json file {0}: {1}".format(file.path, err))
    line = [gene, char, LRT, pval, K]
    lines.append(line)

# write output csv file
# includes columns for the gene name, the character scheme number, LRT value, p-value, and K parameter
with open("relax_results.csv", "w", newline="") as out_csv:
    writer=csv.writer(out_csv, delimiter=",", quotechar=";")
    header=["Gene", "Character scheme", "LRT", "P-value", "K param"]
    writer.writerow(header)
    for line in lines:
        writer.writerow(line)
