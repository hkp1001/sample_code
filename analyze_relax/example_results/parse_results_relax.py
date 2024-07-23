#! /usr/bin/env python3

import json
import os
import csv

lines = []

#grab the json output files from the directory
for file in os.scandir("./results_new_version"):
    gene = file.name.split("_")[0]
    if not gene == "RP1":
        char = file.name.split("_")[4][4]
    else:
        gene = file.name.split("_")[0] + "_" + file.name.split("_")[1]
        char = file.name.split("_")[5][4]
    with open (file.path, "r") as in_json:
        result_json = json.load(in_json)
        LRT = result_json["test results"]["LRT"]
        pval = result_json["test results"]["p-value"]
        K = result_json["test results"]["relaxation or intensification parameter"]
    line = [gene, char, LRT, pval, K]
    lines.append(line)

with open("relax_results_03.csv", "w", newline="") as out_csv:
    writer=csv.writer(out_csv, delimiter=",", quotechar=";")
    header=["Gene", "Character scheme", "LRT", "P-value", "K param"]
    writer.writerow(header)
    for line in lines:
        writer.writerow(line)
