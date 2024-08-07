#! /usr/bin/env python3

__author__ = "Hannah K. Pare"
__email__ = "hannah.pare11@gmail.com"
__version__ = "3.0"

import argparse
import string
import Bio.SeqIO
import os

def main():
    """For each input file: group the seqs for each unique entry and write two output files with one representative 
    seq per group (seq with length closest to median and longest seq per entry)"""

    # script takes the path (string) to a directory containing the fasta files to be condensed
    parser = argparse.ArgumentParser(description= """This script condenses fasta files to include just one 
    sequence per unique entry in the file (e.g. one sequence per species, one sequence per sample, etc.). Each 
    unique entry should have the same identifier in the sequence headers (anything before the first "_"). It 
    writes two output files for each input file. One output file will contain the longest sequences for each 
    unique entry and the other output file will contain the sequences closest to the median sequence length among
    all sequences in the file for each unique entry.""")
    parser.add_argument("-d", type=str, required=True, help="Directory containing fasta files to be condensed")
    args=parser.parse_args()

    for file in os.scandir(args.d):
        
        # store all sequences in the file
        all_seqs = []
        
        # group the sequences for each unique entry in the file 
        # an entry is defined as everything in the header before the first "_"
        # keys = entries, values = seq record(s) for that entry
        grouped_dict = {}

        try:
            for record in Bio.SeqIO.parse(file.path, "fasta"):
                all_seqs.append(record.seq)
                entry = record.id.split("_")[0]
                grouped_dict.setdefault(entry, [])
                grouped_dict[entry].append(record)
        except IOError as err1:
            print(f"Error opening the fasta file {file.path}: {err1}")
        
        # get the median sequence length for all seqs in the file
        median_length = get_median(all_seqs)

        # store one representative seq for each entry
        # seq with the length closest to the median
        median_seqs = []
        
        # longest seq
        longest_seqs = []

        for entry in grouped_dict:
            
            # if the entry only has one representative seq choose that seq
            if len(grouped_dict[entry]) == 1:
                median_seqs.append(grouped_dict[entry][0])
                longest_seqs.append(grouped_dict[entry][0])
            else:
                # choose the seq with a length closest to the median
                keep_index_1 = get_median_seq(grouped_dict[entry], median_length)
                median_seqs.append(grouped_dict[entry][keep_index_1])
                
                # choose the longest seq
                keep_index_2 = longest_seq(grouped_dict[entry])
                longest_seqs.append(grouped_dict[entry][keep_index_2])

        # write a new file with one seq per entry (length is closest to the median seq length)
        tmp_name = file.name.split(".fa")[0]
        new_name = tmp_name + "_condensed_median.fa"
        try:
            Bio.SeqIO.write(median_seqs, new_name, "fasta")
        except IOError as err2:
            print(f"Error writing the output file {new_name}: {err2}")

        # write a new file with one seq per entry (longest seq)
        new_name_2 = tmp_name + "_condensed_long.fa"
        try:
            Bio.SeqIO.write(longest_seqs, new_name_2, "fasta")
        except IOError as err3:
            print(f"Error writing the output file {new_name_2}: {err3}")

def get_median(seqs):
    """Get the median sequence length for all seqs in the file. Takes a list of sequences as input"""
    
    # store each individual sequence length in a list
    ed_seqs_len = []

    for seq in seqs:
        
        # remove any non-sequence characters
        seq = remove_sp_chars(seq)
        ed_seqs_len.append(len(seq))
    ed_seqs_len = sorted(ed_seqs_len)
    num_seqs = len(ed_seqs_len)

    # if there is an even number of seqs, the median is the average of the middle pair
    if num_seqs % 2 == 0:
        index1 = int((num_seqs/2) - 1)
        index2 = int(num_seqs/2)
        median = sum((ed_seqs_len[index1], ed_seqs_len[index2]))/2
    else:
        index = int(num_seqs/2 - 0.5)
        median = ed_seqs_len[index]
    
    return median

def remove_sp_chars(seq):
    """Remove any non-alphabet characters from a sequence. This is particularly helpful for alignment files, in
    order to compare seqs without the gaps ("-"). Takes a sequence as input and returns the sequence without 
    any non-alphabet chars."""

    seq = str(seq).upper()
    letters = set(string.ascii_uppercase)
    non_alphabet_chars = set(seq).difference(letters)
    
    # if the sequence contains non-alphabet characters remove them 
    if len(non_alphabet_chars) > 0:
        for chr in non_alphabet_chars:
            seq = seq.replace(chr, "")

    return seq

def get_median_seq(seq_records, median):
    """Compare seq length to the median and identify the seq with the closest length. Takes a list of sequence 
    records as input and the median value and returns the index number of the chosen sequence"""
    
    # for all seqs, store the index number of the sequence in the list (key) and the difference between the seq length and the median length (value)
    diff_dict = {}
    for num in range(0, len(seq_records)):
        cur_seq = seq_records[num].seq
        
        # remove any non-sequence characters
        cur_seq = remove_sp_chars(cur_seq)
        
        # calculate the difference in length between the current seq and the median
        diff = abs(median - len(cur_seq))
        diff_dict[num] = diff
    
    # identify the seq with a length closest to the median
    # if the difference is the same for all seqs, choose the first seq
    minimum = min(diff_dict.values())
    for index in diff_dict:
        if diff_dict[index] == minimum:
            keep_index = index
            break
    
    return keep_index

def longest_seq(seq_records):
    """Identify the longest seq. Takes a list of sequence records as input and returns the index number of the 
    longest sequence"""
    
    # for all seqs, store the index number of the sequence in the list (key) and the sequence length (value)
    seq_dict = {}
    for num in range(0, len(seq_records)):
        cur_seq = seq_records[num].seq
        
        # remove any non-sequence characters
        cur_seq = remove_sp_chars(cur_seq)
        seq_dict[num] = len(cur_seq)
    
    # get the index number of the longest sequence
    # if there are multiple "longest" seqs choose the first one
    maximum = max(seq_dict.values())
    for index in seq_dict:
        if seq_dict[index] == maximum:
            keep_index = index
            break
    
    return keep_index

if __name__ == "__main__":
    main()