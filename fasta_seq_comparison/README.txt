SCRIPT OBJECTIVES:

1) Condense fasta files to include one sequence per species to allow one-to-one species 
   comparisons across genes. 
   
2) Write a fasta file with one sequence per species where the representative sequence is 
   the longest one available for that species.
   
3) Write a fasta file with one sequence per species where the representative sequence is 
   the one with a length closest to the median for all sequences in the file.
   
FILES:

condense_fastas_v3.py: Main script. I used this script to meet the above objectives for my 
                       specific project, however it can be used to condense fasta files 
                       to include just one sequence per any unique entry specified (does 
                       not have to be one per species). Sequences belonging to the same 
                       entry should have the same identifier in the fasta headers 
                       (anything before the first "_"). This script can handle alignment 
                       files and will compare sequences with gaps ("-") removed.

test_fas:              Sample input files for script testing
