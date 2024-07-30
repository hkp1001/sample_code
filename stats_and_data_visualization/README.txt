OBJECTIVES FOR THE SCRIPTS IN THIS DIRECTORY:

1) Identify genes where the strength of selection has been intensified or relaxed along 
   branches leading to non-nocturnal species using the RELAX analysis (Wertheim et al., 
   2015; https://www.hyphy.org/methods/selection-methods/). 
 
2) Visualize the results including functional classifications for each gene.

3) Test whether the genes experiencing shifts in the strength of natural selection along 
   non-nocturnal lineages are involved in similar functional processes.

FILES:

relax.sh:                       run RELAX for each gene to be tested

input_files:                    sample input files for relax.sh

example_results:                sample result files generated as output of the RELAX 
                                analysis

parse_results_relax.py:         parse the raw output files and condense results into a 
                                .csv file (relax_results.csv). See .txt file for 
                                information on how to interpret the results in the .csv 
                                file

analyze_relax_results.Rmd:      Main script used for data manipulation, visualization and 
                                statistical testing. Written to be run in R.

*relaxed_intensified_genes.csv: Input file for analyze_relax_results.Rmd. Manipulated 
                                results from the RELAX analysis to only include 
                                significant genes and separated by character scheme and 
                                selection intensity (relaxed vs intensified). A value of 1
                                in this file implies significant results.

*gene_cats_simple.csv:          Second input file for analyze_relax_results.Rmd. 
                                Functional classification info for each gene.


 Wertheim, JO et al. "RELAX: detecting relaxed selection in a phylogenetic framework." 
 Mol. Biol. Evol. 32, 820-832 (2015).
