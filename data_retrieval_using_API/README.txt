SCRIPT OBJECTIVE:

1) For each input protein sequence fasta file, write a file containing the corresponding 
   nucleotide sequences which code for the proteins to enable downstream selection 
   analyses relying on DN/DS metrics.
   
FILES:

get_nuc_OGseqs_v5.1.py: Main script. Relies on NCBI's Entrez server to retrieve accession
                        data. Input sequences are downloaded from NCBI's RefSeq database
                        and headers have been edited to be in the format 
                        >Genus_species_accession_prot_accession. 
                        
test_fa:                Sample input fasta files for script testing.
