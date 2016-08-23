#!/bin/bash
set -e

for f in *.fasta;
do

#mkdir chimeric_filtered_uchime_ref
usearch -derep_prefix $f -fastaout "./chimeric_filtered_uchime_ref/${f}_amplicons.fasta" -sizeout -minuniquesize 2
usearch -uchime_ref "./chimeric_filtered_uchime_ref/${f}_amplicons.fasta" -db "./chimeric_filtered_uchime_ref/${f}_amplicons.fasta" -self -nonchimeras "./chimeric_filtered_uchime_ref/${f}.ncref.fasta" -strand plus

#mkdir chimeric_filtered_uchime_denovo
usearch -derep_prefix $f -fastaout "./chimeric_filtered_uchime_denovo/${f}_amplicons.fasta" -sizeout -minuniquesize 2
usearch -uchime_denovo "./chimeric_filtered_uchime_denovo/${f}_amplicons.fasta" -nonchimeras "./chimeric_filtered_uchime_denovo/${f}.ncdenovo.fasta"

#mkdir chimeric_filtered_otu_search
usearch -derep_fulllength $f -sizeout -fastaout "./chimeric_filtered_otu_search/${f}_uniques.fa"
usearch -cluster_otus "./chimeric_filtered_otu_search/${f}_uniques.fa" -minsize 2 -otus "./chimeric_filtered_otu_search/${f}_otus.fa" -uparseout "./chimeric_filtered_otu_search/${f}_uparseOut.txt"

done

