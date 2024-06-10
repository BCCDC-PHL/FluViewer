#!/bin/bash

mkdir -p .github/data/assemblies

curl -o .github/data/assemblies/MK583610.1_segment_1_PB2_H3N2.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MK583610.1&db=nucleotide&rettype=fasta"
curl -o .github/data/assemblies/MK583611.1_segment_2_PB1_H3N2.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MK583611.1&db=nucleotide&rettype=fasta"
curl -o .github/data/assemblies/MK583612.1_segment_3_PA_H3N2.fa  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MK583612.1&db=nucleotide&rettype=fasta"
curl -o .github/data/assemblies/MK583613.1_segment_4_HA_H3N2.fa  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MK583613.1&db=nucleotide&rettype=fasta"
curl -o .github/data/assemblies/MK583614.1_segment_5_NP_H3N2.fa  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MK583614.1&db=nucleotide&rettype=fasta"
curl -o .github/data/assemblies/MK583615.1_segment_6_NA_H3N2.fa  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MK583615.1&db=nucleotide&rettype=fasta"
curl -o .github/data/assemblies/MK583616.1_segment_7_M_H3N2.fa   "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MK583616.1&db=nucleotide&rettype=fasta"
curl -o .github/data/assemblies/MK583617.1_segment_8_NS_H3N2.fa  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MK583617.1&db=nucleotide&rettype=fasta"

cat .github/data/assemblies/MK58361*.fa > .github/data/assemblies/tmp.fa

rm .github/data/assemblies/MK58361*.*_segment_*.fa

mv .github/data/assemblies/tmp.fa .github/data/assemblies/MK58361X-H3N2.fa
