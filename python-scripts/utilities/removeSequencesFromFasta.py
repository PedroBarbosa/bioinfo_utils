#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
remove_file = sys.argv[2] # Input wanted file, one gene name per line
result_file = sys.argv[3] # Output fasta file

remove = set()
with open(remove_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            remove.add(line)

handle=open(fasta_file, "rU")
fasta_sequences = SeqIO.parse(handle,'fasta')

output_handle = open(result_file, "w")
for record in fasta_sequences:
    if record.id not in remove and len(record.seq) > 0:
        SeqIO.write(record,output_handle, "fasta")
