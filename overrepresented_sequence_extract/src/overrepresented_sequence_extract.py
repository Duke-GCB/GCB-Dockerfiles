#!/usr/bin/env python

##################################################
#  overrepresented_sequence_extract.py
#
# @ /data/reddylab/Projects/GGR/analyses/group_general/top_level_scripts/overrepresented_sequence_extract.py
#  This script takes as input a fastqc report and outputs a fasta file of overrepresented/adapter sequences.
#  The pipeline conservatively assumes all overrepresented sequences are "custom adapters", which then fed
#  into trimmomatic for their subsequenct trimming.
##################################################
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import sys
import os

if len(sys.argv) <= 3:
        print "usage: python overrepresented_sequence_extract.py sample_name se/pe qc_report_folder default_adapters adapters_out_folder"
        exit(-1)

sample_name = sys.argv[1]
qc_report_folder = sys.argv[2]
default_adapters = sys.argv[3]
adapters_out_folder = sys.argv[4]


fastqc_data = qc_report_folder + '/' + sample_name + "_fastqc/fastqc_data.txt"
custom_adapter = adapters_out_folder + '/' + sample_name + "_custom_adapters.fasta"

overrep_seq_list = []
rev_comp_list  = []
overrep_sequence_info = []

with open(fastqc_data, "rb") as fh:
    for line in fh:
        if ">>Overrepresented sequences" in line :
            while True:
                overrep_sequence = next(fh).strip("\n")
                if ">>END_MODULE" in overrep_sequence:
                    break
                else:
                    overrep_sequence_info.append(overrep_sequence)

if overrep_sequence_info:
    for line in overrep_sequence_info:
        this_sequence = line.split("\t")[0]
        this_id = line.split("\t")[3]
        this_id = "_".join(this_id.split("(")[0].strip().split(" "))
        this_id.replace(",", "")
        if this_sequence != "#Sequence":
            overrep_seq = SeqRecord(Seq(this_sequence, IUPAC.unambiguous_dna), id=this_id, description ="")
            overrep_seq_list.append(overrep_seq)

if overrep_seq_list:
    for seq_record in overrep_seq_list:
        rev_comp_seq = seq_record.reverse_complement()
        rev_comp_seq.description = ""
        rev_comp_seq.id = seq_record.id + "_rc"
        rev_comp_list.append(rev_comp_seq)

    overrep_seq_list.extend(rev_comp_list)
    output_handle = open(custom_adapter, "w")
    SeqIO.write(overrep_seq_list, output_handle, "fasta")
    output_handle.close()

output_handle = open(custom_adapter, "a")
default_adapters = list(SeqIO.parse(open(default_adapters), "fasta"))
SeqIO.write(default_adapters, output_handle, "fasta")
output_handle.close()
