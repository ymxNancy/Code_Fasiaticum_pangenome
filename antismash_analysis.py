#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from multiprocessing.pool import Pool
from sys import argv 

'''
conda activate antismash

'''

def parse_strain_name(strain_name_file):
    with open(strain_name_file, 'r') as f:
        name_list = []
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list

def parse_antismash(strain_name):
    if strain_name == '180197':
        fasta_file = strain_name + '_pilon4_nomt_masked.fas'
    else:
        fasta_file = strain_name + '_ordered_soft.fas'

    subprocess.check_call("~/bin/run_antismash {0} 13_antismash_result --taxon fungi --cassis --fullhmmer --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --cc-mibig --genefinding-gff3 /input/{1}_without_locus.gff3".format(fasta_file, strain_name), shell=True)

if __name__ == '__main__':
    name_list_file = argv[1]
    getted_name_list = parse_strain_name(name_list_file)

    pool = Pool(7)
    antismash_analysis = pool.map(parse_antismash, getted_name_list)

