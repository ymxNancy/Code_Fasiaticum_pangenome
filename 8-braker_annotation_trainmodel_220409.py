#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from sys import argv
from multiprocessing import Pool

def get_name_list(name_file):
    '''
    name_file not include 180197, SRR2466519
    name_file: strain_names_244.txt
    '''
    with open(name_file, 'r') as f:
        name_list = []
        for line in f:
            line = line.strip()
            name_list.append(line)
        return name_list

def conduct_annotation(genome):
    cmd = f'braker.pl --genome=/home/yangmeixin/asiaticum/soft_masked_genomes_220405/{genome}_ordered_soft.fas --prot_seq=odb_fungi_RR1_tri.fa --etpmode --species={genome}_model --softmasking --bam=/home/yangmeixin/asiaticum/RNAseq_alignment_220405/{genome}-hisat2.sorted.bam --gff3 --cores 65 --fungus --workingdir={genome}'
    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
        print(genome + ' run successfully')
    else:
        print(genome + ' run failed')


if __name__ == '__main__':
    pool = Pool(9)
    name_file = argv[1]

    list_name = get_name_list(name_file)
    #annotate = conduct_annotation(list_name)
    pool.map(conduct_annotation, list_name)
    
