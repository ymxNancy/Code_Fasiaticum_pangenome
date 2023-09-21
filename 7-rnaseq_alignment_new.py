#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re, subprocess
from sys import argv


def get_genome_name(genome_name_file):
    with open(genome_name_file, 'r') as f1:
        genome_list = []
        for line in f1:
            line = line.strip()
            genome_list.append(line)
        return genome_list

def rna_alignment(list):
    for strain_name in list:
        
        cmd  = f'cd /home/yangmeixin/asiaticum/RNAseq_alignment_220329 && mkdir {strain_name} && cd {strain_name} && ln -s /home/yangmeixin/asiaticum/soft_masked_genomes_220405/{strain_name}_ordered_soft.fas ./ && hisat2-build {strain_name}_ordered_soft.fas {strain_name} && hisat2 -p 8 -x {strain_name} --max-intronlen 5000 --dta -1 ../VZ20171225-03-R70644-3-pool_1.fq.gz -2 ../VZ20171225-03-R70644-3-pool_2.fq.gz -S {strain_name}-hisat2.sam && samtools view -buS -@ 8 {strain_name}-hisat2.sam | samtools sort -m 3G -@ 8 -o {strain_name}-hisat2.sorted.bam && ln -s /home/yangmeixin/asiaticum/RNAseq_alignment_220329/{strain_name}/{strain_name}-hisat2.sorted.bam /home/yangmeixin/asiaticum/RNAseq_alignment_220405/'

        check = subprocess.check_call(cmd, shell=True)
        if check == 0:
            print(strain_name + ' run successfully')
        else:
            print(strain_name + ' run failed')

if __name__ == '__main__':
    genome_file = argv[1]
    genome_name_list = get_genome_name(genome_file)
    rnaseq_align = rna_alignment(genome_name_list)

# command line: ./rnaseq_alignment_new.py ../strain_name_140028.txt



    
