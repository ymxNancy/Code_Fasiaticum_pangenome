#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv
import re

def get_contig_dir(strain_name):
    contig_dir = {}
    strain_path = '../../fa_245_sidr_filter_genomes_spades3.15_ntref/' + strain_name + '_clean.fas'
    with open(strain_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                key = line.replace('>', '')
                contig_dir[key] = []
            else:
                contig_dir[key].append(line)
        return contig_dir

def get_GC_content(strain_name, contig_dir):
    def sort_key(s):
        if s:
            try:
                c = re.findall('\d+$', s)[0]
            except:
                c = -1
            return int(c)

    with open(strain_name + '_ctg_GC.txt', 'w') as f:
        #res.write('contig_identifier' + '\t' + 'GC_rate' + '\t' + 'ctg_length' + '\n')
        f.write("contig_identifier\tGC_rate\tctg_length\n")
        for k in sorted(contig_dir.keys(), key = sort_key):
            seq = ''.join(contig_dir[k])
            gc = (seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c'))/len(seq)*100
            #f.write("{}\t{:.2f}\t{}\n").format(k, str(gc), str(len(seq)))
            f.write(k + '\t' + str(gc) + '\t' + str(len(seq)) + '\n')

if __name__ == '__main__':
    strain_name = argv[1]
    getted_contig_dir = get_contig_dir(strain_name)
    getted_GC_content = get_GC_content(strain_name, getted_contig_dir)
