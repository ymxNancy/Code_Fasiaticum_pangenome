#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from sys import argv


def parse_strain_names(name_file):
    name_dir = {}
    with open(name_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            name_dir[line[1]] = line[0]
    return name_dir


def perform_assembly(name_dir):
    for k,v in name_dir.items():
        if k == 'R198':
            cmd = f'spades.py -o fa_211_pre_genomes_3.15_largekmer/{v} -t 8 -m 120 -1 /home/zhanghao/asiaticum/filtered/R198_combine_rename_1.fq.gz -2 /home/zhanghao/asiaticum/filtered/R198_combine_rename_2.fq.gz -k 31,61,91,121'
            res = subprocess.check_call(cmd, shell=True)
            if res == 0:
                print(v + ' run successfully')
            else:
                print(v + ' run failed')

        elif k == 'SRR2466519':
            cmd = f'spades.py -o fa_211_pre_genomes_3.15_largekmer/{v} -t 8 -m 120 -1 /home/zhanghao/asiaticum/filtered/{k}_1.fastq.gz -2 /home/zhanghao/asiaticum/filtered/{k}_2.fastq.gz -k 31,61,91,121'
            res = subprocess.check_call(cmd, shell=True)
            if res == 0:
                print(v + ' run successfully')
            else:
                print(v + ' run failed')
        else:
            cmd = f'spades.py -o fa_211_pre_genomes_3.15_largekmer/{v} -t 8 -m 120 -1 /home/zhanghao/asiaticum/filtered/{k}_filtered_1.fq.gz -2 /home/zhanghao/asiaticum/filtered/{k}_filtered_2.fq.gz -k 31,61,91,121'
            res = subprocess.check_call(cmd, shell=True)
            if res == 0:
                print(v + ' run successfully')
            else:
                print(v + ' run failed')


if __name__ == "__main__":
    strain_file = argv[1]
    getted_name_list = parse_strain_names(strain_file)
    generate_genomes = perform_assembly(getted_name_list)


