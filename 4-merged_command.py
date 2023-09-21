#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess, os
from sys import argv

def get_strain_name(strain_name_file):
    strain_name_list = []
    with open(strain_name_file, 'r') as f:
        for line in f:
            line = line.strip()
            strain_name_list.append(line)
        return strain_name_list

def parse_command(strain_name_list):
    for strain_name in strain_name_list:
        cmd_1 = f'mkdir {strain_name}'
        subprocess.check_call(cmd_1, shell=True)
        os.chdir(strain_name)

        cmd_2 = f'python ../remove_mt.py /home/zhanghao/asiaticum/fa_245_sidr_filter_genomes_spades3.15_ntref/180197_pilon4_mt.fas {strain_name}'
        subprocess.check_call(cmd_2, shell=True)

        cmd_3 = f'python ../calculate_GC.py {strain_name}'
        subprocess.check_call(cmd_3, shell=True)

        cmd_4 = f'python ../remove_ctg.py {strain_name} /home/zhanghao/asiaticum/fa_245_sidr_filter_genomes_spades3.15_ntref/180197_pilon4.fas'
        check = subprocess.check_call(cmd_4, shell=True)
        if check == 0:
            print(strain_name + ' filtered successfully')
        else:
            print(strain_name + ' filtered unsuccessfully')
        os.chdir('../')

if __name__ == '__main__':
    name_file = argv[1]
    getted_strain_name = get_strain_name(name_file)
    filtering_genome = parse_command(getted_strain_name)
