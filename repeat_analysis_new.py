#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from sys import argv
from multiprocessing import Pool

def get_genomes_name(genome_name_file):
    name_list = []
    with open(genome_name_file, 'r') as f:
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list

def conduct_repeat_analysis(strain_name):
    if strain_name == '180197':
        cmd = 'mkdir analysis_fa_245_filtered_genomes_softmask/180197 && cp /home/yangmeixin/asiaticum/fa_245_sidr_filter_genomes_spades3.15_ntref/180197_pilon4_nomt.fas analysis_fa_245_filtered_genomes_softmask/180197/ && cd analysis_fa_245_filtered_genomes_softmask/180197 && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 BuildDatabase -name 180197 -engine ncbi 180197_pilon4_nomt.fas && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 RepeatModeler -database 180197 -pa 5 -LTRStruct && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 RepeatMasker 180197_pilon4_nomt.fas -pa 5 -e ncbi -lib RM*/consensi.fa.classified -dir result_dir -xsmall -gff -html && cp result_dir/*.masked ../../soft_masked_genomes/180197_pilon4_nomt_soft.fas'
        check = subprocess.check_call(cmd, shell=True)
        if check == 0:
            print(strain_name + ' run successfully')
        else:
            print(strain_name + ' run unsuccessfully')

    else:
        cmd = 'mkdir analysis_fa_245_filtered_genomes_softmask/{0} && cp {0}_ordered.fas analysis_fa_245_filtered_genomes_softmask/{0}/ && cd analysis_fa_245_filtered_genomes_softmask/{0} && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 BuildDatabase -name {0} -engine ncbi {0}_ordered.fas && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 RepeatModeler -database {0} -pa 2 -LTRStruct && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 RepeatMasker {0}_ordered.fas -pa 2 -e ncbi -lib RM*/consensi.fa.classified -dir result_dir -xsmall -gff -html && cp result_dir/*.masked ../../soft_masked_genomes/{0}_ordered_soft.fas'.format(strain_name)

        check = subprocess.check_call(cmd, shell=True)
        if check == 0:
            print(strain_name + ' run successfully')
        else:
            print(strain_name + ' run unsuccessfully')


def main(strain_name):
    repeat_analysis = conduct_repeat_analysis(strain_name)

if __name__ == '__main__':
    genome_name_file = argv[1]
    getted_genomes_name = get_genomes_name(genome_name_file)
    
    pool = Pool(62)
    pool.map(main, getted_genomes_name)

    # command line: screen -L ./repeat_analysis_new.py strain_name.txt
    # screen -L ./repeat_analysis_new.py strain_name_140028.txt