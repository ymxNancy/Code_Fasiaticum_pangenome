#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, subprocess

def parse_name_file(strain_name_file):
    name_dir = {}
    with open(strain_name_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            raw_name = line[1]
            name_dir[raw_name] = line[0]
        return name_dir

def parse_qc(name_dir):
    os.chdir('/home/bac/disk4/zhanghao/ymx/Clean_data')
    for raw_name, strain_name in name_dir.items():
        if raw_name == '180243':
            cmd = 'fastp -i {0}-1_FDSW210063682-1r_1.clean.fq.gz -I {0}-1_FDSW210063682-1r_2.clean.fq.gz -o /home/bac/disk4/zhanghao/ymx/qc_data/{1}_1_qc.fastq.gz -O /home/bac/disk4/zhanghao/ymx/qc_data/{1}_2_qc.fastq.gz -w 30 -j /home/bac/disk4/zhanghao/ymx/qc_data/{1}.json -h /home/bac/disk4/zhanghao/ymx/qc_data/{1}.html'.format(raw_name, strain_name)
        elif raw_name == '171465':
            cmd = 'fastp -i {0}-1_FDSW210063680-1r_1.clean.fq.gz -I {0}-1_FDSW210063680-1r_2.clean.fq.gz -o /home/bac/disk4/zhanghao/ymx/qc_data/{1}_1_qc.fastq.gz -O /home/bac/disk4/zhanghao/ymx/qc_data/{1}_2_qc.fastq.gz -w 30 -j /home/bac/disk4/zhanghao/ymx/qc_data/{1}.json -h /home/bac/disk4/zhanghao/ymx/qc_data/{1}.html'.format(raw_name, strain_name)

        elif raw_name == '180104':
            cmd = 'fastp -i P{0}-1_FDSW210063685-1r_1.clean.fq.gz -I P{0}-1_FDSW210063685-1r_2.clean.fq.gz -o /home/bac/disk4/zhanghao/ymx/qc_data/{1}_1_qc.fastq.gz -O /home/bac/disk4/zhanghao/ymx/qc_data/{1}_2_qc.fastq.gz -w 30 -j /home/bac/disk4/zhanghao/ymx/qc_data/{1}.json -h /home/bac/disk4/zhanghao/ymx/qc_data/{1}.html'.format(raw_name, strain_name)
        elif raw_name == 'R198':
            cmd = 'fastp -i {0}_combine_rename_1.fq.gz -I {0}_combine_rename_2.fq.gz -o /home/bac/disk4/zhanghao/ymx/qc_data/{1}_1_qc.fastq.gz -O /home/bac/disk4/zhanghao/ymx/qc_data/{1}_2_qc.fastq.gz -w 30 -j /home/bac/disk4/zhanghao/ymx/qc_data/{1}.json -h /home/bac/disk4/zhanghao/ymx/qc_data/{1}.html'.format(raw_name, strain_name)
        else:
            cmd = 'fastp -i Fungal_H720-{0}_good_1.fq.gz -I Fungal_H720-{0}_good_1.fq.gz -o /home/bac/disk4/zhanghao/ymx/qc_data/{1}_1_qc.fastq.gz -O /home/bac/disk4/zhanghao/ymx/qc_data/{1}_2_qc.fastq.gz -w 30 -j /home/bac/disk4/zhanghao/ymx/qc_data/{1}.json -h /home/bac/disk4/zhanghao/ymx/qc_data/{1}.html'.format(raw_name, strain_name)
        subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    strain_name_file = '/home/bac/disk4/zhanghao/ymx/list_all.txt'
    strain_name_dir = parse_name_file(strain_name_file)
    perform_qc = parse_qc(strain_name_dir)