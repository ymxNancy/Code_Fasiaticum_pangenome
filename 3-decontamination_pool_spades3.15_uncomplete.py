#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess, os
from multiprocessing import Pool

'''
nt_faref (blast database name): nt combined with 180197

fa_210_sidr_filter_genomes_spades3.15_nt: not include 180197 when perform blast

'''

from sys import argv

def build_index(sequence_fas_name):
    '''
    input_fn: fas file of the sequence, such as '140001.fas'
    generate five files: .amb, .ann, .bwt, .pac, .sa
    '''
#    subprocess.check_call('ln -s /home/zhanghao/asiaticum/fa_211_pre_genomes_3.15_largekmer/{0}/contigs.fasta /home/zhanghao/asiaticum/fa_210_sidr_filter_analysis_spades3.15_ntref/{0}/{0}.fas'.format(sequence_fas_name), shell=True)
    cmd = 'bwa index {0}.fas'.format(sequence_fas_name)
    subprocess.check_call(cmd, shell=True)

def parse_strain_names(name_file):
    name_dir = {}
    with open(name_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            name_dir[line[0]] = line[1]
    return name_dir

def mem_alignment(sequence_fas_name, name_dir):
    if sequence_fas_name == '180271':
        cmd = 'bwa mem {0}.fas /home/zhanghao/asiaticum/filtered/{1}_combine_rename_1.fq.gz /home/zhanghao/asiaticum/filtered/{1}_combine_rename_2.fq.gz -t 5 > {0}.sam'.format(sequence_fas_name, name_dir[sequence_fas_name])
        subprocess.check_call(cmd, shell=True)
    elif sequence_fas_name == 'SRR2466519':
        cmd = 'bwa mem {0}.fas /home/zhanghao/asiaticum/filtered/{1}_1.fastq.gz /home/zhanghao/asiaticum/filtered/{1}_2.fastq.gz -t 5 > {0}.sam'.format(sequence_fas_name, name_dir[sequence_fas_name])
        subprocess.check_call(cmd, shell=True)
    else:
        cmd = 'bwa mem {0}.fas /home/zhanghao/asiaticum/filtered/{1}_filtered_1.fq.gz /home/zhanghao/asiaticum/filtered/{1}_filtered_2.fq.gz -t 5 > {0}.sam'.format(sequence_fas_name, name_dir[sequence_fas_name])
        subprocess.check_call(cmd, shell=True)

def transfer_bam(sequence_fas_name):
    cmd = 'samtools view -b -S -@ 5 {0}.sam | samtools sort -@ 5 -m 10G -o {0}.sort.bam'.format(sequence_fas_name)
    subprocess.check_call(cmd, shell=True)


def build_bam_index(sequence_fas_name):
    cmd = 'samtools index {0}.sort.bam'.format(sequence_fas_name)
    subprocess.check_call(cmd, shell=True)


def bulid_blast(sequence_fas_name):
    cmd = 'blastn -task megablast -query {0}.fas -db nt_faref -culling_limit 5 -evalue 1e-25 -num_threads 5 -outfmt "6 std staxids" -out {0}_nt.blast'.format(sequence_fas_name)
    subprocess.check_call(cmd, shell=True)


def build_blast2(sequence_fas_name):
    '''
    cut: print selected parts of lines
    -f: select pointed colcums to print
    '''

    cmd = 'cut -f 1,13 {0}_nt.blast > {0}_nt.blast2'.format(sequence_fas_name)
    subprocess.check_call(cmd, shell=True)


def remove_contaminated(sequence_fas_name):
    cmd = 'sidr default -d /home/zhanghao/software/ncbi_db/taxdump -b {0}.sort.bam -f {0}.fas -r {0}_nt.blast2 -k {0}_tokeep.contigids -x {0}_toremove.contigids -t ascomycota'.format(sequence_fas_name)
    subprocess.check_call(cmd, shell=True)


def get_seq_dic(sequence_fas_name):
    with open(sequence_fas_name + '.fas', 'r') as f:
        total_contigs_dic = {}
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
                total_contigs_dic[key] = []
            else:
                total_contigs_dic[key].append(line)
        return total_contigs_dic

def get_keep_contigs(sequence_fas_name):
    with open('{0}_tokeep.contigids'.format(sequence_fas_name), 'r') as f:
        keep_contigs_list = []
        for line in f:
            line = line.strip()
            keep_contigs_list.append(line)
        return keep_contigs_list

def get_filtered_contig_seq(sequence_fas_name, total_contigs_dic, keep_contigs_list):
    with open('../../fa_210_sidr_filter_genomes_spades3.15_ntref/{0}_clean.fas'.format(sequence_fas_name), 'w') as f:
        for contigs_name in keep_contigs_list:
            seq = ''.join(total_contigs_dic[contigs_name])
            f.write('>' + contigs_name + '\n' + seq + '\n')


def main(sample_name):
    path = os.path.join('fa_210_sidr_filter_analysis_spades3.15_ntref', sample_name)
#    os.mkdir(path)
    os.chdir(path)
    build_index(sample_name)
    mem_alignment(sample_name, getted_strain_dir)
    transfer_bam(sample_name)
    build_bam_index(sample_name)
#    bulid_blast(sample_name)
#    build_blast2(sample_name)
    remove_contaminated(sample_name)
    # save sample fas in dictionary {sample_name:seq}
    fas_seq_dic = get_seq_dic(sample_name)
    keep_contigs_list = get_keep_contigs(sample_name)
    clean_contigs = get_filtered_contig_seq(sample_name, fas_seq_dic, keep_contigs_list)
    os.chdir('../..')

if __name__ == "__main__":
    '''
    fa_210_name_file = argv[2]
    # argv[2]: fa_210_name.txt

    strain_list = []
    # save sample and sequencing name in a dictionary
    #with open('fa_36_name.txt', 'r') as f:
    with open(fa_210_name_file, 'r') as f:
        for line in f:
            line = line.strip()
            strain_list.append(line)
    '''
    list_modified_file = argv[1]
    # argv[1]: list_modified.txt
    getted_strain_dir = parse_strain_names(list_modified_file)
    x = os.listdir('fa_210_sidr_filter_genomes_spades3.15_ntref')
    L1 = [a.replace('_clean.fas', '') for a in x]
    L2 =  os.listdir('fa_210_sidr_filter_analysis_spades3.15_ntref')
    strain_list = list(set(L2) - set(L1))
    pool = Pool(15)
    pool.map(main, strain_list)
