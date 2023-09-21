#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv
import subprocess, re
from collections import defaultdict

class Ddict(defaultdict,dict):
    def __init__(self):
            defaultdict.__init__(self, Ddict)
    def __repr__(self):
            return dict.__repr__(self)

def sort_key(s):
    if s:
        try:
            c = re.findall('\d+$', s)[0]
        except:
            c = -1
        return int(c)

def parse_minimap2(strain_name, ref_file, minimap_result_file):
    cmd = f'minimap2 -xasm5 {ref_file} ../../fa_245_sidr_filter_genomes_spades3.15_ntref/{strain_name}_clean.fas > {minimap_result_file}'
    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
        print('minimap2 was run successfully')
    else:
        print('minimap2 was not run successfully')

def generate_identity_perc(minimap2_result_file, identity_file):
    with open(minimap2_result_file, 'r') as f1, open(identity_file, 'w') as f2:
        f2.write('contig_identifier contig_length start_position end_position ref_contig identity_length matching_length identity_perc(%) \n')
        for line in f1:
            line = line.strip().split()
            identity_perc = round(int(line[9])/int(line[10]) * 100, 2)
            if 500 <= int(line[1]) < 1000:
                f2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t'  + line[5] + '\t'+ line[9] + '\t' + line[10] + '\t' + str(identity_perc) + '\n')
    print('generating "' + identity_file + '" file successfully')
    return

def generate_coverage_perc(identity_result_file, candidate_contigs_file):
    len_dir = {}
    position_dir = Ddict()
    with open(identity_result_file, 'r') as f1, open(candidate_contigs_file, 'w') as f2:
        f2.write('contig_identifier\tref_contig\tcoverage_perc\tcoverage_length\tcontig_length\n')
        for line in f1:
            line = line.strip().split()
            if line[1] != 'contig_length' and line[4] != 'mt' and float(line[7]) >= 70:
                len_dir[line[0]] = int(line[1])
                if not line[4] in position_dir[line[0]]:
                    position_dir[line[0]][line[4]] = []
                position_dir[line[0]][line[4]].append([int(line[2]), int(line[3])])
        #print(position_dir)

        for key1 in sorted(position_dir, key=sort_key):
            merged = defaultdict(list)
            for key2 in position_dir[key1]:
                value = position_dir[key1][key2]
                #value.sort(key=lambda x: x[0])
                value.sort()
                for interval in value:
                    if not merged[key2] or merged[key2][-1][-1] < interval[0]:
                        merged[key2].append(interval)
                    else:
                        merged[key2][-1][-1] = max(merged[key2][-1][-1], interval[-1])
            #print(merged)
                #[[2, 4], [6, 8], [7, 13]]
                #[[2, 4], [6, 13]]

            coverage_perc = 0
            cov_len = 0
            for k in merged.keys():
                coverage_len = 0
                for i in merged[k]:
                    #print(key1 + '\t' + k + '\t' + str(i[1]) +'\t' + str(i[0]))
                    single_len = i[1]-i[0]
                    coverage_len += single_len
                perc = coverage_len/len_dir[key1] * 100
                #print(key1 + '\t' + k + '\t' + str(coverage_len) +'\t' + str(perc))
                if perc > coverage_perc:
                    coverage_perc = perc
                    ref_ctg = k
                    cov_len = coverage_len
            f2.write(key1 + '\t' + ref_ctg + '\t'+ str(coverage_perc) + '\t' + str(cov_len) + '\t' + str(len_dir[key1]) + '\n')
        print('generating "' + candidate_contigs_file + '" file successfully')
        return position_dir

def get_gc_content_dir(gc_content_file):
    with open(gc_content_file, 'r') as f:
        gc_content_dir = {}
        for line in f:
            line = line.strip().split()
            if line[0] != 'contig_identifier' and 500 <= int(line[2]) < 1000:
                gc_content_dir[line[0]] = float(line[1])
        return gc_content_dir

def get_lower_identify_perc_nomatching_contigs(position_dir, gc_content_dir):
    lower_identify_perc_nomatching_contigs = []
    for key in gc_content_dir.keys():
        if key not in position_dir.keys():
            lower_identify_perc_nomatching_contigs.append(key)
    return lower_identify_perc_nomatching_contigs

def get_lower_coverage_perc_contigs(candidate_contigs_file, gc_content_dir, coverage_threshold=80, gc_average=48.19, std=12.71):
    lower_coverage_gc_contigs = []
    with open(candidate_contigs_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[2] != 'coverage_perc' and float(line[2]) <= coverage_threshold:
                lower_coverage_gc_contigs.append(line[0])

            elif  line[2] != 'coverage_perc' and float(line[2]) > coverage_threshold:
                min_gc = gc_average - 2*std
                max_gc = gc_average + 2*std
                if gc_content_dir[line[0]] < min_gc or gc_content_dir[line[0]] > max_gc:
                    lower_coverage_gc_contigs.append(line[0])
        return lower_coverage_gc_contigs

def get_short_length_contigs(strain_name):
    cmd = f'fasta_length ../../fa_245_sidr_filter_genomes_spades3.15_ntref/{strain_name}_clean.fas -o > {strain_name}_ctg_contig_length.txt'
    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
        print('contigs length was generated successfully')
    else:
        print('contigs length was not generated')
    
    contigs_len_all_dir = {}
    short_contigs = []
    with open(strain_name + '_ctg_contig_length.txt', 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
            else:
                contigs_len_all_dir[key] = int(line)
                if contigs_len_all_dir[key] < 500:
                    short_contigs.append(key)
    return short_contigs

def get_mt_contigs(mt_candidate_contigs_file):
    mt_contig = []
    with open(mt_candidate_contigs_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[0] != 'contig_identifier' and float(line[1]) >= 80:
                mt_contig.append(line[0])
        return mt_contig

def get_removed_contigs(lower_identify_perc_nomatching_contigs, lower_coverage_gc_contigs, short_contigs, mt_contig, strain_name):
    print('Success: ' + strain_name + ' lower_identify_perc_nomatching_contigs number is ' + str(len(lower_identify_perc_nomatching_contigs)))
    print('Success: ' + strain_name + ' lower_coverage_gc_contigs number is ' + str(len(lower_coverage_gc_contigs)))
    print('Success: ' + strain_name + ' short_contigs number is ' + str(len(short_contigs)))
    print('Success: ' + strain_name + ' mt_contig number is ' + str(len(mt_contig)))

    removed_contigs = lower_identify_perc_nomatching_contigs + lower_coverage_gc_contigs + short_contigs + mt_contig
    with open(removed_contigs_file, 'w') as f:
        for contig in removed_contigs:
            f.write(contig + '\n')
    return removed_contigs

def filter_genome(removed_contigs, strain_name):
    print('Success: ' + strain_name + ' remove ' + str(len(removed_contigs)) + ' contigs')
    strain_dir = {}
    with open('../../fa_245_sidr_filter_genomes_spades3.15_ntref/' + strain_name + '_clean.fas', 'r') as f1, open('../../fa_245_filtered_genomes/' + strain_name + '.fas', 'w') as f2:
        for line in f1:
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
            else:
                strain_dir[key] = line
                if key not in removed_contigs:
                    f2.write('>' + key + '\n' + strain_dir[key] + '\n')

if __name__ == "__main__":
    ref_file = argv[2]
    strain_name = argv[1]
    minimap_result_file = strain_name + '_ctg_minimap_results.txt'
    identity_result_file = strain_name + '_ctg_identity_result.txt'
    candidate_contigs_file = strain_name + '_ctg_candidate_contigs.txt'
    gc_content_file = strain_name + '_ctg_GC.txt'
    mt_candidate_contigs_file = strain_name + '_mt_candidate_contigs.txt'
    removed_contigs_file = strain_name + '_removed_contigs_identifiers.txt'


    conducted_minimap2 = parse_minimap2(strain_name, ref_file, minimap_result_file)
    identity_result = generate_identity_perc(minimap_result_file, identity_result_file)
    coverage_result = generate_coverage_perc(identity_result_file, candidate_contigs_file)

    getted_gc_content_dir = get_gc_content_dir(gc_content_file)
    getted_lower_identify_perc_nomatching_contigs = get_lower_identify_perc_nomatching_contigs(coverage_result, getted_gc_content_dir)
    getted_lower_coverage_perc_contigs = get_lower_coverage_perc_contigs(candidate_contigs_file, getted_gc_content_dir, coverage_threshold=80, gc_average=48.19, std=12.71)
    getted_short_length_contigs = get_short_length_contigs(strain_name)
    getted_mt_contigs = get_mt_contigs(mt_candidate_contigs_file)
    getted_removed_contigs = get_removed_contigs(getted_lower_identify_perc_nomatching_contigs, getted_lower_coverage_perc_contigs, getted_short_length_contigs, getted_mt_contigs, strain_name)
    get_filtered_genome = filter_genome(getted_removed_contigs, strain_name)

# python3 remove_ctg.py 160205 180197_pilon4.fas
