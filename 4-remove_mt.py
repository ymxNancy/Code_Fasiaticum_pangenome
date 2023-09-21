#!/usr/bin/env python3

from sys import argv
import subprocess

def parse_minimap2(strain_name, mt_file, minimap_result_file):
    cmd = f'minimap2 -xasm5 {mt_file} /home/zhanghao/asiaticum/fa_245_sidr_filter_genomes_spades3.15_ntref/{strain_name}_clean.fas > {minimap_result_file}'
    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
        print('minimap2 was run successfully')
    else:
        print('minimap2 was not run successfully')

def generate_identity_perc(minimap2_result_file, identity_file):
    with open(minimap2_result_file, 'r') as f1, open(identity_file, 'w') as f2:
        f2.write('contig_identifier contig_length start_position end_position identity_length matching_length identity_per(%) \n')
        for line in f1:
            line = line.strip().split()
            identity_per = round(int(line[9])/int(line[10]) * 100, 2)
            f2.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[9] + '\t' + line[10] + '\t' + str(identity_per) + '\n')
    print('generating "' + identity_file + '" file successfully')
    return

def generate_coverage_perc(identity_result_file, candidate_contigs_file):
    len_dir = {}
    position_dir = {}
    with open(identity_result_file, 'r') as f1, open(candidate_contigs_file, 'w') as f2:
        f2.write('contig_identifier\tcoverage_per(%)\tcoverage_length(bp)\tcontig_length(bp)\n')
        for line in f1:
            line = line.strip().split()
            if line[1] != 'contig_length':
                if int(line[1]) > 500 and float(line[6]) >= 70:
                    len_dir[line[0]] = int(line[1])
                    position_dir.setdefault(line[0],[]).append([int(line[2]), int(line[3])])
                    #print(position_dir)

        print(position_dir)

        for key, value in position_dir.items():
            #value.sort(key=lambda x: x[0])
            value.sort()
            merged = []
            #print(value)
            for interval in value:
                if not merged or merged[-1][-1] < interval[0]:
                    merged.append(interval)
                else:
                    merged[-1][-1] = max(merged[-1][-1], interval[-1])
            
            #print(merged)
            #[[2, 4], [6, 8], [7, 13]]
            #[[2, 4], [6, 13]]

            coverage_len = 0
            for i in merged:
                single_len = int(i[1])-int(i[0]) + 1
                coverage_len += single_len
            coverage_per = float(coverage_len/len_dir[key] * 100)
            #print(coverage_per)

            f2.write(key + '\t' + str(coverage_per) + '\t' + str(coverage_len) + '\t' + str(len_dir[key]) + '\n')
        print('generating "' + candidate_contigs_file + '" file successfully')


if __name__ == "__main__":
    mt_file = argv[1]
    strain_name = argv[2]
    minimap_result_file = strain_name + '_mt_minimap_results.txt'
    identity_result_file = strain_name + '_mt_identity_result.txt'
    candidate_contigs_file = strain_name + '_mt_candidate_contigs.txt'

    conducted_minimap2 = parse_minimap2(strain_name, mt_file, minimap_result_file)
    identity_result = generate_identity_perc(minimap_result_file, identity_result_file)
    coverage_result = generate_coverage_perc(identity_result_file, candidate_contigs_file)

# python3 remove_mt.py 180197_pilon4_mt.fas 160205
