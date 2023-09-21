#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
conda activate R4.1

'''

import subprocess, os, re
from sys import argv
from multiprocessing import Pool

def get_name_list(name_file):
    '''
    name_file not include SRR2466519
    name_file: strain_names_245.txt
    '''
    with open(name_file, 'r') as f:
        name_list = []
        for line in f:
            line = line.strip()
            name_list.append(line)
        return name_list

def fun_annotate(strain_name):        
        os.mkdir('/home/yangmeixin/asiaticum/annotations_pre/filtering_annotation_update_new/11_eggnog_annotation/{0}'.format(strain_name))
        os.chdir('/home/yangmeixin/asiaticum/annotations_pre/filtering_annotation_update_new/11_eggnog_annotation/{0}'.format(strain_name))
        cmd = 'emapper.py --cpu 2 -i ../../10_renamed_protein_file/{0}_protein.faa -m diamond --dmnd_iterate yes --pfam_realign realign -o {0} --sensmode sensitive'.format(strain_name)
        check = subprocess.check_call(cmd, shell=True)
        if check == 0:
            print(strain_name + ' run successfully')
        else:
            print(strain_name + ' run failed')
        subprocess.check_call('rm -r emappertmp_*', shell=True)
        os.chdir('/home/yangmeixin/asiaticum/')

def extract_fun_name(strain_name):
    fun_dir = {}
    with open('/home/yangmeixin/asiaticum/annotations_pre/filtering_annotation_update_new/11_eggnog_annotation/{0}/{0}.emapper.annotations'.format(strain_name), 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('#'):
                line = line.split('\t')
                gene_name = line[0]
                if line[20] == '-':
                    continue
                else:
                    #re_fun_name = re.search('([A-Za-z0-9_].+)', line[20])
                    #fun_dir[gene_name] = re_fun_name.group(1)
                    re_fun_name = line[20]
                    fun_dir[gene_name] = re_fun_name
        return fun_dir

def add_fun_name(strain_name, fun_dir):
    with open('/home/yangmeixin/asiaticum/annotations_pre/filtering_annotation_update_new/12_added_functional_annotation/{0}.gff3'.format(strain_name), 'w') as f1, open('/home/yangmeixin/asiaticum/annotations_pre/filtering_annotation_update_new/9_rename_annotation/{0}_final_renamed.gff3'.format(strain_name), 'r') as f2:
        for line in f2:
            line =  line.strip()
            if not line.startswith('#'):
                line_1 = line.split()
                if line_1[2] == 'mRNA' or line_1[2] == 'CDS' or line_1[2] == 'exon':
                    re_id = re.search('ID=(\w+\.t\d)', line)
                    id = re_id.group(1)
                    if id in fun_dir.keys():
                        f1.write(line + ';product=' + fun_dir[id] + '\n')
                    else:
                        f1.write(line + ';product=hypothetical protein\n')
                else:
                    f1.write(line + '\n')
            else:
                f1.write(line + '\n')

def validate_gff(strain_name, functional_gff_file):
    with open('/home/yangmeixin/asiaticum/annotations_pre/filtering_annotation_update_new/validate_error_{0}.txt'.format(strain_name), 'w') as f:

        cmd = 'gt gff3validator annotations_pre/filtering_annotation_update_new/12_added_functional_annotation/{0}'.format(functional_gff_file)
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        result = p.stderr.read().decode()
        # decode: bytes transfer to string

        for list_element in result.split('\n'):
            print(list_element)
            if 'error' in list_element:
                f.write(list_element + '\n')

if __name__ == '__main__':
    name_file = argv[1]
    pool = Pool(40)
    getted_name_list = get_name_list(name_file)
    

    for strain_name in getted_name_list:
        got_protein_seq = get_protein_seq(strain_name)

'''
    pool.map(fun_annotate, getted_name_list)
    for strain_name in getted_name_list:
        functional_gff_file = strain_name + '.gff3'

        extracted_fun_name = extract_fun_name(strain_name)
        added_fun_name = add_fun_name(strain_name, extracted_fun_name)
        validated_gff = validate_gff(strain_name, functional_gff_file)
        

    os.chdir('/home/yangmeixin/asiaticum/annotations_pre/filtering_annotation_update_new/')
    cmd_1 = 'cat validate_error_*.txt > validate_error.function.txt && rm validate_error_*.txt'
    subprocess.check_call(cmd_1, shell=True)
    os.chdir('/home/yangmeixin/asiaticum/')
'''

# command line: screen -L ./functional_annotation_update.py strain_names_245.txt 
# command line: ./functional_annotation_update.py strain_name_140028.txt 



