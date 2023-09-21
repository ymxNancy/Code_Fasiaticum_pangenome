#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from multiprocessing import Pool

'''
work dictionary: /home/zhanghao/asiaticum/call_snp_220409

before running, motify ulimit -u 31308 or larger
command line: screen -L python3 bin/call_snp.py

use docker run gatk: 
    1. docker run -it -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd broadinstitute/gatk:latest bash
    2. docker run -it -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd broadinstitute/gatk:latest bash

    docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx4g' GenomicsDBImport --genomicsdb-workspace-path my_database -R 180197_pilon4.fas -L intervals.list --sample-name-map output/input.list.test

PH-1 Illumina data download from: https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR15032575&display=download

'''

def prepare_ref(ref_fa):
    subprocess.call("bwa-mem2 index {}".format(ref_fa), shell=True)
    print('generate reference index successfully')
    subprocess.call("samtools faidx {}".format(ref_fa), shell=True)
    print('generate reference .fai file successfully')

    subprocess.call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk CreateSequenceDictionary -R {} -O 180197_pilon4_nomt_masked.dict".format(ref_fa), shell=True)
    print('generate reference .dict file successfully')

def generate_strain_list(list_file):
    strains_list = []
    with open(list_file, 'r') as f:
        for line in f:
            line = line.strip()
            strains_list.append(line)
    return strains_list

def get_vcf(strain_names): 
    line = strain_names.strip().split()
    subprocess.call("bwa-mem2 mem -t 80 -M -R '@RG\\tID:{0}\\tLB:{1}\\tPL:illumina\\tSM:{0}\\tPU:{0}' 180197_pilon4_nomt_masked.fas ../filtered/{1}_filtered_1.fq.gz ../filtered/{1}_filtered_2.fq.gz > output/bam/{0}.sam".format(line[0], line[1]), shell=True)

    subprocess.call("samtools view -buS output/bam/{0}.sam | samtools sort -m 500G -o output/bam/{0}.sort.bam".format(line[0]), shell=True)

    subprocess.call("samtools index output/bam/{}.sort.bam".format(line[0]), shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx4g' FixMateInformation -I output/bam/{0}.sort.bam -O output/bam/{0}.fix.bam -MC true".format(line[0]), shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx4g' MarkDuplicates -I output/bam/{0}.fix.bam -O output/bam/{0}.mkd.bam -ASO coordinate -M output/bam/{0}_marked_dup_metrics.txt".format(line[0]), shell=True)

    subprocess.check_call("samtools index output/bam/{}.mkd.bam".format(line[0]), shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx4g' HaplotypeCaller -R 180197_pilon4_nomt_masked.fas -I output/bam/{0}.mkd.bam -ERC GVCF -ploidy 1 -O output/gvcf/{0}.g.vcf.gz".format(line[0]), shell=True)

    #   subprocess.call("rm output/bam/{0}.sam output/bam/{0}.sort.bam output/bam/{0}.fix.bam".format(line[0]), shell=True)
    print(line[0] + ' runs gvcf successfully')

def generate_vcf():
# 1st time use fm and hn9-1 as outgroups:
#     subprocess.check_call("find output/gvcf/ -name '*.g.vcf.gz' > output/add_outgroup/vcf.list", shell=True)

# 2nd time use fm and PH-1 as outgroups:
#    subprocess.check_call("find output/gvcf/ -name '*.g.vcf.gz' > output/add_outgroup_ph1_180197/vcf.list", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v $(pwd):/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' CombineGVCFs -R 180197_pilon4_nomt_masked.fas --variant output/add_outgroup_ph1_180197/vcf.list -O output/add_outgroup_ph1_180197/cohort_outgroup.g.vcf", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v $(pwd):/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' GenotypeGVCFs -R 180197_pilon4_nomt_masked.fas -V output/add_outgroup_ph1_180197/cohort_outgroup.g.vcf -O output/add_outgroup_ph1_180197/fa_com_var_outgroup.vcf", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' VariantFiltration -V output/add_outgroup_ph1_180197/fa_com_var_outgroup.vcf --filter-expression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'my_filter' -O output/add_outgroup_ph1_180197/fa_com_var_outgroup_marked.vcf", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' SelectVariants -R 180197_pilon4_nomt_masked.fas -V output/add_outgroup_ph1_180197/fa_com_var_outgroup_marked.vcf  -O output/add_outgroup_ph1_180197/fa_com_var_filtered_outgroup_gatk.vcf -select 'vc.isNotFiltered()'", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v $(pwd):/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' SelectVariants -R 180197_pilon4_nomt_masked.fas -V output/add_outgroup_ph1_180197/fa_com_var_filtered_outgroup_gatk.vcf  -O output/add_outgroup_ph1_180197/fa_com_snp_filtered_outgroup_gatk.vcf --select-type-to-include SNP", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v $(pwd):/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' SelectVariants -R 180197_pilon4_nomt_masked.fas -V output/add_outgroup_ph1_180197/fa_com_var_filtered_outgroup_gatk.vcf  -O output/add_outgroup_ph1_180197/fa_com_indel_filtered_outgroup_gatk.vcf --select-type-to-include INDEL", shell=True)

    print('finished')

if __name__ == '__main__':
    ref_fa_file = '180197_pilon4_nomt_masked.fas'
# add outgroup do not need to re do this step:
#     prepared_ref = prepare_ref(ref_fa_file)

# 1st use HN9-1:
#    strain_name_file = '/home/zhanghao/asiaticum/qc_name_outgroup.txt'
#    pool = Pool(41)
#    pool.map(get_vcf, generated_strain_list)


# 2nd time replace HN9-1 by PH-1:
#    strain_name_file = '/home/zhanghao/asiaticum/qc_name_outgroup_ph1.txt'
#    generated_strain_list = generate_strain_list(strain_name_file)
#    pool = Pool(1)
#    pool.map(get_vcf, generated_strain_list)

    generated_vcf = generate_vcf()




# CombineGVCFs --variant file name must be vcf.list

# Then using vcfR filter the DP, using vcf_filter.R script in R.
# Generate file: fa_com_snp_filtered_outgroup_gatk_vcfR.vcf

# Then using vcftools filter the max missing:
# vcftools --vcf fa_com_snp_filtered_outgroup_gatk_vcfR.vcf --max-missing 1 --recode --recode-INFO-all --out fa_com_snp_filtered_outgroup_gatk_vcfR_vcftools
#generate file: fa_com_snp_filtered_outgroup_gatk_vcfR_vcftools.recode.vcf

# vcftools --vcf fa_com_snp_filtered_outgroup_gatk_vcfR_noasterisk.vcf  --max-missing 1.0 --recode --recode-INFO-all --out fa_com_snp_filtered_outgroup_gatk_vcfR_noasterisk_maxmissing1

