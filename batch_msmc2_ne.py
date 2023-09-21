#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
from multiprocessing import Pool
import random, subprocess
import pandas as pd

def run_msmc2(input):
    prefix = input[0]
    sel = input[1]
    with open(prefix + "_strains.list", "w") as f2:
        for i in sel:
            f2.write("output/gvcf/" + i + ".g.vcf.gz" + "\n")
    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' CombineGVCFs -R 180197_pilon4_nomt_masked.fas --variant msmc_vcf/{0}_strains.list -O msmc_vcf/{0}_combined.g.vcf".format(prefix), shell = True)
    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' GenotypeGVCFs -R 180197_pilon4_nomt_masked.fas -V msmc_vcf/{0}_combined.g.vcf -O msmc_vcf/{0}.vcf".format(prefix), shell = True)
    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' SelectVariants -R 180197_pilon4_nomt_masked.fas -V msmc_vcf/{0}.vcf -O msmc_vcf/{0}_SNP.vcf --select-type-to-include SNP".format(prefix), shell = True)
    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' VariantFiltration -V msmc_vcf/{0}_SNP.vcf --filter-expression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'my_filter' -O msmc_vcf/{0}_SNP_marked.vcf".format(prefix), shell = True)
    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/asiaticum/call_snp_220409:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' SelectVariants -R 180197_pilon4_nomt_masked.fas -V msmc_vcf/{0}_SNP_marked.vcf -O msmc_vcf/{0}_SNP_filter_gatk.vcf -select 'vc.isNotFiltered()'".format(prefix), shell = True)
    subprocess.check_call("Rscript vcf_filter_outgroup_msmc.R {0}_SNP_filter_gatk.vcf {0}_vcfR_noasterisk.vcf.gz".format(prefix), shell = True)
    df = pd.DataFrame()
    for i in sel:
        subprocess.check_call("bedtools genomecov -ibam /home/zhanghao/asiaticum/call_snp_220409/output/bam/{0}.sort.bam -d -split > mask_bed/{0}.genomecov.bed".format(i), shell= True)
        data = pd.read_csv("mask_bed/{0}.genomecov.bed".format(i), sep="\t", header=None)
        top = data[2].quantile(0.999)
        mask = data[(data[2] < 6) | (data[2] > top)]
        mask = mask.iloc[:,0:2]
        df = pd.concat([df, mask])
    df = df.drop_duplicates()
    with open("{0}.mask".format(prefix), "w") as f3:
        for index,row in df.iterrows():
            f3.write(row[0] + "\t" + str(row[1]) + "\n")
    subprocess.check_call("./prepare_msmc_input2.py {0}.vcf {0}_vcfR_noasterisk.vcf.gz {0}.mask {0}".format(prefix), shell = True)
    subprocess.check_call("msmc2_linux64bit -t 4 -o {0}.msmc2 -I 0,1,2,3,4,5,6,7 {0}_chr01.input {0}_chr02.input {0}_chr03.input {0}_chr04.input".format(prefix), shell = True)

if __name__ == '__main__':
    F = defaultdict(list)
    G = defaultdict(list)
    L = []
    with open("fa_pops.list", "r") as f1:
        for line in f1:
            line = line.strip().split()
            F[line[1]].append(line[0])
    for k,v in F.items():
        for n in range(10):
            sel = random.sample(v, 8)
            prefix = k + "_rep" + str(n)
            L.append([prefix, sel])
    pool = Pool(30)                         # Create a multiprocessing Pool
    pool.map(run_msmc2, L)  # proces strains_list iterable with pool
