#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
from multiprocessing import Pool
import random, subprocess, itertools
import pandas as pd

def run_msmc2(input):
    prefix = input[0]
    pop1 = input[1]
    pop2 = input[2]
    sel = input[3]
    # get the order in vcf which sorted automatically 
    sel2 = sorted(sel)
    K = {}
    for index, sample in enumerate(sel2):
        K[sample] = index
    order = [str(K[x]) for x in sel]
    para = []
    for i1 in range(4):
        for i2 in range(4,8):
            para.append(str(K[sel[i1]]) + "-" + str(K[sel[i2]]))
    parc = ",".join(para)
    par1 = ",".join(order[:4])
    par2 = ",".join(order[4:])
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
    subprocess.check_call("msmc2_linux64bit -t 4 -s -o {2}_sh.msmc2 -I {1} {0}_chr01.input {0}_chr02.input {0}_chr03.input {0}_chr04.input".format(prefix, par1, prefix + "_" + pop1), shell = True)
    subprocess.check_call("msmc2_linux64bit -t 4 -s -o {2}_sh.msmc2 -I {1} {0}_chr01.input {0}_chr02.input {0}_chr03.input {0}_chr04.input".format(prefix, par2, prefix + "_" + pop2), shell = True)
    subprocess.check_call("msmc2_linux64bit -t 4 -s -o {0}_sh.msmc2 -I {1} {0}_chr01.input {0}_chr02.input {0}_chr03.input {0}_chr04.input".format(prefix, parc), shell = True)

if __name__ == '__main__':
    F = defaultdict(list)
    G = defaultdict(list)
    L = []
    pop = []
    with open("fa_pops.list", "r") as f1:
        for line in f1:
            line = line.strip().split()
            F[line[1]].append(line[0])
            pop.append(line[1])
    pop = list(set(pop))
    for i in itertools.combinations(pop,2):
        p1 = F[i[0]]
        p2 = F[i[1]]
        for n in range(10):
            sel1 = random.sample(p1, 4)
            sel2 = random.sample(p2, 4)
            sel = sel1 + sel2
            prefix = i[0] + "_" + i[1] + "_rep" + str(n)
            L.append([prefix, i[0], i[1], sel])
    pool = Pool(30)                         # Create a multiprocessing Pool
    pool.map(run_msmc2, L)
