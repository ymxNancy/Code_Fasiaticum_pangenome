#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess, random
from multiprocessing import Pool

def treemix(m_rep):
    command = "docker run --rm -u $(id -u):$(id -u) -v $(pwd):/root -w /root calkan/treemix:1.13 treemix -i treemix.frq.strat.gz -bootstrap -seed {0} -k 50 -m {1} -o treemix_result/outstemM{1}_rep{2} -noss"
    subprocess.call(command.format(random.randint(0,10000), m_rep[0], m_rep[1]), shell=True)

if __name__ == '__main__':
    L = []
    for i in range(0,10):
        for j in range(0,50):
            L.append([i, j])
    pool = Pool(100)
    pool.map(treemix, L)
