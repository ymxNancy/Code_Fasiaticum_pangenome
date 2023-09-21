#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

pop = ["POP2_POP1", "POP2_POP3", "POP3_POP1"]

for i in pop:
    p1 = i.split("_")[0]
    p2 = i.split("_")[1]
    for n in range(10):
        input1 = i + "_rep" + str(n) + "_sh.msmc2.final.txt"
        input2 = i + "_rep" + str(n) + "_" + p1 + "_sh.msmc2.final.txt"
        input3 = i + "_rep" + str(n) + "_" + p2 + "_sh.msmc2.final.txt"
        output = i + "_rep" + str(n) + "_combined_sh.msmc2.final.txt"
        cmd = " ".join(["./combineCrossCoal.py", input1, input2, input3, ">", output])
        subprocess.check_call(cmd, shell = True)
