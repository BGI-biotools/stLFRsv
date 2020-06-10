#!/bin/bash
./LFR-sv -bam example/L0.sort.rmdup.bam -out example/result -ncpu 20 -phase example/phase_out -bl data/bad_region_hg19_withchr.bllist -cl data/human_hg19_2000_20000_20000_0.9995_0.95_withchr.conlist -human Y
