# -*- coding: utf-8 -*-
# dooguypapua
import sys
import time

#***** Print color *****#
def printcolor(text,style,fg_rgb,colorBool):
    # STYLE : normal (0),bold (1),light (2),italic (3),underline (4),slow blink (5),rapid blink (6), crossed-out (9)
    if colorBool==True: print("\x1b["+style+";38;2;"+fg_rgb+"m"+text+"\x1b[0m",end='')
    else: print(text,end='')

#***** Print title *****#
def printtitle(dicoInit):
    printcolor("\n   __      __              \n  /   _   |__)| _ |_|_ _ _ \n  \__(_)\/|   |(_)|_|_(-|  \n  ____________________________________________________________________\n","0","255;187;108",dicoInit['color'])

#***** Print usage *****#
def printusage(dicoInit):
    printcolor("\n  USAGE: CovPlotter.py [OPTIONS] -g genes.txt -r hg19 -i listBam.txt -o output","0","255;187;108",dicoInit['color'])
    printcolor("\n         CovPlotter.py [OPTIONS] -n ACAD,OPA1 -r hg20 -i listBam.txt -o output","0","255;187;108",dicoInit['color'])
    printcolor("\n         CovPlotter.py [OPTIONS] -l genes.bed -r hg19 -i listBam.txt -o output","0","255;187;108",dicoInit['color'])
    printcolor("\n\n  [OPTIONS]","0","222;220;184",dicoInit['color'])
    printcolor("\n      -gff       FILE     Annotation file  (gff format)","0","222;220;184",dicoInit['color'])
    printcolor("\n                          (else download and save in CovPlotter folder)","0","117;116;97",dicoInit['color'])
    printcolor("\n      -tmp       DIR      Temporary folder (default: /tmp)","0","222;220;184",dicoInit['color'])
    printcolor("\n      -nt        INT      Threads number   (default: 1)","0","222;220;184",dicoInit['color'])    
    printcolor("\n      -color     BOOL     Terminal color   (true or false)","0","222;220;184",dicoInit['color'])
    printcolor("\n      -h,--help","0","222;220;184",dicoInit['color'])
    printcolor("\n      -v,--version","0","222;220;184",dicoInit['color'])
    exit("\n\n")

#***** Download progression *****#
def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = min(int(count*block_size*100/total_size),100)
    sys.stdout.write("\r      ...%d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()
