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
    printcolor("\n   __      __              \n  /   _   |__)| _ |_|_ _ _ \n  \__(_)\/|   |(_)|_|_(-|  \n  ____________________________________________________________________\n\n","0","255;187;108",dicoInit['color'])
    printcolor("  Output","0","255;187;108",dicoInit['color'])
    printcolor(" <> "+dicoInit['out']+"\n","0","222;220;184",dicoInit['color'])
    printcolor("  Temp","0","255;187;108",dicoInit['color'])
    printcolor("   <> "+dicoInit['tmp']+"\n","0","222;220;184",dicoInit['color'])
    printcolor("  ____________________________________________________________________\n\n","0","255;187;108",dicoInit['color'])

#***** Print usage *****#
def printusage(dicoInit):
    printcolor("\n  USAGE: CovPlotter.py [OPTIONS] -g genes.txt -i listBam.txt -o output","0","255;187;108",dicoInit['color'])
    printcolor("\n         CovPlotter.py [OPTIONS] -n ACAD,OPA1 -i listBam.txt -o output","0","255;187;108",dicoInit['color'])
    printcolor("\n         CovPlotter.py [OPTIONS] -l genes.bed -i listBam.txt -o output","0","255;187;108",dicoInit['color'])
    printcolor("\n\n  [OPTIONS]","0","222;220;184",dicoInit['color'])
    printcolor("\n      -ref       STR      Reference hg19/hg20      (default: hg19)","0","222;220;184",dicoInit['color'])
    printcolor("\n      -db        STR      Database refseq/ensembl  (default: refseq)","0","222;220;184",dicoInit['color'])
    printcolor("\n      -lcov      INT      Low coverage threshold   (default: 50)","0","222;220;184",dicoInit['color'])
    printcolor("\n      -hcov      INT      High coverage threshold  (default: 100)","0","222;220;184",dicoInit['color'])
    printcolor("\n      -tmp       DIR      Temporary folder         (default: /tmp)","0","222;220;184",dicoInit['color'])
    printcolor("\n      -nt        INT      Threads number           (default: 1)","0","222;220;184",dicoInit['color'])    
    printcolor("\n      -color     BOOL     Terminal color           (true or false)","0","222;220;184",dicoInit['color'])
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
