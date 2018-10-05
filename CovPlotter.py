# -*- coding: utf-8 -*-
# dooguypapua
import sys,os
import time
import urllib.request
import gzip
from yaspin import yaspin
from yaspin.spinners import Spinners


#***** Print color *****#
def printcolor(text,style,fg_rgb,colorBool):
    # STYLE : normal (0),bold (1),light (2),italic (3),underline (4),slow blink (5),rapid blink (6), crossed-out (9)
    if colorBool==True: print("\x1b["+style+";38;2;"+fg_rgb+"m"+text+"\x1b[0m",end='')
    else: print(text,end='')

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

#***** Download GFF *****#
def get_refseq_gff(ref,tmp_dir,colorBool):
    printcolor("\n    • Download refseq gff file","0","255;187;108",colorBool) ; sys.stdout.write('\x1b7')
    print("\x1b[0;38;2;159;108;102m\n") ; sys.stdout.write("\033[F")
    if ref=="hg19": url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz"
    elif ref=="hg20": url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz"
    path_tmp_download = tmp_dir+"/"+os.path.basename(url)
    urllib.request.urlretrieve(url,path_tmp_download,reporthook)
    # Extract gzip and change chromosomes name
    GZ = gzip.open(path_tmp_download, 'rb')
    GFF = open(os.getcwd()+"/db/"+os.path.basename(url).replace(".gz",""), 'wb')
    GFF.write(GZ.read())
    GZ.close()
    GFF.close()
    sys.stdout.write('\x1b[1K')
    sys.stdout.write('\x1b8')
    printcolor(" [OK]\n","0","222;220;184",colorBool)







spinner = yaspin()
colorBool = True
tmp_dir = "/tmp"

printcolor("\n   __      __              \n  /   _   |__)| _ |_|_ _ _ \n  \__(_)\/|   |(_)|_|_(-|  \n\n","1","255;187;108",colorBool)

# Help argument
printcolor("\n  USAGE: CovPlotter.py [OPTIONS] -g genes.txt -i file.bam file2.bam","0","255;187;108",colorBool)
printcolor("\n         CovPlotter.py [OPTIONS] -n ACAD,OPA1 -i file.bam file2.bam","0","255;187;108",colorBool)
printcolor("\n         CovPlotter.py [OPTIONS] -b genes.bed -i file.bam file2.bam","0","255;187;108",colorBool)
printcolor("\n\n  [OPTIONS]","0","222;220;184",colorBool)
printcolor("\n      -ref     STR      Reference genome (hg19 or hg20)","0","222;220;184",colorBool)
printcolor("\n      -gff     FILE     Genes annotation GFF file name (default: download)","0","222;220;184",colorBool)
printcolor("\n      -tmp     DIR      Temporary folder (default: /tmp)","0","222;220;184",colorBool)

os.makedirs(os.getcwd()+"/db", exist_ok=True)
# get_refseq_gff("hg19",tmp_dir,colorBool)



exit("\n\n")