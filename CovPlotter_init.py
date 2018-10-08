# -*- coding: utf-8 -*-
# dooguypapua
import sys,os
import subprocess
from binaryornot.check import is_binary
from yaspin import yaspin
from yaspin.spinners import Spinners
import urllib.request
import gzip
import uuid
from CovPlotter_display import *


#***** Arg manager *****#
def arg_manager(lst_arg):
    dicoInit = { 'inputType':"", 'lst_bam':[], 'genes':"", 'lst_genes':[], 'bed':"", 'gff':"", "db_gff":None, 'out':"", 'ref':"", 'spinner':yaspin(), 'tmp':"/tmp/"+str(uuid.uuid4())[:8], 'nt':1, 'color':True }
    lstError = []
    bool_i = False ; bool_o = False ; bool_g = False ; bool_n = False ; bool_l = False ; bool_r = False
    if len(lst_arg)==1 or lst_arg[1] in ["-h","--help"]: printtitle(dicoInit) ; printusage(dicoInit)
    if len(lst_arg)==1 or lst_arg[1] in ["-v","--version"]: printtitle(dicoInit) ; printcolor("  version 1.0","0","255;187;108",dicoInit['color']) ; exit("\n\n")
    iterarg = iter(list(range(len(lst_arg)))) ; next(iterarg)
    for i in iterarg:
        # REQUIRED PARAMETERS
        if lst_arg[i]=="-i":
            bool_i = True
            if len(lst_arg)<=i+1: lstError.append("(-i) Any value specified")
            elif not os.path.isfile(lst_arg[i+1]): lstError.append("(-i) File \""+lst_arg[i+1]+"\" was not found")
            else:
                path_lst_bam = lst_arg[i+1]
                if is_binary(path_lst_bam): lstError.append("(-i) Invalid file \""+lst_arg[i+1]+"\" (binary found)")
                else: 
                    try: I = open(path_lst_bam,'r') ; lst_bam = I.read().split("\n") ; I.close()
                    except: lstError.append("(-i) Unable to read file \""+lst_arg[i+1]+"\"")
                    for bam in lst_bam:
                        if os.path.isfile(bam):
                            pipe = subprocess.Popen("samtools quickcheck "+bam, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) ; pipe.communicate()
                            if pipe.returncode!=0: lstError.append("(-i) Invalid file format \""+bam+"\"")
                            else: dicoInit['lst_bam'].append(bam)
                        elif bam!="": lstError.append("(-i) File \""+bam+"\" was not found")
            next(iterarg)
        elif lst_arg[i]=="-o":
            bool_o = True
            if len(lst_arg)<=i+1: lstError.append("(-o) Any value specified")
            else:
                try: os.makedirs(lst_arg[i+1], exist_ok=True) ; dicoInit['out'] = lst_arg[i+1]
                except: lstError.append("(-o) Failed to create folder \""+lst_arg[i+1]+"\"")
            next(iterarg)
        elif lst_arg[i]=="-g":
            bool_g = True
            if len(lst_arg)<=i+1: lstError.append("(-g) Any value genes file specified")
            elif not os.path.isfile(lst_arg[i+1]): lstError.append("(-g) File \""+lst_arg[i+1]+"\" was not found")
            else:
                dicoInit['genes'] = lst_arg[i+1]
                if is_binary(dicoInit['genes']): lstError.append("(-g) Invalid file \""+lst_arg[i+1]+"\" (binary found)")
                else: 
                    try: I = open(dicoInit['genes'],'r') ; dicoInit['lst_genes'] = I.read().split("\n") ; I.close()
                    except: lstError.append("(-g) Unable to read file \""+lst_arg[i+1]+"\"")
            dicoInit['inputType'] = "g"
            next(iterarg)
        elif lst_arg[i]=="-n":
            bool_n = True
            if len(lst_arg)<=i+1: lstError.append("(-n) Any value specified")
            else: dicoInit['lst_genes'] = lst_arg[i+1].split(",")
            dicoInit['inputType'] = "n"
            next(iterarg)
        elif lst_arg[i]=="-l":
            bool_l = True
            if len(lst_arg)<=i+1: lstError.append("(-l) Any value specified")
            elif not os.path.isfile(lst_arg[i+1]): lstError.append("(-l) File \""+lst_arg[i+1]+"\" was not found")
            else:
                dicoInit['bed'] = lst_arg[i+1]
                if is_binary(dicoInit['bed']): lstError.append("(-l) Invalid file \""+lst_arg[i+1]+"\" (binary found)")
            dicoInit['inputType'] = "l"
            next(iterarg)
        elif lst_arg[i]=="-r": 
            bool_r = True
            if len(lst_arg)<=i+1: lstError.append("(-r) Any value specified")
            elif not lst_arg[i+1].lower() in ["hg19","hg20"]: lstError.append("(-r) Invalid value \""+lst_arg[i+1]+"\"")
            else: dicoInit['ref'] = lst_arg[i+1].lower()
            next(iterarg)
        # OPTIONS PARAMETERS
        elif lst_arg[i]=="-gff":
            if len(lst_arg)<=i+1: lstError.append("(-gff) Any value specified")
            elif not os.path.isfile(lst_arg[i+1]): lstError.append("(-gff) File \""+lst_arg[i+1]+"\" was not found")
            else:
                dicoInit['gff'] = lst_arg[i+1]
                if dicoInit['gff'][-3:]!="gff": lstError.append("(-l) Invalid file format \""+lst_arg[i+1]+"\" ("+dicoInit['gff'][-3:]+" found)")
                elif is_binary(dicoInit['gff']): lstError.append("(-l) Invalid file format \""+lst_arg[i+1]+"\" (binary found)")
            next(iterarg)
        elif lst_arg[i]=="-tmp":
            if len(lst_arg)<=i+1: lstError.append("(-tmp) Any value specified")
            else:
                try: os.makedirs(lst_arg[i+1]+"/"+str(uuid.uuid4())[:8], exist_ok=True) ; dicoInit['tmp'] = lst_arg[i+1]
                except: lstError.append("(-tmp) Failed to create folder \""+lst_arg[i+1]+"\"")
            next(iterarg)
        elif lst_arg[i]=="-nt":
            if len(lst_arg)<=i+1: lstError.append("(-nt) Any value specified")
            else:
                try: dicoInit['nt'] = int(lst_arg[i+1])
                except: lstError.append("(-nt) Invalid value \""+lst_arg[i+1]+"\"")
            next(iterarg)
        elif lst_arg[i]=="-color":
            if len(lst_arg)<=i+1: lstError.append("(-color) Any value specified")
            elif not lst_arg[i+1].lower() in ["true","false"]: lstError.append("(-color) Invalid value \""+lst_arg[i+1]+"\"")
            else: dicoInit['color'] = bool(lst_arg[i+1])
            next(iterarg)
        else: lstError.append("Invalid parameters \""+lst_arg[i]+"\"")
    # Check required parameters
    missing = ""
    if not bool_i: missing+="\"-i\" "
    if not bool_o: missing+="\"-o\" "
    if not bool_r: missing+="\"-r\" "
    if not bool_g and not bool_n and not bool_l: missing+="\"-g|-n|-l\" "
    if len(missing)>0: lstError.insert(0,"Missing required parameters "+missing)
    if sum([bool_g, bool_n, bool_l])>1: lstError.insert(0,"Multiple genes input parameters found, select \"-g\" or \"-n\" or \"-l\"")
    # Title
    printtitle(dicoInit)
    # Errors
    if len(lstError)>0:
        sys.stdout.write('\x1b[1K')
        printcolor("\n  ERROR:","0","212;0;0",dicoInit['color'])
        printcolor(" "+lstError[0]+"\n","0","212;0;0",dicoInit['color'])
        for error in lstError[1:]: printcolor("         "+error+"\n","0","212;0;0",dicoInit['color'])
        printusage(dicoInit)
    # Create tmp folder
    os.makedirs(dicoInit['tmp'], exist_ok=True)
    return dicoInit