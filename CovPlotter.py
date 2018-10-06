# -*- coding: utf-8 -*-
# dooguypapua
import sys,os,shutil
import subprocess
from CovPlotter_init import *
from CovPlotter_display import *


#***** INITIALIZATION *****#
# Arguments
dicoInit = arg_manager(sys.argv)
# Check GFF
if dicoInit["gff"]=="": get_refseq_gff(dicoInit)


#***** CREATE GFF DB *****#
printcolor("\n    • Create GFF DB","0","255;187;108",dicoInit['color'])
if len(dicoInit['lst_genes'])>0:
    dicoInit['genes'] = dicoInit['tmp']+"/grep_genes.txt"
    O = open(dicoInit['genes'],'w')
    for gene in dicoInit['lst_genes']: O.write("gene="+gene+";\n")
    O.close()

# fn = gffutils.example_filename(dicoNiourk["refseq_gff"])
# gffutils.create_db(fn, dbfn=dicoNiourk["refseq_gff"].replace(".gff",".db"), force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)


# print(dicoInit)
# shutil.rmtree(dicoInit['tmp'])
exit("\n\n")