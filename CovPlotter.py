# -*- coding: utf-8 -*-
# dooguypapua
import sys,os,shutil
import subprocess
from CovPlotter_init import *
from CovPlotter_display import *
from CovPlotter_bedgff import *



#***** INITIALIZATION *****#
# Arguments
dicoInit = arg_manager(sys.argv)
# Check GFF
if dicoInit["gff"]=="": get_gff(dicoInit)


#***** CREATE & LOAD GFF DB *****#
create_gff_db(dicoInit)
load_gff_db(dicoInit)

#***** GENES to TRANSCRIPTS *****#
retrieve_transcripts(dicoInit)

# print(dicoInit)
# shutil.rmtree(dicoInit['tmp'])
exit("\n\n")