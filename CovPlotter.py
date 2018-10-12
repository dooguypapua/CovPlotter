# -*- coding: utf-8 -*-
# dooguypapua
import sys,os
# Force python3
if sys.version_info<(3,6,0): exit("\nYou need python 3.6 or later to run CovPlotter\n\n")
import subprocess
from CovPlotter_init import *
from CovPlotter_display import *
from CovPlotter_bedgff import *
from CovPlotter_plot import *



#***** INITIALIZATION *****#
# Arguments
dicoInit = arg_manager(sys.argv)
# Check GFF
check_gff(dicoInit)

#***** CREATE & LOAD GFF DB *****#
load_gff_db(dicoInit)

#***** GENES to TRANSCRIPTS *****#
retrieve_transcripts(dicoInit)
create_detailed_bed(dicoInit)

#***** COVERAGE *****#
launch_coverage(dicoInit)
dicoCov = parse_coverage(dicoInit)

#***** PLOTTER *****#
coverview_plot(dicoInit,dicoCov)

#***** END *****#
# cleaning(dicoInit)
