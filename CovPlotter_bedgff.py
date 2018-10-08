# -*- coding: utf-8 -*-
# dooguypapua
import sys,os,shutil
import subprocess
import gffutils
from pybedtools import BedTool,cleanup
from CovPlotter_init import *
from CovPlotter_display import *


#**************************************#
#            DOWNLOAD GFF              #
#**************************************#
def get_gff(dicoInit):
    printcolor("    • Download GFF file\n","0","222;220;184",dicoInit["color"])    
    print("\x1b[0;38;2;222;220;184m") ; sys.stdout.write("\033[F")
    if dicoInit['ref']=="hg19": url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz"
    elif dicoInit['ref']=="hg20": url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz"
    path_tmp_download = dicoInit['tmp']+"/"+os.path.basename(url)
    urllib.request.urlretrieve(url,path_tmp_download,reporthook)
    # Extract gzip and change chromosomes name
    os.makedirs(os.getcwd()+"/db", exist_ok=True) # Create db folder
    dicoInit["gff"] = os.getcwd()+"/db/"+os.path.basename(url).replace(".gz","")
    GZ = gzip.open(path_tmp_download, 'rb')
    GFF = open(dicoInit["gff"], 'wb')
    GFF.write(GZ.read())
    GZ.close()
    GFF.close()
    sys.stdout.write("\033[F")
    printcolor("    • GFF downloaded","0","222;220;184",dicoInit["color"])
    printcolor(" <> "+dicoInit["gff"]+"\n","0","117;116;97",dicoInit["color"])


#**************************************#
#       CREATE gffutils Database       #
#**************************************#
def create_gff_db(dicoInit):
    if not os.path.isfile(dicoInit["gff"].replace(".gff",".db")):
        print("\x1b[0;38;2;222;220;184m") ; sys.stdout.write("\033[F")
        dicoInit["spinner"].text = "    • Create GFF DB"
        dicoInit["spinner"].start()
        fn = gffutils.example_filename(dicoInit["gff"])
        gffutils.create_db(fn, dbfn=dicoInit["gff"].replace(".gff",".db"), force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
        dicoInit["spinner"].stop()
        printcolor("    • GFF DB created","0","222;220;184",dicoInit["color"])
    else: printcolor("    • GFF DB found","0","222;220;184",dicoInit["color"])
    printcolor(" <> "+dicoInit["gff"].replace(".gff",".db")+"\n","0","117;116;97",dicoInit["color"])


#**************************************#
#       LOAD gffutils Database       #
#**************************************#
def load_gff_db(dicoInit):
    print("\x1b[0;38;2;222;220;184m") ; sys.stdout.write("\033[F")
    dicoInit["spinner"].text = "    • Load GFF DB"
    dicoInit["spinner"].start()
    dicoInit["db_gff"] = gffutils.FeatureDB(dicoInit["gff"].replace(".gff",".db"), keep_order=True)
    dicoInit["spinner"].stop()


#**************************************#
#        RETRIEVE genes features       #
#**************************************#
def retrieve_transcripts(dicoInit):
    # print("\x1b[0;38;2;255;187;108m") ; sys.stdout.write("\033[F")
    # dicoInit["spinner"].text = "    • Create GFF DB"
    # dicoInit["spinner"].start()

    #***** FInd target overlapping genes *****#
    if dicoInit['inputType']=="l":
        # Find intersecation between gff and target bed
        bed = BedTool(dicoInit["bed"])
        gff = BedTool(dicoInit["gff"])
        # Search gff exons intresection
        for intersect_elem in gff+bed:
            if intersect_elem.fields[2]=="exon":
                # retrieve correspunding transcript
                for rna in dicoInit["db_gff"].parents(dicoInit["db_gff"][intersect_elem.attrs["ID"]], featuretype='mRNA', order_by='start'):
                    gene = list(dicoInit["db_gff"].parents(rna, featuretype='gene', order_by='start'))[0]
                    print(gene)
        cleanup(remove_all=True) # delete created temp file





    # dicoInit["spinner"].stop()