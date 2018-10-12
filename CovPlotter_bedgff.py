# -*- coding: utf-8 -*-
# dooguypapua
import sys,os,shutil,re,glob
import subprocess
import gffutils
import urllib.request
from pybedtools import BedTool,set_tempdir
from itertools import groupby
from CovPlotter_init import *
from CovPlotter_display import *
from CovPlotter_process import *


#**************************************#
#            DOWNLOAD GFF              #
#**************************************#
def check_gff(dicoInit):
    server = "http://163.172.45.124/share/db/"
    # Create db folder
    os.makedirs(os.getcwd()+"/db", exist_ok=True)
    # Check online files
    listfile = urllib.request.urlopen(server)
    for line in listfile.readlines():
        search_file = re.search(">("+dicoInit["db"]+"_"+dicoInit["ref"]+"_.+\.gff.gz)</a></td>",str(line))
        if search_file:
            urlpath = search_file.group(1)
            urlversion = urlpath.split("_")[2].replace(".gff.gz","")
    # Check local files
    dl = False
    glob_gff = glob.glob(os.getcwd()+"/db/"+dicoInit["db"]+"_"+dicoInit["ref"]+"*.gff")
    glob_gff_db = glob.glob(os.getcwd()+"/db/"+dicoInit["db"]+"_"+dicoInit["ref"]+"*.db")
    # Missing one file
    if len(glob_gff)!=len(glob_gff_db):
        try: os.remove(glob_gff[0])
        except: pass
        try: os.remove(glob_gff_db[0])
        except: pass
        dl = True
    # Missing two files
    elif len(glob_gff)==len(glob_gff_db)==1:
        dicoInit["gff"] = glob_gff[0]
        dicoInit["gff_db"] = glob_gff_db[0]
        version = os.path.basename(dicoInit["gff"]).split("_")[2].replace(".gff","")
        # Update version is available
        if int(urlversion.replace("v",""))>int(version.replace("v","")):
            printcolor("  Found "+dicoInit["db"]+" "+dicoInit["ref"]+" "+version+" but "+urlversion+" is available.\n  Do you want to upgrade (y/n)? ","0","222;220;184",dicoInit["color"])
            choice = input().lower()
            while choice!="y" and choice!="n":
                printcolor("  Do you want to upgrade (y/n)? ","0","222;220;184",dicoInit["color"])
                choice = input().lower()
            # Delete older version
            if choice=="y":
                dl = True
                os.remove(dicoInit["gff"]) ; os.remove(dicoInit["gff_db"])
                sys.stdout.write("\033[F") ; sys.stdout.write("\x1b[2K") ; sys.stdout.write("\033[F") ; sys.stdout.write("\x1b[2K")
    else: dl = True
    # Download files
    if dl==True:
        printcolor("  • Download GFF/DB files\n","0","222;220;184",dicoInit["color"])
        print("\x1b[0;38;2;222;220;184m") ; sys.stdout.write("\033[F")
        dicoInit["gff"] = os.getcwd()+"/db/"+dicoInit["db"]+"_"+dicoInit["ref"]+"_"+urlversion+".gff"
        dicoInit["gff_db"] = os.getcwd()+"/db/"+dicoInit["db"]+"_"+dicoInit["ref"]+"_"+urlversion+".db"
        urllib.request.urlretrieve(server+"/"+urlpath,dicoInit['tmp']+"/"+urlpath,reporthook)
        urllib.request.urlretrieve(server+"/"+urlpath.replace("gff","db"),dicoInit['tmp']+"/"+urlpath.replace("gff","db"),reporthook)        
        # Extract gzip
        GZ1 = gzip.open(dicoInit['tmp']+"/"+urlpath, 'rb')
        GFF = open(dicoInit["gff"], 'wb')
        GFF.write(GZ1.read())
        GZ1.close() ; GFF.close()
        GZ2 = gzip.open(dicoInit['tmp']+"/"+urlpath.replace("gff","db"), 'rb')
        DB = open(dicoInit["gff_db"], 'wb')
        DB.write(GZ2.read())
        GZ2.close() ; DB.close()        
        sys.stdout.write("\x1b[2K") ; sys.stdout.write("\033[F")
        printcolor("  • GFF/DB downloaded","0","222;220;184",dicoInit["color"])
        printcolor(" <> "+dicoInit["gff"].replace(".gff",".(gff/db)")+"\n","0","117;116;97",dicoInit["color"])
   

#**************************************#
#       LOAD gffutils Database       #
#**************************************#
def load_gff_db(dicoInit):
    dicoInit["db_gff"] = gffutils.FeatureDB(dicoInit["gff"].replace(".gff",".db"), keep_order=True)
    printcolor("  • GFF/DB loaded ","0","222;220;184",dicoInit["color"])
    printcolor(" <> "+dicoInit["gff_db"]+"\n","0","117;116;97",dicoInit["color"])



#**************************************#
#        RETRIEVE genes features       #
#**************************************#
def retrieve_transcripts(dicoInit):
    #***** Find target overlapping genes *****#
    if dicoInit['inputType']=="l":
        print("\x1b[0;38;2;222;220;184m") ; sys.stdout.write("\033[F")
        dicoInit["spinner"].text = "• Extract target genes"
        dicoInit["spinner"].start()
        # Find intersecation between gff and target bed
        bed = BedTool(dicoInit["bed"])
        gff = BedTool(dicoInit["gff"])
        set_rna = set()
        set_gene = set()
        # Search gff exons intresection
        for intersect_elem in gff+bed:
            if intersect_elem.fields[2]=="exon":
                # retrieve correspunding transcript
                for rna in dicoInit["db_gff"].parents(dicoInit["db_gff"][intersect_elem.attrs["ID"]], featuretype='mRNA', order_by='start'):
                    gene = list(dicoInit["db_gff"].parents(rna, featuretype='gene', order_by='start'))[0]
                    set_rna.add(rna)
                    set_gene.add(gene)
        dicoInit["spinner"].stop()
        dicoInit['lst_genes'] = list(set_gene)
        printcolor("  • Target genes  ","0","222;220;184",dicoInit["color"])
        printcolor(" <> "+str(len(dicoInit['lst_genes']))+" genes for "+str(len(set_rna))+" transcripts\n","0","117;116;97",dicoInit["color"])


#**************************************#
#   CREATE target genes features BED   #
#**************************************#
def create_detailed_bed(dicoInit):
    printcolor("  • CovPlotter BED\n","0","222;220;184",dicoInit["color"])
    path_target_genes_bed = dicoInit["tmp"]+"/target_genes.bed"
    BED = open(path_target_genes_bed,'w')
    for gene in dicoInit['lst_genes']:
        gene_name = gene.attributes["Name"][0]
        for rna in dicoInit["db_gff"].children(gene, featuretype='mRNA', order_by='start'): # .attributes => "ID","Parent","Dbxref","Name","gbkey","gene","product","transcript_id"
            # Search CDS start and stop
            cds_start = list(dicoInit["db_gff"].children(rna, featuretype='CDS', order_by='start'))[0].start
            cds_end = list(dicoInit["db_gff"].children(rna, featuretype='CDS', order_by='start'))[-1].end
            # Browse exon
            if rna.strand=="+": exon_cpt = 1
            else: exon_cpt = len(list(dicoInit["db_gff"].children(rna, featuretype='exon', order_by='start')))
            for exon in dicoInit["db_gff"].children(rna, featuretype='exon', order_by='start'):
                # Split exon including CDS
                if exon.start<cds_start:
                    if exon.end>cds_start:
                        if cds_end<exon.end: # case of one exon including all the cds
                            if exon.strand=="+":
                                BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(cds_start-1)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR5\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                                BED.write(exon.chrom+"\t"+str(cds_end)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR3\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                            else:
                                BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(cds_start-1)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR3\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                                BED.write(exon.chrom+"\t"+str(cds_end)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR5\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                            BED.write(exon.chrom+"\t"+str(cds_start-1)+"\t"+str(cds_end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                        else: 
                            if exon.strand=="+": BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(cds_start-1)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR5\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                            else: BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(cds_start-1)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR3\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                            BED.write(exon.chrom+"\t"+str(cds_start-1)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                    elif exon.strand=="+": BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR5\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                    else: BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR3\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                elif exon.end>cds_end:
                    if exon.start<cds_end:
                        if exon.strand=="+": BED.write(exon.chrom+"\t"+str(cds_end)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR3\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                        else: BED.write(exon.chrom+"\t"+str(cds_end)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR5\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                        BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(cds_end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                    elif exon.strand=="+":  BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR3\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                    else: BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"_UTR5\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                else:
                    BED.write(exon.chrom+"\t"+str(exon.start-1)+"\t"+str(exon.end)+"\t"+gene_name+"#"+rna.attributes["transcript_id"][0]+"#E"+str(exon_cpt).zfill(2)+"\t"+str(exon.score)+"\t"+str(exon.strand)+"\n")
                if rna.strand=="+": exon_cpt+=1
                else: exon_cpt-=1
    BED.close()


#**************************************#
#     GENE REGION WITHOUT ANY EXON     #
#**************************************#
def gene_exonic_regions(dicoInit,dicoCov,gene,gene_name):
    dico_pos_grp_type = {}
    dico_exon_pos = {}
    exonic_size = 0
    max_exonic_size = 0
    # Initialize dico with empty group type for all positions
    for pos in range(0,gene.end-gene.start+1,1):
        dico_pos_grp_type[pos] = "0"
        dico_exon_pos[pos] = {}
    # Browse transcripts to find position with at least one exon
    for transcript_id in dicoCov[gene_name]:
        for exon in dicoCov[gene_name][transcript_id]:
            relative_exon_start = dicoCov[gene_name][transcript_id][exon]["start"]-gene.start
            relative_exon_end = dicoCov[gene_name][transcript_id][exon]["end"]-gene.start
            exonic_size+=(relative_exon_end-relative_exon_start)
            if exon.__contains__("UTR"): grp_type = "2"
            else: grp_type = "1"
            for pos in range(relative_exon_start,relative_exon_end+1,1):
                dico_pos_grp_type[pos] = grp_type
                try: dico_exon_pos[pos][transcript_id].append(exon)
                except: dico_exon_pos[pos][transcript_id] = [exon]
        if exonic_size>max_exonic_size: max_exonic_size = exonic_size
    # Convert to String and group per type
    str_cover = "".join(list(dico_pos_grp_type.values()))
    groups = groupby(str_cover)
    lst_grp = [(label, sum(1 for _ in group)) for label, group in groups]
    # Define not empty interval
    dico_nonEmpty_interval = {}
    max_interval_length = 0
    current_pos = 0
    for grp in lst_grp: # grp = [grp_type,grp_length]
        if grp[0]!="0":
            dico_nonEmpty_interval[str(current_pos)+"#"+str(current_pos+grp[1])] = dico_exon_pos[current_pos]
            if grp[1]>max_interval_length: max_interval_length = grp[1]
        current_pos+=grp[1]
    return dico_nonEmpty_interval,max_interval_length,max_exonic_size



#**************************************#
#       LAUNCH bedtools coverage       #
#**************************************#
def pybedtoolcoverage(dicoInit,dicoThread,thread_name):
    bed = BedTool(dicoThread[thread_name]['bed'])
    bam = BedTool(dicoThread[thread_name]['bam'])
    cov = bed.coverage(bam,d=True)
    shutil.move(cov.fn,dicoInit['tmp']+"/"+str(dicoThread[thread_name]['bam_num'])+".cov")

def launch_coverage(dicoInit):
    printcolor("  • Compute Depth","0","222;220;184",dicoInit["color"])
    dicoThread = {}
    set_tempdir(dicoInit['tmp'])
    for bam_num in dicoInit['dicoBam'].keys():
        dicoThread["coverage "+dicoInit['dicoBam'][bam_num]] = {"bed":dicoInit["tmp"]+"/target_genes.bed", "bam":dicoInit['dicoBam'][bam_num], "bam_num":bam_num, "returnstatut":None, "returnlines":[]}
    launch_threads(dicoInit,dicoThread,"pybedtoolcoverage",pybedtoolcoverage,1)

def parse_coverage(dicoInit):
    print("\x1b[0;38;2;222;220;184m") ; sys.stdout.write("\033[F")
    dicoInit["spinner"].text = "• CovPlotter Depth"
    dicoInit["spinner"].start()
    dicoCov = {}
    for cov_out in glob.glob(dicoInit["tmp"]+"/*.cov"):
        sample = int(os.path.basename(cov_out).replace(".cov",""))
        COV = open(cov_out,'r')
        lst_lines = COV.read().split("\n")
        COV.close()
        for line in lst_lines:
            if line!="":
                # Read
                split_line = line.split("\t")
                split_feature = split_line[3].split("#")
                gene = split_feature[0]
                transcript = split_feature[1]
                exon = split_feature[2]
                depth = int(split_line[7])
                # Add to dico_cov
                if not gene in dicoCov: dicoCov[gene] = {}
                if not transcript in dicoCov[gene]: dicoCov[gene][transcript] = {}
                if not exon in dicoCov[gene][transcript]: dicoCov[gene][transcript][exon] = { 'start':int(split_line[1])+1, 'end':int(split_line[2]), 'dico_sample':{} }
                if not sample in dicoCov[gene][transcript][exon]['dico_sample']: dicoCov[gene][transcript][exon]['dico_sample'][sample] = [depth]
                dicoCov[gene][transcript][exon]['dico_sample'][sample].append(depth)
    # Compute min depth per region
    for gene in dicoCov:
        for transcript in dicoCov[gene]:
            for exon in dicoCov[gene][transcript]:
                for sample in dicoCov[gene][transcript][exon]['dico_sample']: dicoCov[gene][transcript][exon]['dico_sample'][sample] = min(dicoCov[gene][transcript][exon]['dico_sample'][sample])
    # Return dicoCov
    dicoInit["spinner"].stop()
    return dicoCov
    printcolor("  • CovPlotter Depth\n","0","222;220;184",dicoInit["color"])
