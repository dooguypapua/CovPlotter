# -*- coding: utf-8 -*-
# dooguypapua
# CHU Angers / PBMM
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from PIL import Image, ImageChops
from CovPlotter_bedgff import *
import gc


#**************************************#
#           CoverView PLOTTER          #
#**************************************#
def coverview_plot(dicoInit,dicoCov):
    print("\x1b[0;38;2;222;220;184m") ; sys.stdout.write("\033[F")
    dicoInit["spinner"].text = "• CovPlotter Threads"
    dicoInit["spinner"].start()
    dicoThread = {}
    # To avoid memory warning
    plt.rcParams.update({'figure.max_open_warning': 0})
    #***** THREAD PER GENES *****#
    for gene_id in dicoInit['lst_genes']:
        gene_name = dicoInit["db_gff"][gene_id].attributes["Name"][0]
        # subplots definition based on region without any exon for all transcripts
        dico_nonEmpty_interval,max_interval_length,max_exonic_size = gene_exonic_regions(dicoInit,dicoCov,gene_id,gene_name)
        # matplotlib multithread error => must init fig separately
        dico_plt = {}
        fig_width = max(int((max_exonic_size/1000)),5)
        for bam_num in dicoInit['dicoBam']: dico_plt[bam_num] = plt.figure(gene_name+"_"+str(bam_num),figsize=(fig_width,10))
        # dicoThread
        dicoThread[gene_name] = { "dico_plt":dico_plt, "dicoCov":dicoCov[gene_name], "dico_nonEmpty_interval":dico_nonEmpty_interval, "max_interval_length":max_interval_length, "gene_start":dicoInit["db_gff"][gene_id].start, "max_exonic_size":max_exonic_size, "returnstatut":-1}
        if len(dicoThread)>10: break
    dicoInit["spinner"].stop()
    printcolor("  • CovPlotter Threads\n","0","222;220;184",dicoInit["color"])
    printcolor("  • CovPlotter Plots","0","222;220;184",dicoInit["color"])
    launch_threads(dicoInit,dicoThread,"coverview gene plot",coverview_gene_plot,5)
    printcolor("  • CovPlotter Plots\n","0","222;220;184",dicoInit["color"])


def coverview_gene_plot(dicoInit,dicoThread,gene_name):
    try:
        #***** CREATE ONE PLOT GENE PER SAMPLE *****#
        for bam_num in dicoInit['dicoBam']:
            # plot out path
            path_plot_out = dicoInit["tmp"]+"/"+gene_name+"_"+str(bam_num)+".png"
            # Create gene figure
            fig = dicoThread[gene_name]["dico_plt"][bam_num]
            # Cpt for row and column subplot position
            plot_pos_start = 0
            cpt = 1
            row = 0
            dico_axis = {}
            dico_subplot_col_width = {}
            dicoCov = dicoThread[gene_name]["dicoCov"]
            #***** CREATE A SUBPLOT *****# (for each interval and each transcript)
            # Browse transcript
            for transcript_id in dicoCov:
                # Add a left empty subplot for display transcript_id
                ax = fig.add_subplot(len(dicoCov),len(dicoThread[gene_name]["dico_nonEmpty_interval"])+1,cpt)
                ax.set_title(transcript_id, size=12)
                plot_pos_start = max(plot_pos_start,ax.get_position().x1)
                ax.axis('off')
                dico_axis[row] = {}
                cpt+=1
                col = 0
                # Browse non empty interval position
                for interval in dicoThread[gene_name]["dico_nonEmpty_interval"]:
                    interval_start = int(interval.split("#")[0])
                    interval_end = int(interval.split("#")[1])
                    dico_axis[row][col] = fig.add_subplot(len(dicoCov),len(dicoThread[gene_name]["dico_nonEmpty_interval"])+1,cpt)
                    rectangles = {}
                    if transcript_id in dicoThread[gene_name]["dico_nonEmpty_interval"][interval]:
                        # Browse exons                       
                        for exon in dicoThread[gene_name]["dico_nonEmpty_interval"][interval][transcript_id]:
                            exon_start = dicoCov[transcript_id][exon]['start']
                            exon_end = dicoCov[transcript_id][exon]['end']
                            min_cov = dicoCov[transcript_id][exon]['dico_sample'][bam_num]
                            # Exon rectangle color
                            if min_cov<dicoInit['lcov']: color = ["#FA5858","#BB6060"]
                            elif min_cov<dicoInit['hcov']: color = ["#FAAC58","#BB9060"]
                            else: color = ["#58FA58","#60BB60"]
                            # Create rectangle only if we are in the correct interval and then subplot
                            if exon.__contains__("UTR"): rectangles[exon] = mpatch.Rectangle((exon_start-dicoThread[gene_name]["gene_start"],0.25), exon_end-exon_start, 0.5, color=color[1])
                            else: rectangles[exon] = mpatch.Rectangle((exon_start-dicoThread[gene_name]["gene_start"],0), exon_end-exon_start, 1, color=color[0])
                    # Add rectangles to subplot
                    for r in rectangles:
                        dico_axis[row][col].add_artist(rectangles[r])
                        rx, ry = rectangles[r].get_xy()
                        cx = rx + rectangles[r].get_width()/2.0
                        cy = ry + rectangles[r].get_height()/2.0
                        name = r.replace("E0","").replace("E","").replace("_UTR5","").replace("_UTR3","")
                        dico_axis[row][col].annotate(name, (cx, cy), color='black', weight='bold',fontsize=6, ha='center', va='center')
                    # set axis limit
                    dico_axis[row][col].set_ylim((0, 1))
                    dico_axis[row][col].set_xlim((interval_start, interval_end))
                    dico_axis[row][col].axis('off')
                    cpt+=1
                    col+=1
                row+=1

            #***** RESIZE SUBPLOT COLUMN *****# (for respect exon length)
            # Find real plot size withoptu margin, axis and transcript name
            plot_pos_end = dico_axis[row-1][col-1].get_position().x1
            plot_size = plot_pos_end-plot_pos_start-(0.005*len(dicoThread[gene_name]["dico_nonEmpty_interval"]))
            # For each subplot row
            for row in dico_axis:
                pos_x = dico_axis[row][0].get_position().x0
                previous_interval_end = 0
                # for each subplot column
                for col in dico_axis[row]:
                    pos1 = dico_axis[row][col].get_position()
                    interval = list(dicoThread[gene_name]["dico_nonEmpty_interval"].keys())[col]
                    interval_start = int(interval.split("#")[1])
                    interval_end = int(interval.split("#")[0])
                    length_interval = interval_start-interval_end
                    interval_width = (length_interval*plot_size)/dicoThread[gene_name]["max_exonic_size"]
                    pos2 = [pos_x, pos1.y0, interval_width, pos1.height]
                    dico_axis[row][col].set_position(pos2)
                    # Check if adjacent interval in the case of UTR
                    if previous_interval_end==interval_end or previous_interval_end==interval_start: pos_x+=interval_width
                    else: pos_x+=interval_width+0.005
                    previous_interval_end = interval_end

            #***** SAVE FILE *****#
            fig.savefig(path_plot_out, dpi=100)
            # Close plot and free memory object
            fig.clf()
            plt.close(gene_name+"_"+str(bam_num))
            gc.collect()

            #***** TRIM BLANK *****#
            im = Image.open(path_plot_out)
            im = trim(im)
            im.save(path_plot_out)
            im.close()

        #***** PLOTS to HTML to PDF *****#
        cover_plots_to_html(dicoInit,gene_name)
        cmd_wkhtmltopdf = "wkhtmltopdf -q -O Landscape "+dicoInit["tmp"]+"/"+gene_name+".html "+dicoInit["out"]+"/"+gene_name+".pdf 2>/dev/null"
        os.system(cmd_wkhtmltopdf)
        dicoThread[gene_name]["returnstatut"] = 0
    except: dicoThread[gene_name]["returnstatut"] = 1



def trim(im):
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)


#**************************************#
#          COVER PLOTS TO HTML         #
#**************************************#
def cover_plots_to_html(dicoInit,gene_name):
    lst_bam_num = list(dicoInit['dicoBam'].keys())
    # determine number of image per row
    im = Image.open(dicoInit["tmp"]+"/"+gene_name+"_"+str(lst_bam_num[0])+".png")
    width, height = im.size
    im.close()
    # HTML style
    path_html_out = dicoInit["tmp"]+"/"+gene_name+".html"
    HTML = open(path_html_out,'w')
    HTML.write("<style>\n")
    HTML.write("\t.plottable {border-collapse: collapse;font-weight:bold; font-size:90%;}\n")
    HTML.write("\t.plotrow {border: 1px solid black;padding: 10px;text-align: center;}\n")
    HTML.write("\t.pbi_avoid {page-break-inside: avoid !important;}\n")
    HTML.write("\t.title {word-wrap: break-word;width:100%; max-width: "+str(int(width/2))+"px; font-size:10px;}\n")
    HTML.write("</style>\n")
    HTML.write("<table>\n")
    # Group barcodes list per 3
    split_num = max(int(3500/width),1)
    for sample_grp in [lst_bam_num[n:n+split_num] for n in range(0, len(lst_bam_num), split_num)]:
        HTML.write("\t<tr>\n")
        for i in range(len(sample_grp)):
            HTML.write("\t\t<td>\n")
            HTML.write("\t\t\t<table class='plottable pbi_avoid'>\n")
            HTML.write("\t\t\t\t<tr><td class='plotrow'>"+dicoInit['dicoBam'][sample_grp[i]]+"<br></td></tr>\n")
            HTML.write("\t\t\t\t<tr><td class='plotrow'><img align='center' height='50%' src='"+dicoInit["tmp"]+"/"+gene_name+"_"+str(sample_grp[i])+".png'/></td></tr>\n")
            HTML.write("\t\t\t</table>\n")
            HTML.write("\t\t</td>\n")
        HTML.write("\t</tr>\n")
    HTML.write("</table>\n")
    HTML.close()