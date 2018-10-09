# -*- coding: utf-8 -*-
# dooguypapua
# CHU Angers / PBMM
from tqdm import *
import subprocess
import sys
import threading
import time
from CovPlotter_display import *


#**************************************#
#          SIMPLE SUBPROCESS           #
#**************************************#
def launch_subprocess_thread(dicoThread,thread_name): # dicoThread[thread_name] = {"cmd":cmd, "returnstatut":None, "returnlines":[], "returnlog":""}
    process = subprocess.Popen(dicoThread[thread_name]["cmd"],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,shell=True)
    process.wait()
    for line in process.stdout.readlines():
        dicoThread[thread_name]["returnlines"].append(str(line.decode('utf-8').strip()))
    dicoThread[thread_name]["returnstatut"] = process.returncode


#**************************************#
#         MULTIPLE CLI THREADS         #
#**************************************#
def launch_threads(dicoInit,dicoThread,job_prefix,target_fct,thread_ratio):
    lstThread = []
    lstThreadcopy = []
    PrevNbFinishThread = 0
    # CLI Threads
    if target_fct==None:
        for thread_name in dicoThread.keys(): lstThread.append(threading.Thread(target=launch_subprocess_thread, args=(dicoThread,thread_name), name=job_prefix+"_"+thread_name))
    # Functions Threads
    else:
        for thread_name in dicoThread.keys(): lstThread.append(threading.Thread(target=target_fct, args=(dicoInit,dicoThread,thread_name), name=job_prefix+"_"+thread_name))
    # Tasks variables
    NbTasks = len(lstThread)
    lstThreadcopy = lstThread.copy()
    # Launch threads
    if NbTasks>0:
        print("\x1b[0;38;2;222;220;184m")
        with tqdm(total=NbTasks,ncols=35,position=0,leave=False,bar_format="    {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}") as pbar:
            # Continue while some tasks are availables
            while (len(lstThread)!=0 or job_prefix in str(threading.enumerate())):
                nbThreadActive = str(threading.enumerate()).count(job_prefix)
                if nbThreadActive<dicoInit['nt']*thread_ratio:
                    try: thread = lstThread.pop() ; thread.start() ; nbThreadActive+=1
                    except: pass
                NbFinishThread = NbTasks-len(lstThread)-nbThreadActive
                if PrevNbFinishThread!=NbFinishThread:
                    pbar.update(NbFinishThread-PrevNbFinishThread)
                    PrevNbFinishThread = NbFinishThread
                # check error
                for thread_name in dicoThread.keys():
                    if dicoThread[thread_name]["returnstatut"]!=0: break
                time.sleep(1)
            pbar.close()
        sys.stdout.write("\033[F")
    # Check dicoThread
    lst_error = []
    for name in dicoThread.keys():
        if dicoThread[name]["returnstatut"]!=0: lst_error.append(name)
    if len(lst_error)>0:
        sys.stdout.write('\x1b[1K')
        printcolor("\n  ERROR:","0","212;0;0",dicoInit['color'])
        printcolor(" "+lst_error[0]+"\n","0","212;0;0",dicoInit['color'])
        for error in lst_error[1:]: printcolor("         "+error+"\n","0","212;0;0",dicoInit['color'])