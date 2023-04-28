import glob
import subprocess
import os
import os.path as op
import sys

def runcmdhold(cmd, verbose = True):
    """
    cmd is a string
    """
    process = subprocess.Popen(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        text = True,
        shell = True
    )
    # std_out, std_err = process.communicate()
    # if verbose:
    #     print(cmd,std_out.strip(), std_err, end="")
    return process

def runcmd(cmd,ci, verbose = True):
    """
    cmd is a string
    """
    # out = open("run_many_mss_sim_output({}).out".format(ci),'w')
    err = open("run_many_mss_sim_errors({}).out".format(ci),'w')
    process = subprocess.Popen(
        cmd,
        # stdout = out,
        stderr = err,
        stdout = subprocess.PIPE,
        # stderr = subprocess.PIPE,
        text = True,
        shell = True
    )
    # std_out, std_err = process.communicate()
    # if verbose:
    #     print(cmd,std_out.strip(), std_err, end="")
    return process

def getjobslist(fn):
    lls = open(fn,"r").readlines()
    newlls = []
    for ls in lls:
         if len(ls) > 2:
              newlls.append(ls.strip())
    return newlls

if len(sys.argv) != 3:
     print("run_many_mss_sim.py usage:\n\t\ttwo arguments\n\t\t\t1. the path to the file of mss_sim.py commandlines\n\t\t\t2. the number of jobs to run at a time")
     exit()
jobsfile = sys.argv[1]
jobsatatime = int(sys.argv[2])
jobstrs = getjobslist(jobsfile)
# jobstrs = ["a" for i in range(100)]
# jobsatatime = 20
ji = 0
while ji < len(jobstrs):
    procs_list = [runcmd(cmd,ji + ci) for ci,cmd in enumerate(jobstrs[ji:ji+jobsatatime])]
    for proc in procs_list:
	    proc.wait()
    ji += jobsatatime
    print("done ",ji)
