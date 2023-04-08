import glob
import subprocess
import os
import os.path as op
import sys

def runcmd(cmd, verbose = True):
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
    #     print(std_out.strip(), std_err, end="")
    return process

def getjobslist(fn):
    lls = open(fn,"r").readlines()
    newlls = []
    for ls in lls:
         if len(ls) > 2:
              newlls.append(ls.strip())
    return newlls

jobsfile = sys.argv[1]
jobsatatime = int(sys.argv[2])
jobstrs = getjobslist(jobsfile)
# jobstrs = ["a" for i in range(100)]
# jobsatatime = 20
ji = 0
while ji < len(jobstrs):
    procs_list = [runcmd(cmd) for cmd in jobstrs[ji:ji+jobsatatime]]
    for proc in procs_list:
	    proc.wait()
    ji += jobsatatime
    print("done ",ji)
