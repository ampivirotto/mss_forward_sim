##run in linux
import glob
import subprocess
import os
import random
import string
import os.path as op
import platform
import sys

def run_job(command_list):
    """Run a job using subprocess with the given command list.

    Args:
        command_list (list): A list of strings that make up the command.

    Returns:
        None
    """
    # Redirect the standard output to subprocess.DEVNULL to suppress output
    subprocess.run(command_list, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


##raxml-ng --msa  temp.fa  --msa-format FASTA  --prefix  mytest --model GTR+G --threads 2 --seed 2

fainpath = sys.argv[1]
faoutpath = sys.argv[2]
#make the outpath if it does not exist
curdir = os.getcwd()
dirs = op.split(faoutpath)
for d in dirs:
    if op.exists(d) == False:
        os.mkdir(d)
    os.chdir(d)
os.chdir(curdir)
problemfn = "raxml_problems.txt"
fafns = glob.glob(op.join(fainpath,"*.fa*"))
renameFiles = {}
for fafn in fafns:
    try:
        #import pdb; pdb.set_trace()
        fafn_just_the_name = fafn[fafn.rfind("/")+1:]
        print(fafn_just_the_name)
        cmdlist = ["raxml-ng","--msa",fafn,"--msa-format","FASTA","--prefix","mrout","--model","GTR+G","--threads","2","--seed","123","--redo"]
        # subprocess.call(cmdlist)#,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        run_job(cmdlist)
        treestr = open("mrout.raxml.bestTree","r").readline().strip()
        #pdb.set_trace()
        faflines = open(fafn,'r').readlines()
        faflines.append("\n")
        faflines.append(treestr + "\n")
        try:
            newName = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
            newPath = op.join(faoutpath,newName + '.fa')
            while os.path.isfile(newPath):
                newName = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
                newPath = op.join(faoutpath,newName + '.fa')
            faf = open(newPath,'w')
        except:
            import pdb; pdb.set_trace()
        faf.writelines(faflines)
        faf.close()
        renameFiles[fafn_just_the_name] = [fafn_just_the_name, newName]
        pfiles = glob.glob("mrout*")
        for pf in pfiles:
            os.remove(pf)
    except:
        print("failed")
        fp = open(problemfn,"a")
        fp.write("{}\n".format(fafn))
        fp.close()

with open('{}/anon_match_{}.txt'.format(fainpath, faoutpath.split('/')[-1:][0]), 'w') as o:
    for key in renameFiles.keys():
        o.write("{}\t{}\n".format(renameFiles[key][0], renameFiles[key][1]))



