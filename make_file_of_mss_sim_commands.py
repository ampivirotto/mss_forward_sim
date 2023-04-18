import os
import os.path as op
import argparse
import sys
import numpy as np


def identifyRandomGene(alignment_location):
    """
    picks random gene from fasta alignment files
    """
    geneFiles = []
    for file in os.listdir(alignment_location):
        if file.endswith('fasta'):
            geneFiles.append(file)
    gene = np.random.choice(geneFiles)
    # print("Picked: " + gene)
    return gene

def getlistofRandomGenes(alignment_location,njobs):
    lrg = []
    while len(lrg) < njobs:
        gene = identifyRandomGene(alignment_location)
        if gene not in lrg:
            lrg.append(gene)
    return lrg

def parseargs():
    parser = argparse.ArgumentParser("python make_file_of_mss_sim_commands.py",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-A", help="path for folder containing bacterial alignment files",dest="bacaligndir",required=True,type = str)
    parser.add_argument("-b", help="base filename for putput files (name only, no directories)",dest="basename",required=True,type = str)
    parser.add_argument("-e", help="random number seed for picking alignments",dest="ranseed",type=int)
    parser.add_argument("-F", help="directory path for output fasta file (default is same as for results and log files)",dest="fdir",type = str)
    parser.add_argument("-k", help="Number of species (4,5 or 11)",dest="numSpecies",default=4,type=int)
    parser.add_argument("-L", help="Length of sequence (# amino acids)", dest="aalength",default=300,type=int)
    parser.add_argument("-m", help="Model file path",dest="mssmodelfilename",required = True,type = str)
    parser.add_argument("-N", help="Population size (diploid)",dest="popsize",default=10,type=int)
    parser.add_argument("-R", help="directory path for results and log files",dest="rdir",default = ".",type = str)
    parser.add_argument("-s", help="Synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))", dest="SynSel_s",default=2,type=float)
    parser.add_argument("-y", help="Non-synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))", dest="NonSyn_s",default=10,type=float)
    parser.add_argument("-u", help="expected number of neutral mutations per site, from base of tree", dest="mutationexpectation",default=0.5,type=float)
    parser.add_argument("-c", help="Path to file of list of commands", dest="cmdfn", type=str,required=True)
    parser.add_argument("-j", help="Number of jobs",dest="numjobs",required=True,type=int)


    return parser
def main(argv):
    parser = parseargs()
    if argv[-1] =='':
        argv = argv[0:-1]
    args = parser.parse_args(argv)
    if args.ranseed != None:
        np.random.seed(args.ranseed)
    if args.numSpecies not in [4,5,11]:
        print ("error: -p (# of species) must be 4,5 or 11")
        exit()
    listofrandomgenes = getlistofRandomGenes(args.bacaligndir,args.numjobs)
    cmd = ["python mss_sim.py"]
    temp = ["-A",str(args.bacaligndir)]
    cmd += temp
    temp = ["-b",str(args.basename)]
    cmd += temp
    temp = ["-N",str(args.popsize)]
    cmd += temp
    temp = ["-L",str(args.aalength)]
    cmd += temp
    temp = ["-u",str(args.mutationexpectation)]
    cmd += temp
    temp = ["-s",str(args.SynSel_s)]
    cmd += temp
    temp = ["-y",str(args.NonSyn_s)]
    cmd += temp
    temp = ["-k",str(args.numSpecies)]
    cmd += temp
    temp = ["-m",str(args.mssmodelfilename)]
    cmd += temp
    temp = ["-R",str(args.rdir)]
    cmd += temp
    temp = ["-F",str(args.fdir)]
    cmd += temp
    seeds = set()
    while True:
        seeds.add(np.random.randint(100,100000000))
        if len(seeds) == args.numjobs:
            break
    f = open(args.cmdfn,'w')
    for i,e in enumerate(seeds):
        newcmd = cmd.copy()
        newcmd += ["-e",str(e),"-g",listofrandomgenes[i]]
        f.write("{}\n".format(" ".join(newcmd)))
    f.close()



if __name__ == "__main__":

    if len(sys.argv) < 2:
        #sys.stderr.write("No arguments. Use -h  or --help for help menu")
        main(['-h'])
    else:
        main(sys.argv[1:])
        