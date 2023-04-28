"""
based on  make_file_of_mss_sim_commands.py
but works over a range of parameter values

"""
import os
import os.path as op
import argparse
import sys
import numpy as np
import math


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
    parser.add_argument("-e", help="random integer number seed for picking alignments."
                        "If value is negative,  pick just one gene. ",dest="ranseed",type=int)
    parser.add_argument("-F", help="directory path for output fasta file (default is same as for results and log files)",dest="fdir",type = str)
    parser.add_argument("-k", help="Number of species (4,5 or 11)",dest="numSpecies",default=4,type=int)
    parser.add_argument("-L", help="Length of sequence (# amino acids)", dest="aalength",default=300,type=int)
    parser.add_argument("-m", help="Model file path",dest="mssmodelfilename",required = True,type = str)
    parser.add_argument("-N", help="Population size (diploid)."
                        "If used twice, values are upper and lower bounds."
                        "If used three times, the third value should be 'log' to make it on a log scale", 
                        dest="popsize",default=[],action="append",type=str)
    parser.add_argument("-R", help="directory path for results and log files",dest="rdir",default = ".",type = str)
    parser.add_argument("-s", help="Synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))."
                        " If used twice, values are upper and lower bounds."
                        " If used three times, the third value should be 'log' to make it on a log scale", 
                        dest="SynSel_s",default=[],action="append",type=str)
    parser.add_argument("-y", help="Non-synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))."
                        "If used twice, values are upper and lower bounds."
                        "If used three times, the third value should be 'log' to make it on a log scale",  
                        dest="NonSyn_s",default=[],action='append',type=str)
    parser.add_argument("-u", help="expected number of neutral mutations per site, from base of tree."
                        "If used twice, values are upper and lower bounds." 
                        "If used three times, the third value should be 'log' to make it on a log scale",  
                        dest="mutationexpectation",default=[],action="append",type=str)
    parser.add_argument("-v", help="Number of values (intervals + 1) on each varying dimension",dest="numvalues",required=True,type=int)
    parser.add_argument("-c", help="Path to file of list of commands", dest="cmdfn", type=str,required=True)
    # parser.add_argument("-j", help="Number of jobs",dest="numjobs",required=True,type=int)
    return parser

def create_evenly_spaced_list(lower, upper, num_values,dolog = False):
    """
    Create a list of evenly spaced values between a lower and upper bound, including the bounds in the list.

    Args:
        lower (float): The lower bound.
        upper (float): The upper bound.
        num_values (int): The number of values in the list.

    Returns:
        list: A list of evenly spaced values.
    """
    if dolog:
            if lower <= 0:
                print("can't do log scale with lower bound <= 0")
                exit()
            ratio = upper / lower
            # Calculate the step size between consecutive values
            step = math.pow(ratio, 1/num_values)
            # Initialize the list of values with the minimum value
            values = [lower]

            # Calculate the remaining values
            for i in range(num_values):
                value = values[-1] * step
                values.append(value)
    else:
        step = (upper - lower) / (num_values - 1)  # Calculate the step size
        values = [lower + i * step for i in range(num_values)]  # Generate the list of evenly spaced values
    return values

def varyparamwork(args):
    numjobs = 1
    if len(args.SynSel_s) > 1:
        if len(args.SynSel_s)==3: # do log scale
            args.SynSel_s = list(map(float,args.SynSel_s[:2]))
            args.SynSel_s.sort()
            x = create_evenly_spaced_list(args.SynSel_s[0],args.SynSel_s[1],args.numvalues,dolog=True)
        else:
            args.SynSel_s = list(map(float,args.SynSel_s))
            args.SynSel_s.sort()
            x = create_evenly_spaced_list(args.SynSel_s[0],args.SynSel_s[1],args.numvalues)
        args.SynSel_s = list(map(round,x,[3]*len(x)))
        args.SynSel_s.sort()
        numjobs *= len( args.SynSel_s)
    else:
        assert len(args.SynSel_s)==1
        args.SynSel_s = [float(args.SynSel_s[0])]
    if len(args.NonSyn_s) > 1:
        if len(args.NonSyn_s)==3: # do log scale
            args.NonSyn_s = list(map(float,args.NonSyn_s[:2]))
            args.NonSyn_s.sort()
            x = create_evenly_spaced_list(args.NonSyn_s[0],args.NonSyn_s[1],args.numvalues,dolog=True)
        else:
            args.NonSyn_s = list(map(float,args.NonSyn_s))
            args.NonSyn_s.sort()
            x = create_evenly_spaced_list(args.NonSyn_s[0],args.NonSyn_s[1],args.numvalues)
        args.NonSyn_s = list(map(round,x,[3]*len(x)))
        numjobs *= len( args.NonSyn_s)
    else:
        assert len(args.NonSyn_s)==1
        args.NonSyn_s = [float(args.NonSyn_s[0])]
    if len(args.mutationexpectation) > 1:
        if len(args.mutationexpectation)==3: # do log scale
            args.mutationexpectation = list(map(float,args.mutationexpectation[:2]))
            args.mutationexpectation.sort()
            x = create_evenly_spaced_list(args.mutationexpectation[0],args.mutationexpectation[1],args.numvalues,dolog=True)
        else:
            args.mutationexpectation = list(map(float,args.mutationexpectation))
            args.mutationexpectation.sort()
            x = create_evenly_spaced_list(args.mutationexpectation[0],args.mutationexpectation[1],args.numvalues)
        args.mutationexpectation = list(map(round,x,[3]*len(x)))
        numjobs *= len( args.mutationexpectation)
    else:
        assert len(args.mutationexpectation)==1
        args.mutationexpectation = [float(args.mutationexpectation[0])]

    if len(args.popsize) > 1:
        if len(args.popsize)==3: # do log scale
            args.popsize = list(map(float,args.popsize[:2]))
            args.popsize.sort()
            x = create_evenly_spaced_list(args.popsize[0],args.popsize[1],args.numvalues,dolog=True)
        else:
            args.popsize = list(map(float,args.popsize))
            args.popsize.sort()
            x = create_evenly_spaced_list(args.popsize[0],args.popsize[1],args.numvalues)
        args.popsize = list(map(round,x))
        numjobs *= len( args.popsize)
    else:
        assert len(args.popsize)==1
        args.popsize = [int(args.popsize[0])]

    args.numjobs = numjobs
    return args       


def main(argv):
    parser = parseargs()
    if argv[-1] =='':
        argv = argv[0:-1]
    args = parser.parse_args(argv)
    if args.ranseed != None:
        just1gene = args.ranseed < 0
        args.ranseed = abs(args.ranseed)
    np.random.seed(args.ranseed)
    if args.numSpecies not in [4,5,11]:
        print ("error: -p (# of species) must be 4,5 or 11")
        exit()
    args=varyparamwork(args)
    listofrandomgenes = getlistofRandomGenes(args.bacaligndir,args.numjobs)
    if just1gene:
        gene = listofrandomgenes[0]
    seeds = set()
    while True:
        seeds.add(np.random.randint(100,100000000))
        if len(seeds) == args.numjobs:
            break
    seeds = list(seeds)

    cmd_constant = ["python mss_sim.py"]
    temp = ["-A",str(args.bacaligndir)]
    cmd_constant += temp
    
    temp = ["-L",str(args.aalength)]
    cmd_constant += temp    
    temp = ["-k",str(args.numSpecies)]
    cmd_constant += temp
    temp = ["-m",str(args.mssmodelfilename)]
    cmd_constant += temp
    temp = ["-R",str(args.rdir)]
    cmd_constant += temp    
    temp = ["-F",str(args.fdir)]
    cmd_constant += temp
    f = open(args.cmdfn,'w')
    i = 0
    for popsize in args.popsize:
        for  SynSel_s in args.SynSel_s:
            for NonSyn_s in args.NonSyn_s:
                for mutationexpectation in args.mutationexpectation:
                    cmd = cmd_constant.copy()
                    augname = "{}_N{}_s{}_y{}_u{}".format(args.basename,popsize,SynSel_s,NonSyn_s,mutationexpectation)
                    temp = ["-b",str(augname)]
                    cmd += temp
                    temp = ["-N",str(popsize)]
                    cmd += temp
                    temp = ["-u",str(mutationexpectation)]
                    cmd += temp
                    temp = ["-s",str(SynSel_s)]
                    cmd += temp
                    temp = ["-y",str(NonSyn_s)]
                    cmd += temp
                    if just1gene:
                        cmd += ["-e",str(seeds[i]),"-g",gene]
                    else:
                        cmd += ["-e",str(seeds[i]),"-g",listofrandomgenes[i]]
                    f.write("{}\n".format(" ".join(cmd)))
                    i += 1
    f.close()



if __name__ == "__main__":

    if len(sys.argv) < 2:
        #sys.stderr.write("No arguments. Use -h  or --help for help menu")
        main(['-h'])
    else:
        main(sys.argv[1:])
        