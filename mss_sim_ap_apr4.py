"""
    forward simulation under an MSS model
    officially diploid,  but only in doubling of # of chromosomes
        throughout popsize refers to diploid # and popsize2 to haploid number
    assumes selection at the haplotype/chromosome level

    lifecycle:
        a chromosome is sampled at random from one of the possible parents, based on the parent fitnesses
        after it is sampled new mutations are added and its fitness is calculated


"""
import os
import os.path as op
import argparse
import numpy as np
import pandas as pd
import sys


def identifyRandomGene(alignment_location):
    """
    picks random gene from fasta alignment files
    """
    geneFiles = []
    for file in os.listdir(alignment_location):
        if file.endswith('fasta'):
            geneFiles.append(file)
    gene = np.random.choice(geneFiles)
    print("Picked: " + gene)
    return gene

def readInGeneFile(alignment_location):
    """
    given random gene selected, returns the alignment data for that gene
    """
    gene = identifyRandomGene(alignment_location)
    with open(op.join(alignment_location,gene)) as f:
        species = {}
        dna = ''
        spec = ''
        for line in f:
            if line.startswith('>'):
                assert len(dna) % 3 == 0
                species[spec] = dna
                spec = line.lstrip('>').rstrip('\n')
                dna = ''
            else:
                dna += line.strip('\n')
    return species, gene

def createCodonSequence(alignment_location):
    """
    given alignments for randomly selected gene, it turns all alignments into single strand of DNA excluding codonds with missing bps or stop codons
    """
    species, gene = readInGeneFile(alignment_location)
    allDNA = ''
    stopCodons = ['TAG', 'TAA', 'TGA']
    for x in species:
        dna = species[x]
        codon = ''
        if len(dna) > 0:
            for bp in dna:
                if len(codon) == 3:
                    if '-' in codon:
                        codon = ''
                        continue
                    if codon in stopCodons:
                        codon = ''
                        continue
                    if codon == 'ATG':
                        codon = ''
                        continue
                    allDNA += codon
                    codon = ''
                codon += bp
    return allDNA, gene

def writeCodonsToFile(randomgenebasepath,alignment_location):
    """
    writes out codon (in bp) to file
    """
    allDNA, genefilename = createCodonSequence(alignment_location)
    genename = genefilename[:genefilename.find('_')]
    #checkCodons(allDNA, gene)
    fname = randomgenebasepath + "_" + genename + "_allDNAstrand.txt"
    o = open(fname, 'w')
    #o.write('>' + gene.strip('.fasta') + '\n')
    o.write(allDNA)
    o.close()
    return fname, genename, allDNA

def getCodonProportions(dna):
    """
    calculate proportion of codons for a given strand of DNA, return dictionary
    """
    codons = {}
    codon = ''
    total = 0

    for bp in dna:
        if len(codon) == 3:
            if codon in codons.keys():
                codons[codon] += 1
            else:
                codons[codon] = 1
            codon = ''
            total += 1
        codon += bp

    for key in codons.keys():
        codons[key] = codons[key] / total

    return codons

def makeAncestor(allDNA, aalen):
    """
    take allDNA and take a random assortment of codons based on appearing in allDNA
    """
    props = getCodonProportions(allDNA)

    ancestor = np.random.choice(list(props.keys()), size = aalen, p = list(props.values()))

    return ''.join(ancestor)


def countCodons(dna):
    codonDict = {}
    count = 0
    codon = ''
    for bp in dna:
        codon += bp
        if count == 2:
            if codon in codonDict.keys():
                codonDict[codon] += 1
            else:
                codonDict[codon] = 1
            codon = ''
            count = 0
        else:
            count += 1
    return codonDict

def checkCodons(dna, gene):
    codonDict = countCodons(dna)

    bacterial = pd.read_csv('bacterial_data/bacterial_tbl.csv')
    row = bacterial[bacterial['geneName'] == gene.split("_")[0]]
    expected = {}
    #print(row)
    for c in codonDict.keys():
        expected[c] = list(row[c])[0] / list(row['numSpecies'])[0] / list(row['geneLength'])[0]

    for c in codonDict.keys():
        print(codonDict[c], expected[c])

def readModelFile(df):
    selected = df[df['GROUP'] == 'SELECTED']

    sdict = {}

    for aa in selected['AminoAcid'].unique():
        codonDF = selected[selected['AminoAcid'] == aa]

        for codon in codonDF['CODON']:
            if convertAAformat(aa) in sdict.keys():
                sdict[convertAAformat(aa)].append(codon)
            else:
                sdict[convertAAformat(aa)] = [codon]
    return sdict

def convertAAformat(aa):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLT': 'E', 'TYR': 'Y', 'MET': 'M'}

    if len(aa) == 3:
        return d[aa]
    elif len(aa) == 1:
        for dkey in d.keys():
            if d[dkey] == aa:
                return dkey

def codonInfo():

    codons ={   "I":["ATT", "ATC", "ATA"],
                "L":["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
                "V":["GTT", "GTC", "GTA", "GTG"],
                "F":["TTT", "TTC"],
                "M":["ATG"],
                "C":["TGT", "TGC"],
                "A":["GCT", "GCC", "GCA", "GCG"],
                "G":["GGT", "GGC", "GGA", "GGG"],
                "P":["CCT", "CCC", "CCA", "CCG"],
                "T":["ACT", "ACC", "ACA", "ACG"],
                "S":["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
                "Y":["TAT", "TAC"],
                "W":["TGG"],
                "Q":["CAA", "CAG"],
                "N":["AAT", "AAC"],
                "H":["CAT", "CAC"],
                "E":["GAA", "GAG"],
                "D":["GAT", "GAC"],
                "K":["AAA", "AAG"],
                "R":["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
                "STOP":["TAA", "TAG", "TGA"]}

    aalist = list(codons.keys())
    aalist.remove('STOP')
    codonlist = []

    revCodons  = {}

    for aa in aalist:
        for cd in codons[aa]:
            codonlist.append(cd)
            revCodons[cd] = aa

    return codons, aalist, codonlist, revCodons

def createSelectedDictionary(args):
    """
    if some structure needs to be built that represents codon fitnesses efficiently,  this is the place for it
    """
    selectedDict = {}
    mutDict = {}

    codons, aalist, codonlist, revCodons = codonInfo()

    nonneutral = readModelFile(pd.read_csv(args.mssmodelfilename, sep = '\t'))

    stopCodons = ["TAA", "TAG", "TGA"]

    for codon in codonlist:
        aaDict = {}
        aaMuts = {}

        aa = revCodons[codon]
        synCodons = codons[aa]

        if codon in stopCodons:
            for secondC in codonlist:
                aaDict[secondC] = 0.0  ## stop codon
                aaMuts[secondC] = 0 ## stop codon but nonsyn
        else:
            for secondC in codonlist:
                if secondC == codon:
                    aaDict[secondC] = 1.0  ## same codon exactly
                    aaMuts[secondC] = -10 ## same codon
                elif secondC in synCodons:
                    if aa in nonneutral.keys():
                        aaSelectedCodons = nonneutral[aa]
                        if (secondC in aaSelectedCodons) & (codon in aaSelectedCodons):
                            aaDict[secondC] = args.synNonNetSrescaled  ## both codons found as nonneutral synonymous pairs
                            aaMuts[secondC] = 1 # [-, X, -]
                        else:
                            aaDict[secondC] = 1.0 ## synonymous but these two not nonneutral
                            aaMuts[secondC] = 2 # [-, -, X]
                    else:
                        aaDict[secondC] = 1.0  # synonymous but none of the codons are nonneutral
                        aaMuts[secondC] = 2 # [-, -, X]
                elif secondC in stopCodons:
                    aaDict[secondC] = 0.0  ## stop codon
                    aaMuts[secondC] = 0 ## stop codon but nonsyn
                else:
                    aaDict[secondC] = args.nonsynSrescaled  # nonsynonmous
                    aaMuts[secondC] = 0 # [X, -, - ]

        selectedDict[codon] = aaDict
        mutDict[codon] = aaMuts

    return selectedDict, mutDict



def maketreeshape(numSpecies):
    if numSpecies == 11:
        tree = '((((p1,(p5,p6)),(p4,((p8,p9),p7))),(p3,p10)),(p2,p11));'
        split_generations = {1: ['p1', 1, None],
                            2: ['p2', 0.3431, 'p1'],
                            3: ['p3', 0.4475, 'p1'],
                            4: ['p4', 0.5268, 'p1'],
                            5: ['p11', 0.6062, 'p2'],
                            6: ['p7', 0.6062, 'p4'],
                            7: ['p5', 0.6326, 'p1'],
                            8: ['p10', 0.6458, 'p3'],
                            9: ['p8', 0.6855, 'p7'],
                            10: ['p6', 0.7383, 'p5'],
                            11: ['p9', 0.7648, 'p8']}
    elif  numSpecies==5:
        tree = '((((p1,p5),p4),p3),p2);'
        split_generations = {1: ['p1', 1, None],
                        2: ['p2', 0.3, 'p1'],
                        3: ['p3', 0.5, 'p1'],
                        4: ['p4', 0.7, 'p1'],
                        5: ['p5', 0.8, 'p1']}

    else: #numspecies = 4
        tree = '(((p1,p4),p3),p2);'
        split_generations = {1: ['p1', 1, None],
                        2: ['p2', 0.3, 'p1'],  #0.3
                        3: ['p3', 0.6, 'p1'],
                        4: ['p4', 0.8, 'p1']}
    return tree,split_generations

def makefastafile(samples, filename, loc):

    with open('{}/{}.fa'.format(loc, filename), 'w') as o:
        for pop in samples.keys():
            o.write('>{}\n'.format(pop))
            o.write(str(samples[pop][0]) + "\n")

class chromosome():

    def __init__(self,sequence,fitness,args,mcounts):
        # consider a list of strings (codons)
        self.s = sequence
        self.fitstruct = args.fitnessstructure
        self.mutstruct = args.mutstructure
        self.mrate = args.mutrate
        self.ancestor = args.ancestor
        self.fitness = fitness
        self.mcounts = [mcounts[0],mcounts[1],mcounts[2]]


    #consider identifying the type of mutation and counting them  (maybe in debug mode)
    def mutate(self):
        """
            a function that changes s and recalculates fitness
        """

        # consider using  numpy.random.geometric or numpy.random.poisson
        pos = 0
        while True:
            npos = np.random.geometric(self.mrate)
            if npos < len(self.s):
                ## set position that mutates
                pos = npos

                ## identify old codon
                oldCodon = self.getOldCodon(pos) #write this
                holds = self.s
                while True: # keep sampling at pos until the new codon is not a stop codon
                    ## find new seqence
                    bps =['A', 'G', 'C', 'T']
                    bps.remove(self.s[pos:pos+1])
                    self.s = self.s[:pos] + np.random.choice(bps) + self.s[pos+1:]

                    ## update fitness
                    newCodon = self.fitnessfunction(pos, oldCodon)
                    # identify the mutation type and increment the main counter
                    if newCodon in ['TAG', 'TAA', 'TGA']:
                        self.s = holds
                    else:
                        break
                muttype = self.mutstruct[oldCodon][newCodon]
                # # identify the mutation type and increment the main counter
                # if newCodon in ['TAG', 'TAA', 'TGA']:
                #     muttype = 0
                # else:
                #     muttype = self.mutstruct[oldCodon][newCodon]
                self.mcounts[muttype] += 1
                mainmutationcounter[muttype] += 1
            else:
                break


    def fitnessfunction(self, mut,oldcodon):
        """
        recalculates fitness based on newCodon and ancestral codon just for mutations
        """
        stopCodons = ['TAG', 'TAA', 'TGA']

        anc, newSelf = self.findCodon(mut)

        if newSelf in stopCodons:
            return newSelf
        else:
            self.fitness = self.fitness * self.fitstruct[anc][newSelf] / self.fitstruct[anc][oldcodon]
        return newSelf

    def getOldCodon(self, i):
        """
        identify codon that mutation is in for ancestral and sequence
        """
        position = i % 3

        if position == 0:
            return self.s[i:i+3]
        elif position == 1:
            return self.s[i-1:i+2]
        elif position == 2:
            return self.s[i-2:i+1]
        else:
            print('error')

    def findCodon(self, i):
        """
        identify codon that mutation is in for ancestral and sequence
        """
        position = i % 3

        if position == 0:
            return self.ancestor[i:i+3], self.s[i:i+3]
        elif position == 1:
            return self.ancestor[i-1:i+2], self.s[i-1:i+2]
        elif position == 2:
            return self.ancestor[i-2:i+1], self.s[i-2:i+1]
        else:
            print('error')

    def __str__(self):
        return self.s

class population(list):
    """
    basically a list of chromosomes with some added functionality

    """

    def __init__(self, label, source, args):
        """
        if source is a sequence, then it is the ancestral sequence and all chromosomes are made as copies of it
        if source is a population then the new population is made by sampling from it at random
        """
        self.label = label
        self.mrate = args.mutrate
        self.popsize2 = args.popsize2
        self.args = args
        if isinstance(source,population):
            """
            randomly sample from source to make a new one
            or could duplicated it

            """
            for chrom in source:
                self.append(chrom)
        else:
            for i in range(self.popsize2):
                self.append(chromosome(source,1,args, [0,0,0]))

    def generation(self):
        """
        random sampling of the next generation based on the fitnesses of the chromosomes
        after each chromosome is sampled,  mutations are added and fitness is recalculated
        """
        newpop = []
        randomparentids = self.pickRandomParents()
        for i in randomparentids:
            child = chromosome(self[i].s,self[i].fitness,self.args,self[i].mcounts)
            child.mutate()
            newpop.append(child)
        self.clear()
        for i in range(self.popsize2):
            self.append(newpop[i])

    def pickRandomParents(self):
        fits = []
        selfids = []
        for i, g in enumerate(self):
            if g.fitness > 0.0:
                fits.append(g.fitness)  ## removes fitnesses of 0 from being able to pick ???
                selfids.append(i)

        fits = np.array(fits)

        return np.random.choice(selfids, size =self.popsize2, p = fits/fits.sum())

    """
        ## old code
        fits = np.empty(self.popsize2)
        selfids = np.arange(self.popsize2)
        for g in self:
            np.append(fits, g.fitness)

        return np.random.choice(selfids, size =self.popsize2, p = fits/fits.sum())
    """

    def sample(self, num):
        """
        return random chromosomes of number num from population
        """
        return np.random.choice(self, num)

    def checkpop(self, seqLen, gen):
        """
            check whatever should be true
            e.g. immediately after call to generation()  none should have fitness of 0
        """
        ## make sure no one has a fitness of 0
        for chrom in self:
            if chrom.fitness == 0.0:
                self.args.logfile.write('Individual has a fitness of Zero in Population ' + self.label + ' in genration {} \n'.format(str(gen)))

        ## make sure the popsize is constant
        assert len(self) == self.popsize2

        ## make sure the length of chromosome is right - check random chromsosme
        assert len(self[np.random.randint(self.popsize2)].s) == seqLen*3


class tree():
    """
        somehow represents a tree
        makes initial population

        when run() is called it runs the simulation and at the end returns the sample
    """
    #consider adding something that tracks fitness
    # e.g. every popsize2 gens  print the fitness of a random chromosome

    def __init__(self,args,ancestor):
        self.treestr = args.tree
        self.treedepth = args.treeDepth
        self.args = args
        self.split_generations, self.times = self.translateGens(args.split_generations)
        self.pop0 = population('p1', ancestor,args)
        self.pops = {}

    def translateGens(self, sg):
        newSG = {}
        times = []
        for key in sorted(sg.keys()):
            newTime = round(sg[key][1] * self.treedepth)
            newSG[newTime] = [sg[key][0], sg[key][2]]
            times.append(newTime)
        return newSG, times

    def sample(self):
        """
        samples sequences at the end of the run
        """
        samples = {}
        for pop in self.pops.keys():
            samples[pop] = self.pops[pop].sample(1)
        return samples

    def fitCheck(self):
        """
        picka  random chromosome from the population and output fitness
        """
        for pop in self.pops.keys():
            num = np.random.randint(self.args.popsize2)
            self.args.logfile.write("{}\t{}\t{}\n".format(pop, str(self.pops[pop][num].fitness), str(self.pops[pop][num].mcounts)))


    def run(self):
        """
        runs for treedepth generations
        make root population
        root = population
        """
        self.pops['p1'] = self.pop0
        gen = 1
        while gen < self.treedepth:
            """
            loop over generations,  adding populations as needed
            """
            if gen in self.times:
                splitPop = self.split_generations[gen]

                ## split populations
                self.pops[splitPop[0]] = population(splitPop[0], self.pops[splitPop[1]], self.args)
                #print(splitPop)

            for key in self.pops.keys():
                ## check to make sure fitness not zero
                if gen % (self.args.popsize2 * 4)  == 0:
                    self.pops[key].checkpop(self.args.aalength, gen)

                self.pops[key].generation()

            if self.args.debug == True:
                if gen % (self.args.popsize2 * 4) == 0:
                    self.args.logfile.write(str(gen) + '\n')
                    self.fitCheck()
                    self.args.logfile.write('\n')

            gen += 1
            #if gen % self.args.popsize2 == 0:
            #    print("generation",gen)

        sample = self.sample()
        self.args.logfile.write('\nTotal Mutations: ' + str(mainmutationcounter))
        return sample

def parseargs():
    parser = argparse.ArgumentParser("python makeSlimScript_JH.py",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-N", help="Population size (diploid)",dest="popsize",default=10,type=int)
    parser.add_argument("-L", help="Length of sequence (# amino acids)", dest="aalength",default=300,type=int)
    parser.add_argument("-Q", help="Theta = 4Nu", dest="theta",default=0.004,type=float)
    parser.add_argument("-s", help="Synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))", dest="synNonNetS",default=2,type=float)
    parser.add_argument("-y", help="Non-synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))", dest="nonsynS",default=10,type=float)
    parser.add_argument("-D", help="Tree depth in generations and in multiple of population size",dest="treeDepth",default=3500,type=int)
    parser.add_argument("-k", help="Number of species (4,5 or 11)",dest="numSpecies",default=4,type=int)
    parser.add_argument("-m", help="Model file name",dest="mssmodelfilename",default = "neutral-sim.tsv",type = str)
    parser.add_argument("-P", help="Prefix label for output files",dest="basename",required=True,type = str)
    parser.add_argument("-e", help="random number seed for picking alignment",dest="ranseed",type=int)
    parser.add_argument("-d", help="Debug mode", dest="debug", default=True, type=bool)
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
    args.popsize2 = args.popsize*2
    args.mutrate = args.theta/args.popsize2
    args.treeDepth = args.treeDepth * args.popsize
    #rescale the selection coefficients from 2Ns values to Slim values
    args.synNonNetSrescaled = synNonNetS = max(0.0,1.0 - (args.synNonNetS/(args.popsize2)))
    args.nonsynSrescaled = max(0.0,1.0 - (args.nonsynS/(args.popsize2)))
    if args.nonsynSrescaled <= 0.0:
        print("fitness error")
        exit()

    #update this when mutating
    global mainmutationcounter
    mainmutationcounter = [0,0,0]

    rundir  = "runs"
    #rundir  = "../runs"
    alignmentsdir = "bacterial_data/alignments"
    #alignmentsdir = "../alignments"
    randomgenebasepath = op.join(rundir,args.basename)

    # get ancestral sequence
    randomgenebasepath,genename, dnaStrand = writeCodonsToFile(randomgenebasepath,alignmentsdir)
    args.basename = args.basename + "_" + genename
    args.ancestor = makeAncestor(dnaStrand, args.aalength)
    args.logfile = open(rundir + '/' + args.basename + '_log.txt', 'w')
    #print(ancestor)

    # set tree shape
    args.tree, args.split_generations = maketreeshape(args.numSpecies)

    ## create selected dictionary
    args.fitnessstructure, args.mutstructure = createSelectedDictionary(args)

    # run the simulation
    sim = tree(args, args.ancestor)

    sampledsequences = sim.run()

    makefastafile(sampledsequences, args.basename, rundir)

    args.logfile.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        #sys.stderr.write("No arguments. Use -h  or --help for help menu")
        main(['-h'])
    else:
        main(sys.argv[1:])