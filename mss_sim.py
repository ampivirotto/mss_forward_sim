"""
    forward simulation under an MSS model
    officially diploid,  but only in doubling of # of chromosomes
        throughout popsize refers to diploid # and popsize2 to haploid number
    assumes selection at the haplotype/chromosome level

    lifecycle:
        a chromosome is sampled at random from one of the possible parents, based on the parent fitnesses
        after it is sampled new mutations are added and its fitness is calculated

    discrete generations
    
    selection model - stabilizing 
        a random chromosome is sampled and assigned a fitness value of 1
        this is set as the ancestral chromosome
        all mutations away from this cause lower fitness
        all mutations toward this cause higher fitness
        fitness can never be greater than 1
    
    population initiation
        the root population is started with a bunch of chromosomes,  all with copies of the ancestral chromosome 
        burnin period 1 proceeds until population mean fitness decline has slowed or stopped
        burnin period 2 proceeds until the entire population is descended from a single chromosome that was present 
            in the generation that burnin period 2 started. 
            this establishes a time that should actually be the tmrca for ree for the sequences sampled at the end 

    phylogeny
        a fixed phylogeny over 100000 generations is simulated
        population splitting happens by randomly sampling (with replacement) a population that becomes the sister to the population on the new branch of the tree
        

"""
import os
import os.path as op
import argparse
import numpy as np
import sys
import time

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

def readInGeneFile(alignment_location,gene = None):
    """
    given gene selected, returns the alignment data for that gene
    """
    if gene==None:
        gene = identifyRandomGene(alignment_location)
    with open(op.join(alignment_location,gene)) as f:
        species = {}
        dna = ''
        spec = ''
        for line in f:
            if line.startswith('>'):
                assert len(dna) % 3 == 0
                if spec != '':
                    species[spec] = dna
                spec = line.lstrip('>').rstrip('\n')
                dna = ''
            else:
                dna += line.strip('\n')
    return species, gene

def createCodonSequence(alignment_location,gene = None):
    """
    given alignments for randomly selected gene,
    it turns all alignments into single strand of DNA excluding codons with missing bps or stop codons
    """
    species, gene = readInGeneFile(alignment_location,gene = gene)
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
                    # if codon == 'ATG':
                    #     codon = ''
                    #     continue
                    allDNA += codon
                    codon = ''
                codon += bp
    return allDNA, gene


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

def getModelCodonPairs(lls):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLT': 'E', 'TYR': 'Y', 'MET': 'M'}
    sdict = {}
    assert len(lls) == 87, "missing 1 or more codon pairs,  should be 87 of them "
    for ls in lls:
        [aa,codon1,codon2,selneu] = ls.strip().split()
        if len(aa) == 3:
            A1 = d[aa]
        else:
            A1 = aa.upper()
            assert len(A1) == 1
        if A1 in sdict:
            sdict[A1].append([codon1,codon2,selneu])
        else:
            sdict[A1] = [[codon1,codon2,selneu]]
    return sdict

def readModelFile(fn):
    """
        fn is the model file
        if "CODON1" and "CODON2" are in the header, it is a file of pairs of codons, in which case call getModelCodonPairs()
        otherwise it is a file of codons grouped by amino acid
    """
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLT': 'E', 'TYR': 'Y', 'MET': 'M'}
    lls = open(fn,'r').readlines()
    while len(lls[-1]) < 2:
        lls = lls[:-1]
    if "CODON1" in lls[0] and "CODON2" in lls[0]:
        sdict = getModelCodonPairs(lls[1:])
        return sdict,"codonpairs"
    lls = lls[1:] # skip first line
    sdict = {}
    for ls in lls:
        [aa,codon,selneu] = ls.strip().split()
        A1 = d[aa]
        if selneu == "SELECTED":
            if A1 in sdict:
                sdict[A1].append(codon)
            else:
                sdict[A1] = [codon]
    return sdict,"aminoacidsets"

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
    # aalist.remove('STOP')
    codonlist = []

    revCodons  = {}

    for aa in aalist:
        for cd in codons[aa]:
            codonlist.append(cd)
            revCodons[cd] = aa

    return codons, aalist, codonlist, revCodons

def createSelectedDictionary(args):
    """
    if some structure needs to be built that represents codon fitnesses efficiently,  this is the place for
    mutDict :  0,1 or 2  for nonsynonymous,  synonymous-selected, synonymous-neutral
    """
    selectedDict = {}
    mutDict = {}

    codons, aalist, codonlist, revCodons = codonInfo()
    nonneutral,modeltype = readModelFile(args.mssmodelfilename)

    stopCodons = ["TAA", "TAG", "TGA"]
    if modeltype=="aminoacidsets":
        for codon in codonlist:
            aaDict = {}
            aaMuts = {}

            aa = revCodons[codon]
            synCodons = codons[aa]

            if codon in stopCodons:
                for secondC in codonlist:
                    aaDict[secondC] = 0.0  ## stop codon
                    aaMuts[secondC] = 3 ## stop codon but nonsyn
            else:
                for secondC in codonlist:
                    if secondC == codon:
                        aaDict[secondC] = 1.0  ## same codon exactly
                        aaMuts[secondC] = -10 ## same codon
                    elif secondC in synCodons:
                        if aa in nonneutral.keys():
                            aaSelectedCodons = nonneutral[aa]
                            if (secondC in aaSelectedCodons) & (codon in aaSelectedCodons):
                                aaDict[secondC] = args.SynSel_s_rescaled  ## both codons found as nonneutral synonymous pairs
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
                        aaDict[secondC] = args.NonSyn_s_rescaled  # nonsynonmous
                        aaMuts[secondC] = 0 # [X, -, - ]

            selectedDict[codon] = aaDict
            mutDict[codon] = aaMuts
    else: #"codonpairs"
        tempd = {}
        for aa in nonneutral:
            for [codon1,codon2,val] in nonneutral[aa]:
                if codon1 in tempd:
                    tempd[codon1][codon2] = val
                else:
                    tempd[codon1] = {codon2:val}
                if codon2 in tempd:
                    tempd[codon2][codon1] = val
                else:
                    tempd[codon2] = {codon1:val}
        for codon1 in codonlist:
            aaDict = {}
            aaMuts = {}
            aa1 = revCodons[codon1]
            synCodons = codons[aa1]
            for codon2 in codonlist:
                aa2 = revCodons[codon2]
                if codon1 in stopCodons or codon2 in stopCodons:
                    aaDict[codon2] = 0.0  ## stop codon
                    aaMuts[codon2] = 3 ## stop codon
                elif codon2 == codon1:
                    aaDict[codon2] = 1.0  ## same codon exactly
                    aaMuts[codon2] = -10 ## same codon
                elif aa1 != aa2:
                    aaDict[codon2] = args.NonSyn_s_rescaled  # nonsynonmous
                    aaMuts[codon2] = 0 # [X, -, - ]
                else:
                    assert codon2 in synCodons
                    if tempd[codon1][codon2] == "NEUTRAL":
                        aaDict[codon2] = 1.0  # synonymous but none of the codons are nonneutral
                        aaMuts[codon2] = 2 # [-, -, X]
                    else:
                        assert tempd[codon1][codon2] == "SELECTED"
                        aaDict[codon2] = args.SynSel_s_rescaled  ## both codons found as nonneutral synonymous pairs
                        aaMuts[codon2] = 1 # [-, X, -]


            selectedDict[codon1] = aaDict
            mutDict[codon1] = aaMuts

    return selectedDict, mutDict


def maketreeshape(numSpecies):
    """
    use fixed phylogeny, given number of species 
    """
    if numSpecies == 11:
        tree = '((((p1,(p5,p6)),(p4,((p8,p9),p7))),(p3,p10)),(p2,p11));'
        split_generations = {1: ['p1', 1, None],
                            2: ['p2', 0.0, 'p1'],
                            3: ['p3', 0.1475, 'p1'],
                            4: ['p4', 0.3268, 'p1'],
                            5: ['p11', 0.40, 'p2'],
                            6: ['p7', 0.4062, 'p4'],
                            7: ['p5', 0.4326, 'p1'],
                            8: ['p10', 0.5458, 'p3'],
                            9: ['p8', 0.5855, 'p7'],
                            10: ['p6', 0.8383, 'p5'],
                            11: ['p9', 0.9148, 'p8']}
        mean_branches_root_to_tip = 4.091
    elif  numSpecies==5:
        tree = '((((p1,p5),p4),p3),p2);'
        split_generations = {1: ['p1', 1, None],
                        2: ['p2', 0.0, 'p1'],
                        3: ['p3', 0.25, 'p1'],
                        4: ['p4', 0.5, 'p1'],
                        5: ['p5', 0.75, 'p1']}
        mean_branches_root_to_tip = 2.8
    else: #numspecies = 4
        tree = '(((p1,p4),p3),p2);'
        split_generations = {1: ['p1', 1, None],
                        2: ['p2', 0.0, 'p1'],  #0.3
                        3: ['p3', 0.4, 'p1'],
                        4: ['p4', 0.8, 'p1']}
        mean_branches_root_to_tip = 2.25
    return tree,split_generations,mean_branches_root_to_tip 

def makefastafile(samples, fn):
    f = open(fn,'w')
    for pop in samples.keys():
        f.write('>{}\n'.format(pop))
        f.write(str(samples[pop]) + "\n")
    f.close()
    return

class chromosome():
    """
    a chromosome is an 'individual' in the population
    it has a DNA sequence from which its fitness is determined 
    """

    def __init__(self,sequence,fitness,args,mcounts,ancestornumber):
        """
        mcounts :  nummutationtypes positions 
            0,1,2 or 3 for nonsynonymous,  synonymous-selected, synonymous-neutral, and STOP 
        re ancestornumber:
          is just the index of the chromosome in the list
          it is reset at the beginning of burn2
          during burn2 ancestornumber is updated for each chromosome, by copying the value from the ancestor the chromosome was copied from
          when all values of ancestornumber are the same,  all the chromosomes are descended from that ancestor in that generation
        """
        global nummutationtypes
        self.s = sequence
        self.fitstruct = args.fitnessstructure
        self.mutstruct = args.mutstructure
        self.mrate = args.mutrate
        self.mrateinverse = 1.0/self.mrate  # use with exponential random variables for mutation locations because numpy.random.exponential() takes inverse of the rate as the parameter for some reason
        self.ancestor = args.ancestor
        self.debugmutloc = args.debug
        self.fitness = fitness
        self.mcounts = [mcounts[i] for i in range(nummutationtypes)]
        self.ancestornumber = ancestornumber

    def mutate(self):
        """
            a function that changes the sequence s and recalculates fitness
            it uses exponential to get the distance to the next mutation as an approximation for geometric
            sites are selected by jumping along the chromosome (distance sampled from exponential)
            fitness is not updated until all of the sites in a codon have been changed
                usually each mutation is in its own codon, but sometimes 2 or 3 changes occur in the same codon
                code is a kind of kludgy, in order to keep a list of all changes in the current codon
            
        """
        global mutationlocations # use in debug mode 
        pos = 0 # a position that mutates  (if not past the end of the sequence)	
        lastcodonpos = -1
        
        while True:
            # distance_to_mut = np.random.geometric(self.mrate) # geometric a bit slower than exponential
            while True: # gets one or more changes in a codon
                # exponential is faster than geometric,  but can return 0 
                distance_to_mut = int(np.random.exponential(self.mrateinverse)) 
            ## set position that mutates
                pos += distance_to_mut
                codonpos = pos //3
                if lastcodonpos == -1:
                    mutlocs = [pos]
                    lastcodonpos = codonpos
                    pos += 1 # increment to value that is next possible position to mutate 
                elif codonpos != lastcodonpos:
                    lastcodonpos = codonpos
                    pos += 1 # increment to value that is next possible position to mutate 
                    break
                else:
                    mutlocs.append(pos)
                    pos += 1 # increment to value that is next possible position to mutate 
            if mutlocs[-1] < len(self.s):
                ## identify old codon
                oldCodon = self.getOldCodon(mutlocs[0]) 
                for ml in mutlocs:
                    bps =['A', 'G', 'C', 'T']
                    bps.remove(self.s[ml:ml+1])
                    self.s = self.s[:ml] + np.random.choice(bps) + self.s[ml+1:]
                ## update fitness
                newCodon = self.fitnessfunction(mutlocs[0], oldCodon)
                if self.debugmutloc:
                    for ml in mutlocs:
                        mutationlocations[ml] += 1
                muttype = self.mutstruct[oldCodon][newCodon]
                self.mcounts[muttype] += 1
                mainmutationcounter[muttype] += 1
                if pos > len(self.s):
                    break
                else: # reset mutlocs to contain only the last exponential jump position (before it was incremented)
                    mutlocs = [pos-1] # the next mutation location (pos was previously incremented to the start of the next interval)
            else:
                break
    def resetmcounts(self):
        global nummutationtypes
        self.mcounts = [0 for i in range(nummutationtypes)]


    def fitnessfunction(self, mut,oldcodon):
        """
        recalculates fitness based on newCodon and ancestral codon just for mutations
        the ancestor was defined as having a fitness of 1
        for x  codons  fitness is a product of x values 
        the ancestor was defined as having a fitness of 1
        for a change from oldcodon to newcodon the fitness is updated by multiplying 
        fitness times the fitness associated with a change from ancestor to the new codon (i.e. the absolute fitness of the new codon at that position)
        and by dividing by the fitness associate with a change from the ancestor to the old codon 
        """
        stopCodons = ['TAG', 'TAA', 'TGA']
        anc, newSelf = self.findCodon(mut)
        assert newSelf != oldcodon


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
            print('findcodon() error',i)
            exit()

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
        global nummutationtypes
        self.label = label
        self.mrate = args.mutrate
        self.popsize2 = args.popsize2
        self.args = args
        if isinstance(source,population):
            """
            copy each chromosome from source and put it in the new population
            """
            for chrom in source:
                self.append(chrom)
        else:# at the beginning, fill up the pouplation with chromosomes made from the ancestor 
            for i in range(self.popsize2):
                self.append(chromosome(source,1,args, [0 for j in range(nummutationtypes)],i))

    def generation(self):
        """
        random sampling of the next generation based on the fitnesses of the chromosomes
            make array of fitnesses
            get list of unique values and indices for these values
            get expected frequencies
            sample randomparents using multinomial
        after each chromosome is sampled,  mutations are added and fitness is recalculated
        replace the old population with the new sampled chromosomes
        """
        fits = np.array([c.fitness for c in self],dtype=float)
        unique_vals,indices,counts = np.unique(fits, return_counts=True,return_inverse=True)
        unique_indices = []
        for i in range(len(unique_vals)):
            unique_indices.append(np.where(indices == i)[0])
        numfits = unique_vals.shape[0]
        if numfits == 1:
            randomparentids = [np.random.choice(unique_indices[0],size = self.args.popsize2,replace = True)]
        else:
            expfreqs = (unique_vals*counts)/len(fits)/fits.mean()
            expfreqs = [v if v <= 1.0 else 1.0 for v in expfreqs] # can get a value slightly greater than 1  when a fitness is 0
            samplecounts = np.random.multinomial(self.args.popsize2,expfreqs,1)
            randomparentids = [np.random.choice(unique_indices[i],size=samplecounts[0,i],replace=True) for i in range(numfits)]
            if self.args.debug:
                for parentgroup in randomparentids:
                    for i in parentgroup:
                        assert self[i].fitness > 0 
        newpop = []       
        for parentgroup in randomparentids:
            for i in parentgroup:
                #copy the chromosome.  this is much, much faster than deepcopy
                child = chromosome(self[i].s,self[i].fitness,self.args,self[i].mcounts,self[i].ancestornumber)
                child.mutate()
                newpop.append(child)
        self.clear()
        for child in newpop:
            self.append(child)
        return numfits
    
    def checkancestors(self):
        anc0 = self[0].ancestornumber
        allthesame = all(anc0 == c.ancestornumber for c in self)
        return allthesame
    
    def sampleindividual(self, num):
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
        # for chrom in self:
        #     if chrom.fitness == 0.0:
        #         self.args.logfile.write('Individual has a fitness of Zero in Population ' + self.label + ' in genration {} \n'.format(str(gen)))

        ## make sure the popsize is constant
        assert len(self) == self.popsize2

        ## make sure the length of chromosome is right - check random chromsosme
        assert len(self[np.random.randint(self.popsize2)].s) == seqLen*3

    def reset_mutation_counts(self):
        for c in self:
            c.resetmcounts()

class tree():
    """
        Represents a tree of populations
        makes initial population
        when tree.run() is called it runs the simulation and at the end returns the sample
    """

    def __init__(self,args,ancestor):
        self.treestr = args.tree
        self.args = args
        self.split_generations, self.times = self.translateGens(args.split_generations)
        self.pop0 = population('p1', ancestor,args)
        self.pops = {}

    def translateGens(self, sg):
        newSG = {}
        times = []
        for key in sorted(sg.keys()):
            # newTime = round(sg[key][1] * self.args.treeDepth + self.args.burnin)
            newTime = round(sg[key][1] * self.args.treeDepth )
            newSG[newTime] = [sg[key][0], sg[key][2]]
            times.append(newTime)
        return newSG, times


    def samplefrompops(self):
        """
        samples sequences at the end of the run
        """
        stopCodons = ['TAG', 'TAA', 'TGA']
        samples = {}
        for pop in self.pops.keys():
            while True: # avoid sampling an individual with a stop codon, will hang if all individuals have a stop codon
                temp =self.pops[pop].sampleindividual(1)[0]
                notok = True in [codon in stopCodons for codon in [temp.s[i:i+3] for i in range(0, len(temp.s), 3)]]
                if notok is False:
                    break
            samples[pop] = temp
        return samples

    def fitCheck(self):
        """
        picka  random chromosome from the population 
        #write fitnesses to log file
        return mean fitness
        """
        meanfit = 0
        popkeys = self.pops.keys()
        for pop in popkeys:
            num = np.random.randint(self.args.popsize2)
            temp = self.pops[pop][num].fitness
            # self.args.logfile.write("{}\t{}\t{}\n".format(pop, str(temp), str(self.pops[pop][num].mcounts)))
            meanfit += temp
        return meanfit/len(popkeys)

    def fitmutsummary(self):
        """
        pick a  random chromosome from each population
        return a list  [meanfitness, list of fitness,  list of mcount lists]
        """
        meanfit = 0
        popkeys = self.pops.keys()
        fitlist = []
        mcountlist = []
        assert len(popkeys)== self.args.numSpecies
        for pop in popkeys:
            num = np.random.randint(self.args.popsize2)
            temp = self.pops[pop][num].fitness
            # self.args.logfile.write("{}\t{}\t{}\n".format(pop, str(temp), str(self.pops[pop][num].mcounts)))
            meanfit += temp
            fitlist.append(temp)
            mcountlist.append(self.pops[pop][num].mcounts)
        return meanfit/len(popkeys),fitlist,mcountlist

    def summarize_results(self,starttime):
        global mainmutationcounter
        global nummutationtypes
        mnames = ["NonSynonymous ","Synonymous_Sel","Synonymous_Neu","Stop"]
        meanfit,fitlist,mcountlist = self.fitmutsummary()
        rf = open(self.args.resultsfilename,'w')
        rf.write("mss_sim\n\narguments:\n")
        # for arg in vars(self.args):
        #     rf.write("\t{}: {}\n".format(arg, getattr(self.args, arg)))
        for arg in vars(self.args):
            if arg=="fitnessstructure" or  arg=="mutstructure":
                if self.args.debug:
                    rf.write("\t{}: {}\n".format(arg, getattr(self.args, arg)))
                else:
                    rf.write("\t{}: {}\n".format(arg, " - printed only in debug mode")) # quite large
            else:
                rf.write("\t{}: {}\n".format(arg, getattr(self.args, arg)))

        rf.write("\nFinal Mean Fitness: {:.4g}\n".format(meanfit))
        rf.write("\nMean number of fitness values each generation: {:.1f}\nMean number of individuals per fitness value: {:.1f}\n".format(self.args.meannumfits,self.args.popsize2/self.args.meannumfits))
        rf.write("\nSampled Individual Fitnesses: {}\n".format(fitlist))
        rf.write("\nSampled Individual Mutation Counts ([NonSyn,Syn-Sel,Syn-Neu,STOP]): {}\n".format(mcountlist))
        rf.write("\nMutation Total Counts/Rates (per effective bp)\n")
        rf.write("\tmutation_type\tcounts\teffective_#bp\tproportions\tmutations_per_effective_bp:\n")
        totsum = sum(mainmutationcounter)
        effectivenumbp = []
        mutperebp = []
        propebp = []
        for i in range(nummutationtypes):
            propebp.append(mainmutationcounter[i]/totsum)
            effectivenumbp.append(3*self.args.aalength*propebp[i])
            mutperebp.append(np.nan if effectivenumbp[i] == 0 else mainmutationcounter[i]/effectivenumbp[i])
            rf.write("\t{}\t{:d}\t{:.0f}\t{:.3g}\t{:.3g}\n".format(mnames[i],mainmutationcounter[i],effectivenumbp[i],propebp[i],mutperebp[i]))
        rf.write("\nSubstitutions per gene copy Counts/Rates (per effective bp)\n")
        rf.write("\tsubstitution_type\tmean counts per gene\tsubstitutions_per_effective_bp\tsubstitutions_per_effective_bp_per_branch\tsubstitutions_per_effective_bp_per_generation\n")
        subsum = [0 for i in range(nummutationtypes)]
        for mc in mcountlist:
            for i in range(nummutationtypes):
                subsum[i] += mc[i]
        for i in range(nummutationtypes): # take the mean count per sampled chromosome
            subsum[i] /= self.args.numSpecies
        subrates = []
        totalnumgen = self.args.treeDepth + self.args.burn2_generation_time
        for i in range(nummutationtypes):
            subrates.append(np.nan if effectivenumbp[i] == 0 else subsum[i]/effectivenumbp[i])
            self.args.mean_branches_root_to_tip 
            rf.write("\t{}\t{:.1f}\t{:.3g}\t{:.3g}\t{:.3g}\n".format(mnames[i],subsum[i],subrates[i],subrates[i]/self.args.mean_branches_root_to_tip ,subrates[i]/totalnumgen))
        #positions in mutation, substitution arrays of length 4 
        NonSpos = 0
        SynSelpos = 1
        SynNeupos = 2
        STOPpos = 3
        rf.write("\nRate Ratios:\n")
        total_syn_effectivenumbp = effectivenumbp[SynSelpos]+ effectivenumbp[SynNeupos]
        total_syn_sub_rate = (subsum[SynSelpos] + subsum[SynNeupos])/total_syn_effectivenumbp
        rf.write("\tdN/dS (Nonsynonymous/Synonymous (selected and neutral))\t{:.3g}\n".format(np.nan if total_syn_sub_rate == 0 else subrates[NonSpos]/total_syn_sub_rate))
        rf.write("\tdN/dSn (Nonsynonymous/Synonymous_Neu))\t{:.3g}\n".format(np.nan if subrates[NonSpos] == 0 else subrates[NonSpos]/subrates[SynNeupos]))
        rf.write("\tdSs/dSn (Synonymous_Sel/Synonymous_Neu)\t{:.3g}\n".format(np.nan if subrates[SynSelpos] == 0 else subrates[SynSelpos]/subrates[SynNeupos]))
        totaltime = time.time()-starttime
        rf.write("\ntotal time: {}\n".format(time.strftime("%H:%M:%S",time.gmtime(totaltime))))
        rf.close()

    def run_burn1(self):
        """
        loop over generations until meanfitness appears to stop going down
        this is very crude and does not really find the point when mutation stops declining (if at all)
        but allows process to accumulate a number of mutations so that the mutation counter does not include counts on the ancestral chromosome with fitness 1 
        """
        maxburn1gen = 10000 # fail stopping point
        gen = 0
        numgens_check_burn1 = 300 # just picked a value 
        numgenssplit = numgens_check_burn1//3
        meanfitlist = [] # will include numgens_check_burn1 values.  compare the mean in first third to mean in last third 
        while True:
            self.pop0.generation()
            meanfit = sum([c.fitness for c in self.pop0])/len(self.pop0)
            meanfitlist.append(meanfit)
            gen += 1
            if len(meanfitlist) > numgens_check_burn1:
                meanfitlist.pop(0)
                # compare first 3rd of list to last 3rd of list
                if sum(meanfitlist[:numgenssplit])/numgenssplit < sum(meanfitlist[numgens_check_burn1-numgenssplit:])/numgenssplit:
                    break
            if gen > maxburn1gen:
                print('burn 1 has exceeded {} generations. fitness continues to decline'.format(maxburn1gen))
                break
                # exit()
        self.pop0.reset_mutation_counts()      
        return gen,sum(meanfitlist)/numgens_check_burn1
    
    def run_burn2(self):
        """
        burn to the point that all chromosomes are descended from a common ancestor
        This sets the generation number at the base of the tree
        total time will then be self.args.burn2_generation_time + self.args.treeDepth
        """
        gen = 0
        for i,c in enumerate(self.pop0):
            c.ancestornumber = i
        while self.pop0.checkancestors()== False:
            self.pop0.generation()
            gen += 1        
        return gen

    def run(self):
        """
        calls run_burn1(), run_burn2() and then  runs for treedepth generations
        """

        self.args.burn1_generation_time,self.args.burn1_mean_fitness = self.run_burn1()
        self.args.burn2_generation_time = self.run_burn2()
        self.args.mutationexpectation_adjusted_for_burn2 = (self.args.mutationexpectation * self.args.treeDepth)/(self.args.treeDepth + self.args.burn2_generation_time)
        gen = 0
        countpopgens = 0
        self.pops['p1'] = self.pop0
        while gen < self.args.treeDepth:
            """
            loop over generations,  adding populations as needed
            """
            if gen in self.times:
                splitPop = self.split_generations[gen]
                ## split populations
                self.pops[splitPop[0]] = population(splitPop[0], self.pops[splitPop[1]], self.args)

            for key in self.pops.keys():
                numdifferentfitnessvalues = self.pops[key].generation()
                self.args.meannumfits += numdifferentfitnessvalues
                countpopgens +=1

            
            if self.args.debug == True:
                if gen % (self.args.popsize2 * 4) == 0:
                    self.pops[key].checkpop(self.args.aalength, gen)
                    # self.args.logfile.write(str(gen) + '\n')
                    meanfit = self.fitCheck()
                    # self.args.logfile.write('\n')

            gen += 1
            if self.args.debug and gen % self.args.popsize2 == 0:
                meanfit = self.fitCheck()
                print("generation {} ({:.1f}%)  # populations: {}  mean fitness: {:.4f}".format(gen,100*gen/self.args.treeDepth,len(self.pops.keys()),meanfit))
        self.args.meannumfits /= countpopgens
        sample = self.samplefrompops()
        # self.args.logfile.write('\nTotal Mutation counts [NonSyn,Syn-Sel,Syn-Neu]: ' + str(mainmutationcounter))
        return sample

def parseargs():
    parser = argparse.ArgumentParser("python makeSlimScript_JH.py",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-A", help="path for folder containing bacterial alignment files",dest="bacaligndir",required=True,type = str)
    parser.add_argument("-b", help="base filename for output files (name only, no directories)",dest="basename",required=True,type = str)
    parser.add_argument("-d", help="Debug mode", dest="debug", default=False,action="store_true")
    parser.add_argument("-e", help="random number seed for simulation (and for picking alignment if needed)",dest="ranseed",type=int)
    parser.add_argument("-F", help="directory path for output fasta file (default is same as for results file)",dest="fdir",type = str)
    parser.add_argument("-g", help="bacterial gene name, optional - if not used a random gene is selected",dest="genename",type=str)
    parser.add_argument("-k", help="Number of species (4,5 or 11)",dest="numSpecies",default=4,type=int)
    parser.add_argument("-L", help="Length of sequence (# amino acids)", dest="aalength",default=300,type=int)
    parser.add_argument("-m", help="Model file path",dest="mssmodelfilename",required = True,type = str)
    parser.add_argument("-N", help="Population size (diploid)",dest="popsize",default=10,type=int)
    parser.add_argument("-R", help="directory path for results file",dest="rdir",default = ".",type = str)
    parser.add_argument("-s", help="Synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))", dest="SynSel_s",default=2,type=float)
    parser.add_argument("-y", help="Non-synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))", dest="NonSyn_s",default=10,type=float)
    parser.add_argument("-u", help="expected number of neutral mutations per site, from base of tree", dest="mutationexpectation",default=0.5,type=float)
    return parser

def main(argv):
    starttime = time.time()
    parser = parseargs()
    if argv[-1] =='':
        argv = argv[0:-1]
    args = parser.parse_args(argv)
    if args.ranseed != None:
        np.random.seed(args.ranseed)
    if args.numSpecies not in [4,5,11]:
        print ("error: -p (# of species) must be 4,5 or 11")
        exit()
    args.meannumfits = 0
    args.popsize2 = args.popsize*2
    args.treeDepth = 100000 # fixed at a specific value # previously scaled by population size  args.treeDepth * args.popsize
    args.mutrate = args.mutationexpectation/args.treeDepth  # got rid of using theta 4Nu,  as not really relevant here 
    #rescale the selection coefficients from 2Ns values to Slim values
    args.SynSel_s_rescaled = max(0.0,1.0 - (args.SynSel_s/(args.popsize2)))
    args.NonSyn_s_rescaled = max(0.0,1.0 - (args.NonSyn_s/(args.popsize2)))
    if args.NonSyn_s_rescaled <= 0.0:
        print("fitness error")
        exit()
    curdir = os.getcwd()
    #if path is a string that can be separated into folders,  and one or more of them do not exist,  this will create them
    try:# create folder(s) if needed. When running lots of jobs, a dir may not exist at one moment, but then does exist in the next because of another job running, so use try/except
        curdir = os.getcwd()
        normalized_path = op.normpath(args.rdir)
        dirs = normalized_path.split(os.sep)
        for d in dirs:
            if op.exists(d) == False:
                os.mkdir(d)
            os.chdir(d)
        os.chdir(curdir)
    except:
        pass


    if args.fdir== None:
        args.fdir=args.rdir
    else:
        #if path is a string that can be separated into folders,  and one or more of them do not exist,  this will create them
        try:# create folder(s) if needed. When running lots of jobs, a dir may not exist at one moment, but then does exist in the next because of another job running, so use try/except
            curdir = os.getcwd()
            normalized_path = op.normpath(args.fdir)
            dirs = normalized_path.split(os.sep)
            for d in dirs:
                if op.exists(d) == False:
                    os.mkdir(d)
                os.chdir(d)
            os.chdir(curdir)        
        except:
            pass  
        
    
    #update this when mutating
    global mainmutationcounter
    global nummutationtypes
    nummutationtypes = 4
    mainmutationcounter = [0 for i in range(nummutationtypes)]  #positions 0,1,2 or 3 for nonsynonymous,  synonymous-selected, synonymous-neutral, and STOP

    # when debugging for checking distribution of mutation locations
    global mutationlocations
    mutationlocations = [0 for i in range(3*args.aalength)]

    # get ancestral sequence
    if args.genename != None:
        singleGene = True
    dnaStrand,genefilename = createCodonSequence(args.bacaligndir,gene=args.genename)# if args.genename is None,  then a random gene is picked
    args.genename = genefilename[:genefilename.find('_')]
    args.ancestor = makeAncestor(dnaStrand, args.aalength)

    #set file names 
    # args.logfilename = op.join(args.rdir,args.basename +  "_" + args.genename + '_log.txt')
    # if os.path.exists(args.logfilename):
    #     args.logfilename = '{}(1)_log.txt'.format(args.logfilename[:-8])
    # args.logfile = open(args.logfilename, 'w')
    args.resultsfilename = op.join(args.rdir,args.basename +  "_" + args.genename + '_results.txt')
    if (os.path.exists(args.resultsfilename)) or (singleGene == True):
        args.resultsfilename =  '{}_{}_results.txt'.format(args.resultsfilename[:-12], str(args.ranseed))
    args.fastafilename = op.join(args.rdir, '{}_{}.fa'.format(args.basename, args.genename))
    if (os.path.exists(args.fastafilename)) or (singleGene == True):
        args.fastafilename = op.join(args.rdir, '{}_{}_{}.fa'.format(args.basename, args.genename, str(args.ranseed)))

    # set tree shape
    args.tree, args.split_generations,args.mean_branches_root_to_tip = maketreeshape(args.numSpecies)

    ## create selected dictionary
    args.fitnessstructure, args.mutstructure = createSelectedDictionary(args)

    # run the simulation
    sim = tree(args, args.ancestor)
    sampledsequences = sim.run()
    sim.summarize_results(starttime)
    if args.debug:
        print("mutation counts by base position\n",mutationlocations)

    totaltime = time.time()-starttime
    # args.logfile.write("\ntotal time: {}\n".format(time.strftime("%H:%M:%S",time.gmtime(totaltime))))
    # args.logfile.close()
    makefastafile(sampledsequences, args.fastafilename)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        main(['-h'])
    else:
        main(sys.argv[1:])
