"""
    5/22/2023
    redid some organization and chunks of code to tidy things up

    5/19/2023

    realized that the model files we were trying where internally inconsistent. 
    e.g. if all 4 codons of an AA are selected,  but there is only 1 selection coefficient
    then there is no way to specify the selected directions in a way that makes sense
    e.g. if  x> y and y> z  by the same amounts,  how can x > z by that same amount ?? 

    to avoid this,  the following rules must apply:
    any AA set can be all neutral
    k=2
        can be a selected pair, without any neutral pairs 
        or can have a neutral pair 
    for k=3,4,6
        can have a single neutral set of k-1 neutrals and all rates from the kth codon to any of the neutrals is selected 
        
    for k=4,6
        can have two sets of neutrals, each with 2 (or 3 in case of k=6) codons 
        all rates between the two sets are selected 
    5/14/2023
    record info about mutations that accumulate
    5/10/2023
        tried coding it with lists of codons for the sequence and ancestor  but this did not speed things up. 

    5/5/2023 

    implemented a method for starting with an ancestral sequence that is close the equilibrium for preferred
    and unpreferred codons
    see makeAncestor()

    5/1/2023
    Handling selection for synonymous changes
    If all syn changes have the same value in the selection structure dictionary then
    if anc is a different codon for the same amino acid this ratio will often be 1. 
    i.e. self.fitstruct[anc][newcodon] / self.fitstruct[anc][oldcodon] == 1
    Also if anc is for a different amino acid than oldcodon and newcodon,  then the ratio will be 1.  
    But we want every type 1 mutation to give us a fitness change. 
    And we want some of these type 1 mutations to be favored. 

    
    To deal with this,  the codons as ordered in codondic and codonlist are assumed to be
      in fitness order,  such that the high fitness ones come first. 
      Any pair of syn codons, codon1 and codon2  with indices in codondic of c1 and c2,  
      And if c1 < c2,  then they represent a pair for which codon1 has the higher fitness. 
      So if the mutation is from codon1 to codon2, then it  is to a lower fitness with factor  args.SynSelDel_s_rescaled
      However if the change is from codon2 to codon1,  then it is to a higher fitness with factor  args.SynSelFav_s_rescaled

      to make this work,  we have to set the ancestral sequence with the high fitness codons to begin with
      otherwise there is a lot of evolution in which the random low fitness codons get replaced with high fitness ones
      and the fitness goes up and up for awhile. 
    

      
    The alternative to such schemes would be to have a ranking among the synonymous codons for an amino acid 
        which would require having more fitness values for synonymous changes , and thus more fitness values among chromosomes 
        
    adding adaptation 
        option -w gives the rate at which the population experiences a change in the optimal sequence
        e.g. -w 1e-3  means that the probability each generation of a change is 1e-3
        when there is a change,  a random position is selected, and a new codon for a different amino acid is 
        put in the population ancestor at that position
        
        need to keep a new data structure that has a list of all these. by position 
            dictionary keys are codon positions,  values are the new beneficial amino acid 

        if a mutation type is nonsynonymous (class 0) then look up to see if it is an adaptive change and count it as such

        sites change optimal amino acid 

        need to be able to count adaptive mutations so we can count them as we do other mutations 
        

        for each population,  simulate the generation numbers in which there is a change 

        each population gets its own ancestor 
        when there is a change,  a random codon is changed in the ancestor 
"""
"""
    before 4/30/2023
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
            this establishes a time that should actually be the tmrca for sequences sampled at the end 

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
from copy import copy
import itertools


def codonInfo():

    stopCodons = ['TAG', 'TAA', 'TGA']
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
    codonlist = []

    revCodons  = {}
    for aa in aalist:
        for ci,cd in enumerate(codons[aa]):
            codonlist.append(cd)
            revCodons[cd] = aa
    
    return codons, aalist, codonlist, revCodons,stopCodons


def makeAncestor(args):
    """
        Try to make an ancestral sequence that is at the equilibrium for preferred and unpreferred codons 
        1.take allDNA and take a random assortment of codons based on appearing in allDNA
        2. replace each codon with a random draw from the the most fit neutral set for that same amino acid
        3. replace favored codons randomly by making random draws fro mthe less fit neutral set at the approximate equilibrium rate 
    """
    
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
        given random gene selected, returns the alignment data for that gene
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
        global stopCodons
        species, gene = readInGeneFile(alignment_location,gene = gene)
        allDNA = ''
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

    global codondic, aalist, codonlist, revcodondic,stopCodons, SynNeuX,SynSelX

    dnaStrand,genefilename = createCodonSequence(args.bacaligndir,gene=args.genename)# if args.genename is None,  then a random gene is picked 
    props = getCodonProportions(dnaStrand)
    ancestor = np.random.choice(list(props.keys()), size = args.aalength, p = list(props.values()))
    bestcodonancestor = []
    fitness = 1.0 
    for codon in ancestor:
        aa = revcodondic[codon]
        if aa in args.neutralsets:
            newcodon = np.random.choice(list(args.neutralsets[aa][0]),replace = True) # random codon in first neutral set
            bestcodonancestor.append(newcodon)
        else:
            bestcodonancestor.append(codon)

    expterm = np.exp(2*args.SynSel_s) # crude attempt at equilibrium setting of some synonymous codons to less fit values 
    expprobs = [0.0,1/(1+expterm),1/(1+(3/4)*expterm),1/(1+(1/3)*expterm)]

    ancestor = []
    nc = 0 
    for codon in bestcodonancestor:
        aa = revcodondic[codon]
        if aa in args.neutralsets:
            k = len(args.neutralsets[aa])
            if aa in args.neutralsets and np.random.random() < expprobs[k-1]: # change to a less favored codon
                ri = np.random.randint(k-1)
                newset = args.neutralsets[aa][1:][ri]
                newcodon = np.random.choice(list(newset))
                assert revcodondic[newcodon]==aa
                fitness *= args.fitnessstructure[codon][newcodon] 
                ancestor.append(newcodon)
                nc += 1
            else:
                ancestor.append(codon)
        else:
            ancestor.append(codon)
    return ''.join(ancestor),fitness,genefilename


def createSelectedDictionary(args):
    """

    """
    def readModelFile(mssmodelfilename):
        """
        format:  Any line beginning with a '#' is a comment
        each line specifies a neutral set of two or more codons for an amino acid 
        the line begins with either the one letter code for an amino acid
        after that comes a series of one or more codons that are respectively neutral to each other (if there is just one, then it is neutral only with respect to itself)
        If all codons for an amino acid are neutral,  then that amino acid has only one line. 
        If an amino acid has multipe sets of neutral codons,  then each set is given on a separate line. 
        The first line for an amino acid gives the neutral set that has the highest fitness
        The second line for an amino acid gives the neutral set that has the next lower fitness
        The third line is lower still.   etc 
        Example:
            # A,E,F,P, S, Y all neutral
            # G   GGG selected with respect to others
            # I ATT selected with respect to others
            # L  4-fold set neutral, 2 fold set neutral,  between them is selected 
            # R  CGT, CGC :neutral  AND  "CGA", "CGG", "AGA", "AGG": neutral,  between them selected
            # T  ACT, ACC :neutral    ACA, ACG : neutral     between them selected 
            # V  GTT, GTC: neutral    GTA, GTG: neutral  between them selected 
            I ATT
            I ATC ATA
            L CTT CTA CTG CTC 
            L TTG TTA
            V GTT GTC 
            V GTA GTG
            F TTT TTC
            C TGT 
            C TGC
            A GCG GCT GCA GCC
            G GGA GGC GGT 
            G GGG
            P CCA CCC CCT CCG
            T ACC ACT 
            T ACG ACA
            S AGT TCG TCT AGC TCA TCC
            Y TAC TAT
            Q CAA 
            Q CAG
            N AAT 
            N AAC
            H CAT 
            H CAC
            E GAA GAG
            D GAT 
            D GAC
            K AAA 
            K AAG
            R CGC CGT 
            R CGA AGA AGG CGG        
        """
        global codonlist
        checkcodons = {codon:False for codon in codonlist if codon not in ["TAA", "TAG", "TGA","ATG","TGG"]} #avoid STOP, Met, Trp
        temp = open(mssmodelfilename,'r').readlines()
        lines = [line for line in temp if line[0] != '#' and len(line) > 1]
        d = {}
        for line in lines:
            ls = line.strip().split()
            aa = ls[0]
            neutralcodonlist = ls[1:]
            for codon in neutralcodonlist:
                assert checkcodons[codon] == False,"codon '{}' specified more than once in MSS model file".format(codon)
                checkcodons[codon] = True
            if aa in d:
                d[aa].append(ls[1:])
            else:
                d[aa] = [ls[1:]]
        assert len(d)==18,"not all degenerate amino acids specified in MSS model file"
        assert all(codonpresent for codonpresent in checkcodons.values()), "missing codons {} in MSS model file".format([key for key, value in checkcodons.items() if value ==False])
        return d


    global codondic,codonlist,revcodondic,stopCodons
    selectedDict = {}
    mutDict = {}
    neutralsetDict = readModelFile(args.mssmodelfilename)
    
    #initialize all with nonsense and missense 
    selectedDict = {}
    mutDict = {}
    for c1,codon1 in enumerate(codonlist):
        tempaaDict = {}
        tempaaMuts = {}
        aa1 = revcodondic[codon1]
        synCodons = codondic[aa1]            
        for c2,codon2 in enumerate(codonlist):
            aa2 = revcodondic[codon2]
            if codon1 in stopCodons or codon2 in stopCodons:
                tempaaDict[codon2] = 0.0  ## stop codon
                tempaaMuts[codon2] = STOPX ## stop codon 
            elif codon2 == codon1:
                tempaaDict[codon2] = 1.0  ## same codon exactly
                tempaaMuts[codon2] = -10 ## same codon  
            elif aa1 != aa2:
                tempaaDict[codon2] = args.NonSyn_s_rescaled  # nonsynonmous
                tempaaMuts[codon2] = NonSynDelX  # [X, -, - ]   
        selectedDict[codon1] = tempaaDict
        mutDict[codon1] = tempaaMuts           
    #now go through and set all synonymous changes 
    for aa in neutralsetDict:
        for ni in range(len(neutralsetDict[aa])):
            codonpairs = list(itertools.combinations(neutralsetDict[aa][ni], 2))
            for (c1,c2) in codonpairs:
                selectedDict[c1][c2] = selectedDict[c2][c1]  = 1.0  # synonymous but neutral
                mutDict[c1][c2] = mutDict[c2][c1] = SynNeuX 
            if len(neutralsetDict[aa]) > 1:
                for nj in range(ni+1,len(neutralsetDict[aa])):
                    codonpairs = list(itertools.product(neutralsetDict[aa][ni],neutralsetDict[aa][nj]))
                    for (c1,c2) in codonpairs:
                        mutDict[c1][c2] = mutDict[c2][c1] = SynSelX         
                        if len(neutralsetDict[aa])==2: # 2 sets,  a full selective difference apart from each other
                            selectedDict[c1][c2] = args.SynSelDel_s_rescaled
                            selectedDict[c2][c1] = args.SynSelFav_s_rescaled
                        elif len(neutralsetDict[aa])== 3: # 3 sets, 2 steps is a full difference,  sets either 1 or 2 steps apart
                            if nj-ni==2:
                                selectedDict[c1][c2] = args.SynSelDel_s_rescaled
                                selectedDict[c2][c1] = args.SynSelFav_s_rescaled
                            else:
                                assert nj-ni==1
                                selectedDict[c1][c2] = args.SynSelDel_half_s_rescaled
                                selectedDict[c2][c1] = args.SynSelFav_half_s_rescaled
                        elif len(neutralsetDict[aa])== 4: # 4 sets, 3 steps is a full difference,  sets either 1, 2 or 3 steps apart
                            if nj-ni == 3:
                                selectedDict[c1][c2] = args.SynSelDel_s_rescaled
                                selectedDict[c2][c1] = args.SynSelFav_s_rescaled
                            elif nj - ni ==2:
                                selectedDict[c1][c2] = args.SynSelDel_twothird_s_rescaled
                                selectedDict[c2][c1] = args.SynSelFav_twothird_s_rescaled
                            else:
                                assert nj - ni == 1
                                selectedDict[c1][c2] = args.SynSelDel_third_s_rescaled
                                selectedDict[c2][c1] = args.SynSelFav_third_s_rescaled
                        else:
                            print("something wrong with neutral sets for {}".format(aa))
                            exit()

    if args.debug and args.debug.usesubtimeinfo:
        codon_substitution_time_info_dict = create_codon_pair_substitution_analysis_structure(selectedDict,mutDict)
    else:
        codon_substitution_time_info_dict = None
    return selectedDict, mutDict, neutralsetDict, codon_substitution_time_info_dict

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
    elif numSpecies == 4:
        tree = '(((p1,p4),p3),p2);'
        split_generations = {1: ['p1', 1, None],
                        2: ['p2', 0.0, 'p1'],  #0.3
                        3: ['p3', 0.4, 'p1'],
                        4: ['p4', 0.8, 'p1']}
        mean_branches_root_to_tip = 2.25
    else: #numSpecies == 1
        tree = '(p1);'
        split_generations = {1: ['p1', 1, None]}
        mean_branches_root_to_tip = 1

    return tree,split_generations,mean_branches_root_to_tip 


def create_codon_pair_substitution_analysis_structure(selectiondic,mtypedic):
    """
        only used when debugging, requires debug.usesubtimeinfo
        make a dict of a dict that will contain generation times for different synonymous substitutions 
    """
    global codondic,revcodondic,SynSelX,  SynNeuX
    cccs = {}

    for codon1 in selectiondic:
        for codon2 in selectiondic[codon1]:
            if codon1 != codon2  and revcodondic[codon1] == revcodondic[codon2] and mtypedic[codon1][codon2] != STOPX:
                aa = revcodondic[codon1]
                if mtypedic[codon1][codon2] ==SynNeuX: # NEUTRAL
                    status = "NEUTRAL"
                    holdcodon2 = codon2
                    holdcodon1 = codon1 
                    if holdcodon1 > holdcodon2:# use alphabetic order, if necessary swap them so codon1 is the alphabetically first 
                        temp = holdcodon2
                        holdcodon2 = holdcodon1
                        holdcodon1 = temp 
                    if holdcodon1 in cccs:
                        cccs[holdcodon1][holdcodon2] = [aa,status,0,0,0,0] #status, number of forward subs,  total time of forward subs,  number of backward subs, total time of backward subs
                    else:
                        cccs[holdcodon1] = {holdcodon2:[aa,status,0,0,0,0]}
                else: #SELECTED,  get index positions of both codons,  the one to the left of the other in codondic is the favored one,  put that on the left in this dictionary
                    assert  mtypedic[codon1][codon2] == SynSelX
                    status = "SELECTED"
                    index1 = codondic[aa].index(codon1)
                    index2 = codondic[aa].index(codon2)
                    holdcodon2 = codon2
                    holdcodon1 = codon1 
                    if index1 > index2: # use index order,  if necessary swap so codon1 is to address cccs directly 
                        temp = holdcodon2
                        holdcodon2 = holdcodon1
                        holdcodon1 = temp 
                    if holdcodon1 in cccs:
                        cccs[holdcodon1][holdcodon2] = [aa,status,0,0,0,0] #status, number of forward (unfavored) subs,  total time of forward subs,  number of backward (disfavored) subs, total time of backward subs
                    else:
                        cccs[holdcodon1] = {holdcodon2:[aa,status,0,0,0,0]}
    return cccs


def makefastafile(samples, fn):
    f = open(fn,'w')
    for pop in samples.keys():
        f.write('>{}\n'.format(pop))
        f.write(str(samples[pop]) + "\n")
    f.close()
    return

class debuginfo():
    """
        holds some booleans for switching on/off different types of debugging
    """
    def __init__(self,useturnoff,usesubtimeinfo,userescaler):
        self.useturnoff = useturnoff
        self.usesubtimeinfo = usesubtimeinfo
        self.userescaler = userescaler

class chromosome():
    """
    a chromosome is an 'individual' in the population
    it has a DNA sequence from which its fitness is determined 
    """

    def __init__(self,sequence,fitness,args,mcounts,ancestornumber,fits = None,subtimeinfo=None):
        """
        mcounts :  nummutationtypes positions 
            0,1,2,3 or 4  for nonsynonymous deleterious, nonsynonymous favored,  synonymous-selected, synonymous-neutral, and STOP 
        re ancestornumber:
          is just the index of the chromosome in the list
          it is reset at the beginning of burn2
          during burn2 ancestornumber is updated for each chromosome, by copying the value from the ancestor the chromosome was copied from
          when all values of ancestornumber are the same,  all the chromosomes are descended from that ancestor in that generation
        re subtimeinfo
            a list of values about the mutations that accumulate on a chromosome
            for each mutation there are 5 values, in order:
                codon position
                generation number
                mutation type
                old codon
                new codon
            these are all appended as they occur,  so the mutations are ordered in time,  and the full list is a length that is a multiple of 5
        """
        global nummutationtypes
        self.s = sequence
        self.fitstruct = args.fitnessstructure
        self.mutstruct = args.mutstructure
        self.mrate = args.mutrate
        self.mrateinverse = 1.0/self.mrate  # use with exponential random variables for mutation locations because numpy.random.exponential() takes inverse of the rate as the parameter for some reason
        self.fitness = fitness
        self.mcounts = copy(mcounts)
        self.ancestornumber = ancestornumber
        self.SynSelFav_s_rescaled = args.SynSelFav_s_rescaled
        if args.debug and args.debug.usesubtimeinfo:
            if subtimeinfo == None:
                self.subtimeinfo = []
            else:
                self.subtimeinfo = copy(subtimeinfo)
            if fits is None:
                self.fits = np.empty(shape=0,dtype=np.float128)
            else:
                self.fits = np.copy(fits)                
        else:
            self.fits = None
    
    def chromosomecopy(self,args):
        if args.debug and args.debug.usesubtimeinfo:
            newchromosome = chromosome(self.s,self.fitness,args,self.mcounts,self.ancestornumber,fits=self.fits, subtimeinfo = self.subtimeinfo)
        else:
            newchromosome = chromosome(self.s,self.fitness,args,self.mcounts,self.ancestornumber)
        return newchromosome

    def mutate(self,popancestor,gennum,adaptivechanges,turnoffmut,debug):
        """
            a function that changes the sequence s and recalculates fitness
            it uses exponential to get the distance to the next mutation as an approximation for geometric
            sites are selected by jumping along the chromosome (distance sampled from exponential)
            fitness is not updated until all of the sites in a codon have been changed
                usually each mutation is in its own codon, but sometimes 2 or 3 changes occur in the same codon
                code is a kind of kludgy, in order to keep a list of all changes in the current codon
            
        """
        global mutationlocations # use in debug mode 
        global  NonSynDelX,  NonSynFavX, SynSelX,  SynNeuX, STOPX 
        global revcodondic
        global args
        pos = 0 # a position that mutates  (if not past the end of the sequence)	
        lastcodonpos = -1
        muttype = -1
        oldCodon = ''
        anc = ''
        newCodon = ''
        fitratio = 0.0
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
                holds  = self.s
                codonpos = mutlocs[0] //3
                oldCodon = self.getOldCodon(mutlocs[0]) 
                # assert oldCodon == self.s[3*codonpos:3*codonpos+3]
                # assert (len(mutlocs)==1 or (len(mutlocs)==2 and (mutlocs[1] -mutlocs[0] <= 2)  )  or (len(mutlocs)==3 and (mutlocs[2] -mutlocs[0] <= 2) and (mutlocs[0] < mutlocs[1] < mutlocs[2]) ))
                for ml in mutlocs:
                    bps =['A', 'G', 'C', 'T']
                    bps.remove(self.s[ml:ml+1])
                    self.s = self.s[:ml] + np.random.choice(bps) + self.s[ml+1:]
                ## update fitness
                oldfitness = self.fitness
                anc,newCodon,muttype = self.fitnessfunction(mutlocs[0], oldCodon,popancestor)
                if muttype == NonSynDelX:
                    if codonpos in adaptivechanges and revcodondic[newCodon] == adaptivechanges[codonpos]:
                        muttype = NonSynFavX
                self.mcounts[muttype] += 1
                if debug:
                    if (turnoffmut == 1 and muttype in (NonSynDelX,NonSynFavX)) or (turnoffmut==2 and muttype in (SynNeuX,SynSelX)):
                        self.s = holds
                    else:
                        for ml in mutlocs:
                            mutationlocations[ml] += 1
                        if gennum >= 0 and debug and debug.usesubtimeinfo:
                            if len(self.subtimeinfo) >= 7:
                                assert gennum >= self.subtimeinfo[-6], "subtimeinfo genum error {}".format(self.subtimeinfo)
                            fitratio = 0.0 if oldfitness <= 0 else self.fitness/oldfitness
                            self.subtimeinfo.append(codonpos)
                            self.subtimeinfo.append(gennum)
                            self.subtimeinfo.append(muttype)
                            self.subtimeinfo.append(fitratio)
                            self.subtimeinfo.append(anc)
                            self.subtimeinfo.append(oldCodon)
                            self.subtimeinfo.append(newCodon)
                            if oldfitness > 0.0: # if it is zero it wont end up getting sampled 
                                self.fits = np.append(self.fits,fitratio)
                mainmutationcounter[muttype] += 1
                if pos > len(self.s):
                    break
                else: # reset mutlocs to contain only the last exponential jump position (before it was incremented)
                    mutlocs = [pos-1] # the next mutation location (pos was previously incremented to the start of the next interval)
            else:
                break
        return
        
    def resetmcounts(self):
        global nummutationtypes
        self.mcounts = [0 for i in range(nummutationtypes)]
        self.subtimeinfo = []

    def fitnessfunction(self, mut,oldcodon,popancestor):
        """
        recalculates fitness based on newCodon and ancestral codon just for mutations
        the ancestor was defined as having a fitness of 1
        for x  codons  fitness is a product of x values 
        the ancestor was defined as having a fitness of 1
        for a change from oldcodon to newcodon the fitness is updated by multiplying 
        fitness times the fitness associated with a change from ancestor to the new codon (i.e. the absolute fitness of the new codon at that position)
        and by dividing by the fitness associate with a change from the ancestor to the old codon 
        """
        global SynSelX, NonSynDelX,NonSynFavX,STOPX 
        global stopCodons,revcodondic
        anc, newSelf = self.findCodon(mut,popancestor)
        # assert newSelf != oldcodon
        
        muttype = self.mutstruct[oldcodon][newSelf]
        if muttype == SynSelX:
            self.fitness *=  self.fitstruct[oldcodon][newSelf]
        elif muttype == NonSynDelX:
            if revcodondic[anc] == revcodondic[newSelf]:
                # assert revcodondic[anc] != revcodondic[oldcodon]
                self.fitness *= 1.0 / self.fitstruct[anc][oldcodon]
            elif revcodondic[anc] == revcodondic[oldcodon]:
                # assert revcodondic[anc] != revcodondic[newSelf]
                self.fitness *= self.fitstruct[anc][newSelf]
            else:
                revcodondic[anc] not in (revcodondic[oldcodon],revcodondic[newSelf])
                self.fitness *= self.fitstruct[anc][newSelf] / self.fitstruct[anc][oldcodon]
        elif muttype == STOPX:
            self.fitness  = 0.0
        return anc,newSelf, muttype

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
            print('getOldCodon error position {} in sequence {}'.format(i,self.s))
            exit()

    def findCodon(self, i,popancestor):
        """
        identify codon that mutation is in for ancestral and sequence
        """
        position = i % 3

        if position == 0:
            return popancestor[i:i+3], self.s[i:i+3]
        elif position == 1:
            return popancestor[i-1:i+2], self.s[i-1:i+2]
        elif position == 2:
            return popancestor[i-2:i+1], self.s[i-2:i+1]
        else:
            print('findcodon() error position {} in ancestor {}'.format(i,popancestor))
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
        self.popancestor = self.args.ancestor
        if isinstance(source,population):
            """
            copy each chromosome from source and put it in the new population
            """
            for chrom in source:
                # self.append(chrom)
                self.append(chrom.chromosomecopy(args))
        else:# at the beginning, fill up the pouplation with chromosomes made from the ancestor 
            for i in range(self.popsize2):
                self.append(chromosome(source,args.ancestorfitness,args, [0 for j in range(nummutationtypes)],i))
        self.adaptivechanges = {}

    def generation(self,gennum):
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

        # this is not faster than above 3 lines of code 
        # unique_vals, inverse_indices,counts = np.unique(fits,return_counts=True, return_inverse=True)
        # unique_indices = np.split(np.argsort(inverse_indices), np.cumsum(np.unique(inverse_indices, return_counts=True)[1])[:-1])
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
                child = self[i].chromosomecopy(self.args)
                child.mutate(self.popancestor,gennum,self.adaptivechanges,self.args.turnoffmut,self.args.debug)
                newpop.append(child)
                
        self.clear()
        for child in newpop:
            self.append(child)
        return numfits
    
    def changeancestor(self):
        global codondic,revcodondic,aalist
        pos = np.random.randint(self.args.aalength)
        codon = self.popancestor[3*pos:3*pos+3]
        aa = revcodondic[codon]
        while True:
            newaa = np.random.choice(aalist)
            if newaa not in (aa,'STOP'):
                break
        newcodon = np.random.choice(codondic[newaa])
        self.popancestor = self.popancestor[:3*pos] + newcodon + self.popancestor[3*(pos+1):]
        assert len(self.popancestor) == 3*self.args.aalength
        for c in self:
            if revcodondic[c.s[3*pos:3*pos+3]] != newaa:
                c.fitness *= self.args.NonSyn_s_rescaled
        self.adaptivechanges[pos] = newaa
        return pos,newaa
        
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

    def __init__(self,args):
        self.treestr = args.tree
        self.args = args
        self.split_generations, self.times = self.translateGens(args.split_generations)
        self.pop0 = population('p1', args.ancestor,args)
        self.pops = {}

    def translateGens(self, sg):
        newSG = {}
        times = []
        for key in sorted(sg.keys()):
            newTime = round(sg[key][1] * self.args.treeDepth )
            newSG[newTime] = [sg[key][0], sg[key][2]]
            times.append(newTime)
        return newSG, times


    def samplefrompops(self):
        """
        samples sequences at the end of the run
        """
        global stopCodons
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
        meanfits = []
        popkeys = self.pops.keys()
        for pop in popkeys:
            meanfits.append(sum([c.fitness for c in self.pops[pop]])/len(self.pops[pop]))
        return meanfits,sum(meanfits)/len(meanfits)

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
            meanfit += temp
            fitlist.append(temp)
            mcountlist.append(self.pops[pop][num].mcounts)
        return meanfit/len(popkeys),fitlist,mcountlist
    
    def checkcodonsubstitutions(self,sampledsequences,cccs):
            global nummutationtypes
            for p in sampledsequences:
                c = sampledsequences[p]
                d = {}
                prevgens = 0
                for i in range(0, len(c.subtimeinfo), 7):
                    pos = c.subtimeinfo[i]
                    sub = [c.subtimeinfo[i], c.subtimeinfo[i + 5], c.subtimeinfo[i + 6]]  
                    gens = c.subtimeinfo[i+1]
                    assert gens >= 0, "{} {}".format(gens, c.subtimeinfo[i:i+7])
                    muttype = c.subtimeinfo[i+2]
                    codon1 = sub[1]
                    codon2 = sub[2]         
                    k = i
                    while True:
                        k -= 5
                        if k >= 0:
                            if c.subtimeinfo[k] == pos:
                                prevcodon1 = c.subtimeinfo[k+5]
                                prevcodon2 = c.subtimeinfo[k+6]
                                prevgens = c.subtimeinfo[k+1]
                                useit = True
                                break
                        else:
                            useit = False
                            break
                    if useit and  codon1 == prevcodon2 and codon2 == prevcodon1:
                        if codon1 in cccs and codon2 in cccs[codon1]:
                            cccs[codon1][codon2][2] += 1
                            cccs[codon1][codon2][3] += gens - prevgens
                        elif codon2 in cccs and codon1 in cccs[codon2]:
                            cccs[codon2][codon1][4] += 1
                            cccs[codon2][codon1][5] += gens - prevgens
                        if sub[0] in d:
                            d[sub[0]].append(codon1)
                            d[sub[0]].append(codon2)
                        else:
                            d[sub[0]] = [codon1,codon2]
                for pos in d: #need to debug why this crashes when this is on when debugging 
                    if len(d[pos])>3:
                        for k in range(1, len(d[pos]), 2):
                            if k <  len(d[pos])-1:
                                assert d[pos][k] == d[pos][k+1],"{} {} {} {}".format(i,pos,d[pos][k],d[pos][k+1] )
            return cccs
    
    def makesubtimeinfoliststring(self,sampledsequences):
        global revcodondic,codondic
        s = ["\npop#\tAApos\tgennum\ttype\tfitratio\tancAA\tancCodon\tfromAA\tfromCodon\tfromIndex\ttoAA\ttoCodon\ttoIndex\n"] 
        for pi,p in enumerate(sampledsequences):
            c = sampledsequences[p]
            i = 0
            subtimeinfo = [c.subtimeinfo[i:i + 7] for i in range(0, len(c.subtimeinfo), 7)]
            for sti in subtimeinfo:
                s.append("{}\t{}\t{}\t{}\t{:.5f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(pi,sti[0],sti[1],sti[2],sti[3],revcodondic[sti[4]],sti[4],revcodondic[sti[5]],sti[5],codondic[revcodondic[sti[5]]].index(sti[5]),revcodondic[sti[6]],sti[6],codondic[revcodondic[sti[6]]].index(sti[6])))
        return ''.join(s)
    
    def calc_fitness_from_changes(self,sampledsequences):
        fitness_from_changes = []
        numchanges = [len(sampledsequences[p].fits)-1 for p in sampledsequences]
        for p in sampledsequences:
            c = sampledsequences[p]
            fitness_from_changes.append(np.prod(c.fits))
        return fitness_from_changes,numchanges

    def makefitchangeinfostring(self,sampledsequences):
        global revcodondic,codondic
        s = ["\npop#\tsub#\tfitchange\tAApos\tgen\ttype\tfromAA\tfromCodon\tfromIndex\ttoAA\ttoCodon\ttoIndex\n"] 
        for pi,p in enumerate(sampledsequences):
            c = sampledsequences[p]
            i = 0
            subtimeinfo = [c.subtimeinfo[i:i + 5] for i in range(0, len(c.subtimeinfo), 5)]
            for f,m in zip(c.fits[1:],subtimeinfo):
                if True:#f> 1:
                    s.append("{}\t{}\t{:.5f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(pi,i,f,m[0],m[1],m[2],revcodondic[m[3]],m[3],codondic[revcodondic[m[3]]].index(m[3]),revcodondic[m[4]],m[4],codondic[revcodondic[m[4]]].index(m[4])))
                i += 1
        return ''.join(s)
    
    def makeancestoraminoacidstring(self):
        global revcodondic,aalist
        sum = 0
        aacountdic = {aa:0 for aa in aalist}
        codons = [self.args.ancestor[i:i+3] for i in range(0,len(self.args.ancestor),3)]
        for codon in codons:
            aacountdic[revcodondic[codon]] += 1
            sum += 1
        ancestoraminoacidstring  = ["\nAncestor Amino Acid Counts and Frequencies:\n"]
        for aa in aacountdic:
            ancestoraminoacidstring.append ("{}\t{}\t{:.3f}\n".format(aa,aacountdic[aa],aacountdic[aa]/sum))
        ancestoraminoacidstring.append("\n")
        return ''.join(ancestoraminoacidstring)

        
    def write_debug_information(self,sampledsequences,rf):
        rf.write("\n\nDEBUGGING STUFF\n\n")

        fitness_from_changes,numchanges = self.calc_fitness_from_changes(sampledsequences)
        rf.write("\nSample chromosome fitnesses calculated from changes\n")
        rf.write("Pop\t#changes\tfitness\n")
        for i in range(len(numchanges)):
            rf.write("{}\t{}\t{:.10f}\n".format(i,numchanges,fitness_from_changes[i]))
        self.args.codon_substitution_time_info_dict = self.checkcodonsubstitutions(sampledsequences,self.args.codon_substitution_time_info_dict)
        rf.write("\tMean Synonymous Substitution Times: NEUTRAL Pairs\n\nAA\tCodon1\tCodon2\t#Frwd\tMeanFrwdGens\t#Bkwd\tMeanBkwdGens\n")
        fc=bc=fsum = bsum = 0
        for codon1 in self.args.codon_substitution_time_info_dict:
            for codon2 in self.args.codon_substitution_time_info_dict[codon1]:
                x = self.args.codon_substitution_time_info_dict[codon1][codon2] # fewer characters
                if x[1] == "NEUTRAL":
                    fc += x[2]
                    fsum += x[3]
                    bc += x[4]
                    bsum += x[5]
        rf.write("Overall\t#Frwd\tMeanFrwdGens\t#Bkwd\tMeanBkwdGens\n")
        rf.write("mean   \t{}\t{:.0f}\t{}\t{:.0f}\n".format(fc,np.nan if fc==0 else fsum/fc,bc,np.nan if bc == 0 else bsum/bc))
        for codon1 in self.args.codon_substitution_time_info_dict:
            for codon2 in self.args.codon_substitution_time_info_dict[codon1]:
                x = self.args.codon_substitution_time_info_dict[codon1][codon2] # fewer characters
                if x[1] == "NEUTRAL":
                    rf.write("{}\t{}\t{}\t{}\t{:.0f}\t{}\t{:.0f}\n".format(x[0],codon1,codon2,x[2],np.nan if x[2]==0 else x[3]/x[2],x[4],np.nan if x[4]==0 else x[5]/x[4]))
        rf.write("\nMean Synonymous Substitution Times: SELECTED Pairs\n\nAA\tCodon1\tCodon2\t#Frwd\tMeanFrwdGens\t#Bkwd\tMeanBkwdGens\n")
        fc=bc=fsum = bsum = 0
        for codon1 in self.args.codon_substitution_time_info_dict:
            for codon2 in self.args.codon_substitution_time_info_dict[codon1]:
                x = self.args.codon_substitution_time_info_dict[codon1][codon2] # fewer characters
                if x[1] == "SELECTED":
                    fc += x[2]
                    fsum += x[3]
                    bc += x[4]
                    bsum += x[5]
        rf.write("Overall\t#Frwd\tMeanFrwdGens\t#Bkwd\tMeanBkwdGens\n")
        rf.write("mean   \t{}\t{:.0f}\t{}\t{:.0f}\n".format(fc,np.nan if fc==0 else fsum/fc,bc,np.nan if bc == 0 else bsum/bc))
        for codon1 in self.args.codon_substitution_time_info_dict:
            for codon2 in self.args.codon_substitution_time_info_dict[codon1]:
                x = self.args.codon_substitution_time_info_dict[codon1][codon2] # fewer characters
                if x[1] == "SELECTED":
                    rf.write("{}\t{}\t{}\t{}\t{:.0f}\t{}\t{:.0f}\n".format(x[0],codon1,codon2,x[2],np.nan if x[2]==0 else x[3]/x[2],x[4],np.nan if x[4]==0 else x[5]/x[4]))
        ancestoraminoacidstring = self.makeancestoraminoacidstring()
        rf.write(ancestoraminoacidstring)
        # changeinfostring = self.makefitchangeinfostring(sampledsequences)
        # rf.write(changeinfostring)
        subtimeinfostring = self.makesubtimeinfoliststring(sampledsequences)
        rf.write(subtimeinfostring)
        return rf

    def summarize_results(self,starttime,sampledsequences):

        class   substitution_info():
            def __init__(self,label, count,subs,totalcount,totalnumgen,args,parent=None):
                if parent is not None:
                    self.label = parent + "_" + label
                else:
                    self.label = label
                while len(self.label) < 18:
                    self.label += ' '
                self.count = count 
                self.subs = subs
                self.mutproportion = count/totalcount
                self.ebp = 3*args.aalength*self.mutproportion
                self.mutperebp = np.nan if self.ebp == 0 else self.count/self.ebp
                self.subrate = np.nan if self.ebp == 0 else self.subs/self.ebp
                self.totalnumgen = totalnumgen
                self.args = args

            def tablestr(self,tabletype,withheader):
                if withheader == False:
                    header = ""
                else:
                    if tabletype == "Mutation":
                        header = "\nMutation Total Counts/Rates\n\tmutation_type    total_count  effective_#bp  proportions	mutations_per_effective_bp:\n"
                    else: # tabletype == "Substitution"
                        header = "\nSubstitution Counts/Rates\n\tsubstitution_type	count_per_gene	per_effective_bp	per_effective_bp_per_branch	per_effective_bp_per_generation\n"
                if tabletype == "Mutation":
                    return "{}\t{}\t{:>12d}\t{:>3.1f}\t{:>6.3g}\t{:.3g}\n".format(header,self.label,self.count,self.ebp,self.mutproportion,self.mutperebp)
                else:
                    assert tabletype == "Substitution"
                    return "{}\t{}\t{:>4.1f}\t{:>4.3g}\t{:>3.3g}\t{:.3g}\n".format(header,self.label,self.subs,self.subrate,self.subrate/self.args.mean_branches_root_to_tip ,self.subrate/self.totalnumgen)

            def ratiostr(self,denominator,label,withheader=False):
                if withheader == False:
                    header = ""
                else:
                    header = "\nRate Ratios:\n"
                while len(label) < 65:
                    label += ' '
                ratio = np.nan if denominator.subrate == 0 else self.subrate/denominator.subrate 
                return "{}\t{}\t{:.3g}\n".format(header,label,ratio)
                        
        global mainmutationcounter
        global nummutationtypes
        global NonSynDelX, NonSynFavX, SynSelX, SynNeuX, STOPX
        mnames = ["NonSyn_Del","NonSyn_Fav","Synon_Sel","Synon_Neu","Stop"]
        meanfit,fitlist,mcountlist = self.fitmutsummary()
        rf = open(self.args.resultsfilename,'w')
        rf.write("mss_sim\n\narguments:\n")
        for arg in vars(self.args):
            if arg=="fitnessstructure" or  arg=="mutstructure" or arg == "codon_substitution_time_info_dict":
                if self.args.debug:
                    rf.write("\t{}: {}\n".format(arg, getattr(self.args, arg)))
                else:
                    rf.write("\t{}: {}\n".format(arg, " - printed only in debug mode")) # quite large 
            else:
                rf.write("\t{}: {}\n".format(arg, getattr(self.args, arg)))

        rf.write("\nFinal Mean Fitness: {:.4g}\n".format(meanfit))
        rf.write("\nMean number of fitness values each generation: {:.1f}\nMean number of individuals per fitness value: {:.1f}\n".format(self.args.meannumfits,self.args.popsize2/self.args.meannumfits))
        rf.write("\nSampled Individual Fitnesses: {}\n".format(fitlist))
        rf.write("\nSampled Individual Mutation Counts ({}): {}\n".format(mnames,mcountlist))  

        totalnumgen = self.args.treeDepth + self.args.burn2_generation_time
        allmuttotsum = sum(mainmutationcounter)
        subsum = [0 for i in range(nummutationtypes)]
        for mc in mcountlist:
            for i in range(nummutationtypes):
                subsum[i] += mc[i]
        for i in range(nummutationtypes): # take the mean count per sampled chromosome
            subsum[i] /= self.args.numSpecies

        #populate substitution info list
        
        # relies on global vars, in order  NonSynDelX, NonSynFavX, SynSelX, SynNeuX, STOPX
        subinfolist = []
        # all nonsynonymous subinfolist[0]
        allmuts = mainmutationcounter[NonSynDelX] + mainmutationcounter[NonSynFavX]
        meansubs = subsum[NonSynDelX] + subsum[NonSynFavX]
        subinfolist.append(substitution_info("NonSyn",allmuts,meansubs,allmuttotsum,totalnumgen,self.args))
        # nonsynonymous deleterious subinfolist[1]
        subinfolist.append(substitution_info("Deleterious",mainmutationcounter[NonSynDelX],subsum[NonSynDelX],allmuttotsum,totalnumgen,self.args,parent="NonSyn"))
        # nonsynonymous favored subinfolist[2]
        subinfolist.append(substitution_info("Favored",mainmutationcounter[NonSynFavX],subsum[NonSynFavX],allmuttotsum,totalnumgen,self.args,parent="NonSyn"))
        # all synonymous subinfolist[3]
        allmuts = mainmutationcounter[SynSelX] + mainmutationcounter[SynNeuX]
        meansubs = subsum[SynSelX] + subsum[SynNeuX]
        subinfolist.append(substitution_info("Synon",allmuts,meansubs,allmuttotsum,totalnumgen,self.args))
        # synonymous selected  subinfolist[4]
        subinfolist.append(substitution_info("Selected",mainmutationcounter[SynSelX],subsum[SynSelX],allmuttotsum,totalnumgen,self.args,parent="Synon"))
        # synonymous neutral subinfolist[5]
        subinfolist.append(substitution_info("Neutral",mainmutationcounter[SynNeuX],subsum[SynNeuX],allmuttotsum,totalnumgen,self.args,parent="Synon"))
        # STOP  subinfolist[6]
        subinfolist.append(substitution_info("STOP",mainmutationcounter[STOPX],subsum[STOPX],allmuttotsum,totalnumgen,self.args))

        #print tables
        tabletype = "Mutation"
        for i,sb in enumerate(subinfolist):
            rf.write(sb.tablestr(tabletype,i==0))
        tabletype = "Substitution"
        for i,sb in enumerate(subinfolist):
            rf.write(sb.tablestr(tabletype,i==0))
        #rate ratios
        rf.write(subinfolist[1].ratiostr(subinfolist[3],"\tdN*/dS (Nonsynonymous_deleterious/Synonymous (selected and neutral)",withheader=True))
        rf.write(subinfolist[0].ratiostr(subinfolist[3],"\tdN/dS (Nonsynonymous_total/Synonymous (selected and neutral)"))
        rf.write(subinfolist[1].ratiostr(subinfolist[5],"\tdN*/dSn (Nonsynonymous_deleterious/Synonymous_Neu)"))
        rf.write(subinfolist[0].ratiostr(subinfolist[5],"\tdN/dSn (Nonsynonymous_total/Synonymous_Neu)"))
        rf.write(subinfolist[4].ratiostr(subinfolist[5],"\tdSs/dSn (Synonymous_Sel/Synonymous_Neu)"))
        if self.args.debug and self.args.debug.usesubtimeinfo:
            rf = self.write_debug_information(sampledsequences,rf)
                
        totaltime = time.time()-starttime
        rf.write("\ntotal time: {}\n".format(time.strftime("%H:%M:%S",time.gmtime(totaltime))))            

        rf.close()
        

    def run_burn1(self):
        """
        run for a multiplier,  burn1X, x popsize2 generations
        200x  seems good enough,  with the 5/5/ version of makeAncestor()
        """
        defaultgens = 200
        if self.args.debug and self.args.debug.userescaler:
                burn1X = int(defaultgens * self.args.treerescaler) 
        else:
               burn1X = defaultgens
        for i in range(burn1X*self.args.popsize2):
            self.pop0.generation(-1)
            if self.args.debug and  i % 1000 == 0:
                    meanfit = sum([c.fitness for c in self.pop0])/len(self.pop0)
                    print(i, meanfit)                
        meanfit = sum([c.fitness for c in self.pop0])/len(self.pop0)
        self.pop0.reset_mutation_counts()  
        if self.args.debug and self.args.debug.usesubtimeinfo:
            for c in self.pop0:
                c.fits = [c.fitness]
        
        return i+1,meanfit
    
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
            self.pop0.generation(gen)
            gen += 1        
        return gen

    def run(self):
        """
        calls run_burn1(), run_burn2() and then  runs for treedepth generations
        """
        self.args.burn1_generation_time,self.args.burn1_mean_fitness = self.run_burn1()
        self.args.burn2_generation_time = self.run_burn2()
        self.args.mutationexpectation_adjusted_for_burn2 = (self.args.mutationexpectation * self.args.treeDepth)/(self.args.treeDepth + self.args.burn2_generation_time)
        allgen = self.args.burn2_generation_time
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
                if np.random.random() < self.args.adaptchangerate:
                    aapos,newaa = self.pops[key].changeancestor()
                numdifferentfitnessvalues = self.pops[key].generation(allgen)
                self.args.meannumfits += numdifferentfitnessvalues
                countpopgens +=1
            
            gen += 1
            allgen += 1
            if self.args.debug and gen % self.args.popsize2 == 0:
                meanfits,meanmeanfit = self.fitCheck()
                print("generation {} ({:.1f}%)  # populations: {} mean fitnesses: {} overall mean: {:.4f}  sample mutation counts: {}".format(gen,100*gen/self.args.treeDepth,len(self.pops.keys()),meanfits,meanmeanfit,self.pop0[0].mcounts))
                
        self.args.meannumfits /= countpopgens
        sample = self.samplefrompops()
        return sample

def parseargs():
    parser = argparse.ArgumentParser("python makeSlimScript_JH.py",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-A", help="path for folder containing bacterial alignment files",
                        dest="bacaligndir",required=True,type = str)
    parser.add_argument("-b", help="base filename for output files (name only, no directories)",
                        dest="basename",required=True,type = str)
    parser.add_argument("-d", help="Debug mode. Generates screen output during run, with mutation counts by position at the end",
                        dest="debug", default=False,action="store_true")
    parser.add_argument("-e", help="random number seed for simulation (and for picking alignment if needed)",
                        dest="ranseed",required=True,type=int)
    parser.add_argument("-F", help="directory path for output fasta file (default is same as for results file)",
                        dest="fdir",type = str)
    parser.add_argument("-g", help="bacterial gene name, optional - if not used a random gene is selected",
                        dest="genename",type=str)
    parser.add_argument("-k", help="Number of species (1, 4, 5 or 11)",
                        dest="numSpecies",required=True,type=int)
    parser.add_argument("-L", help="Length of sequence (# amino acids)",
                        dest="aalength",required=True,type=int)
    parser.add_argument("-m", help="Model file path",
                        dest="mssmodelfilename",required = True,type = str)
    parser.add_argument("-N", help="Population size (diploid)",
                        dest="popsize",required=True,type=int)
    parser.add_argument("-o", help="turn off a type of mutation: 1 for nonsynonymous, 2 for synonymous (requires -d))",
                        dest="turnoffmut",default = 0, type=int)
    parser.add_argument("-R", help="directory path for results file",
                        dest="rdir",default = ".",type = str)
    parser.add_argument("-q", help="compress/expand run time, useful when debugging to get a short run,requires -d",
                        dest="treerescaler",default = 1.0, type=float)
    parser.add_argument("-s", help="Synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))",
                        dest="SynSel_s",default=2,type=float)
    parser.add_argument("-y", help="Non-synonymous population selection coefficient, 2Ns (Slim uses 1-(2Ns/2N))",
                        dest="NonSyn_s",default=10,type=float)
    parser.add_argument("-u", help="expected number of neutral mutations per site, from base of tree",
                        dest="mutationexpectation",default=0.5,type=float)
    parser.add_argument("-w", help="rate per generation of adaptive amino acid change",
                        dest="adaptchangerate",default=0.0,type=float)
    parser.add_argument("-x", help="save detailed substitution info (requires -d, slows things down a good bit)",
                        dest="savesubtimeinfo",default = False,action="store_true")
    return parser

def main(argv):
    global args

    starttime = time.time()
    parser = parseargs()
    if argv[-1] =='':
        argv = argv[0:-1]
    args = parser.parse_args(argv)
    if args.ranseed != None:
        np.random.seed(args.ranseed)
    if args.numSpecies not in [1, 4,5,11]:
        print ("error: -p (# of species) must be 1, 4,5 or 11")
        exit()
    args.commandstring = " ".join(argv)
    if args.debug:
        args.debug = debuginfo(args.turnoffmut != 0,args.savesubtimeinfo == True,args.treerescaler != 1.0)
    else: # args.debug==False
        if args.turnoffmut != 0:
            args.commandstring += "\n\t -o ignored "
            print("-o ignored")
            args.turnoffmut = 0
        if args.treerescaler != 1.0:
            print("-q ignored")
            args.commandstring += "\n\t -q ignored "
            args.treerescaler = 1.0
        if args.savesubtimeinfo:
            print("-x ignored")
            args.commandstring += "\n\t -x ignored "
            args.savesubtimeinfo = False 
    
    args.meannumfits = 0
    args.popsize2 = args.popsize*2
    args.defaulttreeDepth = 100000 # fixed at a specific value # previously scaled by population size  args.treeDepth * args.popsize
    args.mutrate = args.mutationexpectation/args.defaulttreeDepth  # got rid of using theta 4Nu,  as not really relevant here 
    args.treeDepth = round(args.defaulttreeDepth * args.treerescaler)
    if args.debug and args.debug.userescaler:
        args.treeDepth = round(args.defaulttreeDepth * args.treerescaler)
    else:
        args.treeDepth = args.defaulttreeDepth
    #rescale the selection coefficients from 2Ns values to Slim values
    args.SynSelDel_s_rescaled = max(0.0,1.0 - (args.SynSel_s/(args.popsize2)))
    args.SynSelDel_half_s_rescaled = np.power(args.SynSelDel_s_rescaled,1/2)
    args.SynSelDel_third_s_rescaled = np.power(args.SynSelDel_s_rescaled,1/3)
    args.SynSelDel_twothird_s_rescaled = np.power(args.SynSelDel_s_rescaled,2/3)
   
    args.SynSelFav_s_rescaled = 1.0/args.SynSelDel_s_rescaled  # makes more sense to use the reciprocal for a favored change
    args.SynSelFav_half_s_rescaled = 1.0/args.SynSelDel_half_s_rescaled
    args.SynSelFav_third_s_rescaled = 1.0/args.SynSelDel_third_s_rescaled
    args.SynSelFav_twothird_s_rescaled = 1.0/args.SynSelDel_twothird_s_rescaled
    
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
        
    #create global vars
    global mainmutationcounter #when debug update this when mutating
    global nummutationtypes
    global  NonSynDelX,  NonSynFavX, SynSelX,  SynNeuX, STOPX # all refer to positions in the mutation counter arrays,  all end in 'X' 
    
    nummutationtypes = 5
    #positions in mutation counter lists
    NonSynDelX = 0 # nonsynonymous deleterious
    NonSynFavX = 1 # nonsynonymous favored (0's unless adaptive changes can occur via -w  )
    SynSelX = 2  # synonymous selected 
    SynNeuX = 3  # synonymous neutral
    STOPX = 4  # nonsense 
    global codondic, aalist, codonlist, revcodondic,stopCodons
    codondic, aalist, codonlist, revcodondic,stopCodons = codonInfo()
    mainmutationcounter = [0 for i in range(nummutationtypes)]  #positions 0,1,2 or 3 for nonsynonymous,  synonymous-selected, synonymous-neutral, and STOP 
    global mutationlocations # used when debugging for checking distribution of mutation locations
    mutationlocations = [0 for i in range(3*args.aalength)]

    # get ancestral sequence
    
    args.fitnessstructure, args.mutstructure,args.neutralsets,args.codon_substitution_time_info_dict = createSelectedDictionary(args)
    args.ancestor,args.ancestorfitness,genefilename = makeAncestor(args)
    args.genename = genefilename[:genefilename.find('_')]

    args.resultsfilename = op.join(args.rdir,args.basename +  "_" + args.genename + '_results.txt')
    if os.path.exists(args.resultsfilename):
        args.resultsfilename = '{}(1)_results.txt'.format(args.resultsfilename[:-12])
    args.fastafilename = op.join(args.fdir, args.basename +  "_" + args.genename +".fa")
    if os.path.exists(args.fastafilename):
        args.fastafilename = '{}(1).fa'.format(args.fastafilename[:-3])

    # set tree shape
    args.tree, args.split_generations,args.mean_branches_root_to_tip = maketreeshape(args.numSpecies)

    # run the simulation
    sim = tree(args)
    sampledsequences = sim.run()
    sim.summarize_results(starttime,sampledsequences)
    if args.debug:
        print("mutation counts by base position\n",mutationlocations)
    makefastafile(sampledsequences, args.fastafilename)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        main(['-h'])
    else:
        main(sys.argv[1:])
    