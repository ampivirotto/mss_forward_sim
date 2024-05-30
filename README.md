# mss_forward_sim
Forward Simulator 

Running mss_sim.py:

	python mss_sim.py -A <location> -b <label> -N <number> -L <number> -u <number> -s <number> -y <number> -k <number> -m <model file> -R <location> -e <number>

	where:
	-A is the location of the alignments folder
	-b is the label for the simulation
	-N is the population size 
	-L is the number of amino acids
	-u is the expected number of mutations per site
	-s is the synonymous selection coefficient 
	-y is the nonsynonymous selection coefficient
	-k is the number of species 
	-m is the path to the model file (model files can be found in the github repo)
	-R is the path to results folder 
	-e is the random seed 


Generating multiple simulation commands:

	python make_file_of_mss_sim_commands.py <PARAMETERS>
	python make_file_of_varying_mss_sim_commands.py <PARAMETERS> 

Running multiple simulations simultaneously: 

	python run_many_mss_sim.py <FILE LOCATION> <NUMBER to run at once>

	
Add newick trees and anonymize FASTA files:

	python anon_run_raxml_on_mss_sim_fasta_files.py <FILE FOLDER INPUT> <FILE FOLDER OUTPUT>
	

	