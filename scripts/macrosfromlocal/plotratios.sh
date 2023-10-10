#!bin/bash
#for PRIOR in "pythia" #"powlaw45" "powlaw5" "powlaw55" "tsalis_1" "tsalis_2" "tsalis_3" "tsalis_4" "tsalis_5" "tsalis_6" "tsalis_7" "tsalis_8" "tsalis_9"
for PRIOR in "Pythia6" #"pythia" "powlaw45" "powlaw5" "powlaw55" "tsalis_1" "tsalis_2" "tsalis_3" "tsalis_4" "tsalis_5" "tsalis_6" "tsalis_7" "tsalis_8" "tsalis_9"
do
	export PRIOR
	for R in 0.2 0.4 #0.3
	do
		export R
		for CENTRALITY in "central" "peripheral" 
		do
			export CENTRALITY
			for UTYPE in "Bayes" "SVD"
			do
				export UTYPE
				for PTLEAD in 0.0 5.0 7.0 #9.0
				do
					export PTLEAD
					root unfoldedratios.cxx+ -b -q
				done #ptlead
			done #unfolding type
		done	#centrality
	done	#R
done #prior
