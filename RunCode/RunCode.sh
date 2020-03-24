#Directories
#/disco1/SIMULACIONES/w+jets/w+jets_460/Events/run_01/m_delphes_events.root
#--------------------------------------------------------------
ANALYZERFOLDER="/home/n.cardonac/AnalysisCode_SUSY_VBF_Higgsino"
PHENOANALYZERFOLDER="/home/n.cardonac/AnalysisCode_SUSY_VBF_Higgsino/PhenoAnalyzer/"
PHENOANALYZEREXE="PhenoAnalyzer"
INROOTFILE="Events/run_01/m_delphes_events.root"
TEMPORALFOLDER="/home/n.cardonac/AnalysisCode_SUSY_VBF_Higgsino/RunCode/outPuts/temporal"
OUTPUTFOLDER="/home/n.cardonac/AnalysisCode_SUSY_VBF_Higgsino/RunCode/outPuts/higgsino_complete"
#--------------------------------------------------------------
# Processes
#---------------------------------------------------------------
PROCESSFOLDER[1]="/disco2/SIMULACIONES/ttbar"
PROCESSSSUBFOLDER[1]="ttbar"
RUNS[1]=500
TIMES[1]=1

PROCESSFOLDER[2]="/disco1/SIMULACIONES/z+jets"
PROCESSSSUBFOLDER[2]="z+jets"
RUNS[2]=450
TIMES[2]=1

PROCESSFOLDER[3]="/disco1/SIMULACIONES/w+jets"
PROCESSSSUBFOLDER[3]="w+jets"
RUNS[3]=550
TIMES[3]=1

PROCESSFOLDER[4]="/disco3/SIMULACIONES/ww"
PROCESSSSUBFOLDER[4]="ww"
RUNS[4]=250
TIMES[4]=1

PROCESSFOLDER[5]="/disco3/SIMULACIONES/zz"
PROCESSSSUBFOLDER[5]="zz"
RUNS[5]=200
TIMES[5]=1

PROCESSFOLDER[6]="/disco3/SIMULACIONES/wz"
PROCESSSSUBFOLDER[6]="wz"
RUNS[6]=200
TIMES[6]=1

PROCESSFOLDER[7]="/disco3/SIMULACIONES/SUSY/VBF/Higgsino/m_n2_200_c1_175_n1_150"
PROCESSSSUBFOLDER[7]="m_n2_200_c1_175_n1_150"
RUNS[7]=2
TIMES[7]=1

PROCESSFOLDER[8]="/disco3/SIMULACIONES/SUSY/VBF/Higgsino/m_n2_400_c1_385_n1_370"
PROCESSSSUBFOLDER[8]="m_n2_400_c1_385_n1_370"
RUNS[8]=2
TIMES[8]=1

PROCESSFOLDER[9]="/disco3/SIMULACIONES/SUSY/VBF/Higgsino/m_n2_100_c1_75_n1_50"
PROCESSSSUBFOLDER[9]="m_n2_100_c1_75_n1_50"
RUNS[9]=2
TIMES[9]=1

PROCESSFOLDER[10]="/disco3/SIMULACIONES/SUSY/VBF/Higgsino/m_n2_100_c1_80_n1_60"
PROCESSSSUBFOLDER[10]="m_n2_100_c1_80_n1_60"
RUNS[10]=2
TIMES[10]=1

#---------------------------------------------------------------

# index process to run
#-----------------------
INDEX[1]=1
INDEX[2]=10
#----------------------

# Cut value
#--------------
VARIABLE="nocut"
#--------------

typeset -i start_cut=1
typeset -i end_cut=1
typeset -i delta=1
#-------------------------
# loop over the cut values
#------------------------
for ((i = $start_cut; i <= $end_cut; i = i + $delta)); do
	parameters="sed s/mm/$i/ $PHENOANALYZERFOLDER/initial_parameters.in"
	$parameters >$PHENOANALYZERFOLDER/config.in
	cd $PHENOANALYZERFOLDER
	make compile_ROOT_Delphes
	cp *.h $ANALYZERFOLDER
	cp *.in $ANALYZERFOLDER
	cd -
	#-----------------------
	#loop over the processes
	#----------------------
	for j in $(seq ${INDEX[1]} ${INDEX[2]}); do
		typeset -i start_times=1
		typeset -i end_times=${TIMES[$j]} #Number of repetitions
		typeset -i start_runs=1
		typeset -i end_runs=${RUNS[$j]} #Number of runs
		#loop over repetitions
		for k in $(seq $start_times $end_times); do
			#loop over runs
			for l in $(seq $start_runs $end_runs); do
				$PHENOANALYZERFOLDER/$PHENOANALYZEREXE ${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_$l.root &
				echo Running ${PROCESSSSUBFOLDER[$j]}_$l.root file
				echo ${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE >> allProcesses.dat
			done
			wait
			start_runs=start_runs+${RUNS[$j]}
			end_runs=end_runs+${RUNS[$j]}
		done
		hadd -f $OUTPUTFOLDER/${PROCESSSSUBFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$VARIABLE-$i.root $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_*.root
		cd $TEMPORALFOLDER
		rm *
		cd -
		echo FINISH ${PROCESSSSUBFOLDER[$j]}
	done
done
