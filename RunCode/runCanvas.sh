#Directories
#/disco1/SIMULACIONES/w+jets/w+jets_460/Events/run_01/m_delphes_events.root
#--------------------------------------------------------------
ANALYZERFOLDER="/home/n.cardonac/AnalysisCode"
PHENOANALYZERFOLDER="/home/n.cardonac/AnalysisCode/PhenoAnalyzer_SUSY_VBF_Higgsino/"
PHENOANALYZEREXE="PhenoAnalyzer"
INROOTFILE="Events/run_01/m_delphes_events.root"
TEMPORALFOLDER="/home/n.cardonac/RunPhenoCodes/temporal"
OUTPUTFOLDER="/home/n.cardonac/RunPhenoCodes/outputfiles2"
#--------------------------------------------------------------
# Processes
#---------------------------------------------------------------

PROCESSSSUBFOLDER[1 ]="ttbar"
PROCESSSSUBFOLDER[2 ]="z+jets"
PROCESSSSUBFOLDER[3 ]="w+jets"
PROCESSSSUBFOLDER[4 ]="ww"
PROCESSSSUBFOLDER[5 ]="zz"
PROCESSSSUBFOLDER[6 ]="wz"
PROCESSSSUBFOLDER[7 ]="m_n2_200_c1_175_n1_150"
PROCESSSSUBFOLDER[8 ]="m_n2_400_c1_385_n1_370"
PROCESSSSUBFOLDER[9 ]="m_n2_100_c1_75_n1_50"
PROCESSSSUBFOLDER[10]="m_n2_100_c1_80_n1_60"

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

#-----------------------
#loop over the processes
#----------------------
for j in $(seq ${INDEX[1]} ${INDEX[2]}); do

	echo /home/n.cardonac/RunPhenoCodes/outputfiles2${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE
	# $PHENOANALYZERFOLDER/$PHENOANALYZEREXE ${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_$l.root &
	# echo Running ${PROCESSSSUBFOLDER[$j]}_$l.root file
	# echo ${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE >>allProcesses.dat

	echo FINISH ${PROCESSSSUBFOLDER[$j]}
done
