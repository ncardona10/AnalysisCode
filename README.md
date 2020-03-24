# AnalysisCode
Root data analyzer to perform phenomenological analysis.  

Project Organization
------------

├── LICENSE.txt                            <- MIT license 
├── PhenoAnalyzer                          <- Analysis code
│   ├── Cuts.h
│   ├── DelphesFunctions.h
│   ├── HepMCFunctions.h
│   ├── LeptonCounter.h
│   ├── Makefile
│   ├── Overlaps.h
│   ├── PhenoAnalyzer
│   ├── PhenoAnalyzer.cc
│   ├── PhenoAnalyzer.h
│   ├── Physics.h
│   ├── ROOTFunctions.h
│   ├── Test.root
│   ├── compile
│   ├── counts.csv
│   ├── countsTemplate.csv
│   ├── run
│   └── run_signal
├── README.md                              <- 
├── RunCode                                <- Run the analysis code in PhenoAnalyzer and storage the outputs
│   ├── DelphesFunctions.h         
│   ├── HepMCFunctions.h
│   ├── PhenoAnalyzer.h
│   ├── ROOTFunctions.h
│   ├── RunCode.sh
│   ├── outPuts
│   └── runCanvas.sh
├── effiencies
│   ├── DelphesFunctions.h
│   ├── Makefile
│   ├── PhenoAnalyzer.h
│   ├── ROOTFunctions.h
│   ├── efficiencies_errors
│   ├── efficiencies_errors.cc
│   ├── error_bayesian.h
│   ├── error_binomial.h
│   ├── errors.csv
│   └── numEvents.csv
├── histograms
├── normalizeLuminosity
│   ├── config_xsec.json
│   ├── joinHistos.cc
│   ├── jsonLibrary.hpp
│   ├── readConfig.cpp
│   └── stackHistos.cc
└── plots
    ├── DelphesFunctions.h
    ├── Makefile
    ├── MyHistograms.h
    ├── PhenoAnalyzer.h
    ├── ROOTFunctions.h
    ├── compileJoin2Histos
    ├── compileJoinBGHistos
    ├── compileJoinHistos
    ├── compileJoinSignalHistos
    ├── compileStackHistos
    ├── join2histos
    ├── join2histos.cc
    ├── joinBGHistos
    ├── joinBGHistos.cc
    ├── joinHistos
    ├── joinHistos.cc
    ├── joinSignalHistos
    ├── joinSignalHistos.cc
    ├── results
    ├── stackHistos
    └── stackHistos.cc
    
------------
