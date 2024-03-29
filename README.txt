MetaboPAC: Metabolomics Prediction of Absolute  Concentrations
Leveraging the mass balances of cellular metabolism to infer absolute concentrations from relative abundance metabolomics data
------------------------------------------------------------------------------------------------------------------------------

Folders:
-------
calcFluxWithKinetics: Files to calculate flux reaction rates based on known kinetic equations.
GenData: Files used to simulate in silico data for the four systems.
miscfunctions: Miscellaneous functions
modelInfo: Contains basic information about each of the four systems.
penalty_functions: Contains penalty functions used in the optimization approach.
plots: Code to generate plot figures found in the manuscript.
results: Results generated for and presented in the manuscript.
solveRF_KEapproach: Files used to identify response factors using the kinetic equations approach.
trueRF: Randomly generated response factors (20 replicate sets for each system) that MetaboPAC tries to identify.

Main File Descriptions:
----------------------
MetaboPAC_*.m: Run MetaboPAC on noiseless data with true response factors sampled from a uniform distribution between 1 and 1000.
MetaboPAC_*_noisy.m: Run MetaboPAC on noise-added data with different sampling frequencies with true response factors sampled from a uniform distribution between 1 and 1000.
MetaboPAC_*_log.m: Run MetaboPAC on noiseless data with true response factors sampled from a log uniform distribution between 1 and 1000.
MetaboPAC_generic.m: A generalized version of MetaboPAC. If known, kinetic equations will need to be hard-coded where specified in this file.

Instructions for generating results in manuscript:
-------------------------------------------------
1) Generate noiseless data using first cell of GenData/driver_genDatasets*.m code files.
2) Generate noise-add data with different sampling frequencies using second cell of GenData/driver_genDatasets*.m code files.
3) Run MetaboPAC_*.m for repetitions [1:20] and percKnownKinetics values of [0 20 40 60 80 100].
4) Run MetaboPAC_*_noisy.m for timepoints [50 15], CoV [5 15], repetitions [1:20], and percKnownKinetics values of [0 20 40 60 80 100].
5) Run MetaboPAC_*_log.m for repetitions [1:20] and percKnownKinetics values of 100.
6) Generate figures using files in plots folder.

Instructions for using generic version of MetaboPAC:
---------------------------------------------------
1) generic_data.mat should include relative abundance data (named concMatrix) and information about when each datapoint was measured (timeVec).
2) Edit modelInfo_generic to include the stoichiometric matrix of the system of interest and any know concentration and flux bounds.
3) Edit calcFluxWithKinetics_generic.m and solveRF_KEapproach_generic.m to include known kinetic equations.
4) Create penalty functions within the penalty_functions folder.
5) Run MetaboPAC_generic_noisy.m