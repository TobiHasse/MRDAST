# MRDAST
Meandering River Dynamics and Storage Time
This repository was created by Tobias Hasse (tobiack@udel.edu) to host source code related to my thesis research.

### Function and use
- MRDAST uses input files which are output files from LOMR (https://github.com/TobiHasse/LOMR)
- The input files containt a MATLAB struct the river planform through time as well as other parameters for the meandring river model run
- MRDAST outputs include various metrics of the dynamic evolution of the channel planform
- and: the storage time distributions for sediment deposited in the floodplains

### Input and output files
- Input files will be similar to the files from my dissertation on https://doi.org/10.5281/ZENODO.5651841
- output files will be similar to the files from my dissertation on https://doi.org/10.5281/ZENODO.5651874

The code in this repository is required for my dissertation: Hasse, Tobias Raphael. "Storage Time Dynamics of Meandering River Floodplain Sediments: A Modeling Study." PhD diss., University of Delaware, 2021.

The dissertation includes the code as appendices and that code runs to completion but is dependant on two MATLAB toolboxes:
- The Statistics and Machine Learning Toolbox
- The Image Processing Toolbox
- The Fit Toolbox

### This repository contains additional code to eliminate the Statistics and Machine Learning Toolbox dependency:
- pctile_TRH.m is a homemade approximation using the algorithm described in the MATLAB online documentation for pctile.m
- additionally this repository updates ....
- _____ compatibility for MATLAB 2021a

### This repository remains incomplete 
- it is current and completes through to creating most of the meander dynamis analysis figures and the storage time outputs
- the code for creating curve fits to the storage time distributions will be posted later

A package of code files will be added to this repository at a later time.
