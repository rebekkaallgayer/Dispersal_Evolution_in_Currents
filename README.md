# Dispersal_Evolution_in_Currents
This repository provides source code needed to run the model which produced the results presented in Dispersal Evolution in Currents by Allgayer et al , submitted to Proceedings of the Royal Society B

Authors: Allgayer, Rebekka L 1,*; Scarpa, Alice 1; Fernandes, Paul G 1; Wright, Peter J 2; Lancaster, Lesley 1; Bocedi, Greta 1; Travis, Justin MJ 1
1: Institute of Biological and Environmental Sciences, University of Aberdeen, Aberdeen, UK; 
2: Marine Scotland Science, Marine Laboratory, Victoria Road, AB11 9DB, Aberdeen, UK; 
*Corresponding Author; Email: r.allgayer.17@abdn.ac.uk


The coding language used is C++, run on the Visual Studio 2017 Compiler
Please note that the code itself requires user-specific editing to set the working directory as well as boolean True/False statements to determine which types of outputs should be produced (ie to track genotype distributions, the "bool geno_dist" should be put to "true" at line 35). Settings as they have been submitted here will produce the baseline results for current strength values 0.5-1 (intervals of 0.1) with reflective boundaries for 100 reps, once the user has inputted working directory details.

Note also that where files are saved to investigate key parameters, filenames would have to be changed before running simulations to ensure files are saved to the correct repository. Not all combinations of parameters are accounted for within this code as is BUT RUNNING THE SIMULATIONS ONLY REQUIRES CHANGING OF FILENAMES AND PARAMETERS BEFORE THE FOR LOOP.

Outputs produced will require the following structure:

in wd//:

- Results
  - change_curr #the baseline values will save in this repository 
  - track_geno #if genotype distribution is outputted, it will be saved here
  - change_mut #if mutation rate is investigated 
  - change_ldd
  
  
