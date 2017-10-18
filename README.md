# FEL_code_PBPL
1d FEL code
PERAVE 1-D FEL code basic documentation

What is Perave?

Perave is a matlab based code package which solves the 1-D FEL equations for a helical undulator system. The motivation for writing this code was to have a fast, simple tool which allowed us to study high efficiency FELs and the physics of undulator tapering. 

How to run the code

First you will need to place each function and script in the same folder. Next, the script you need to run to start the code from the Matlab Command line is “perave_MainCode_7h.m”. You will first need to have specified the input parameters in a separate Matlab script called “Perave_User_Input_7h.m”. The input parameters which you can specify are commented in the sample user input file which is provided.

Breakdown of the code work-flow

Opening perave_MainCode_7h.m shows a basic breakdown of the code which is the following:

1)	The script loads the user initial conditions from “Perave_User_Input_7h.m”
2)	Some FEL parameters are calculated (Pierce parameter, gain length) 
3)	The script computes the undulator field profile as specified by the user in the input conditions. The code currently supports either a constant undulator parameter (no tapering) or a tapered undulator field, whose taper profile is determined by the resonant phase specified by the user in the input file.
4)	The script generates the initial electron beam input distribution
5)	The script performs the integration of the FEL equations
6)	The output is post-processed and some plots of interest (power,bunching etc.) are shown.


Example output

The plots below show typical results for a soft X-ray FEL run (1 nm). The (helical) undulator period is 3 cm, the RMS undulator parameter is 2.475. The e-beam energy is 5.3 GeV, the relative energy spread is 5.7*10-4 and the peak current is 4 kA. The radiation field starts from zero and the amplification is from noise.

The power grows exponentially until the saturation power (Psat~Pbeam) is reached after roughly 20 gain lengths. The gain length from simulation is 51.75 cm, in good agreement with theoretical value of 51.51 cm. The average electron energy loss and energy spread at saturation are approximated given by the Pierce parameter.


