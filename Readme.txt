Readme for ASCA analysis of experimental design LC-MS data
Code for article: Vikram Mitra, Natalia Govorukhina, Gooitzen Zwanenburg, Inge Westra, Age Smilde, Theo Reijmers, Ate G. J. van der Zee, Frank Suits, Rainer Bischoff, Péter Horvatovich, Identification of analytical factors affecting complex proteomics profiles acquired in a factorial design study with ANOVA – simultaneous component analysis.

All source code was developed using Matlab 7.14.0.739 (R2012a), and tested on 64 bit Windows 8. To access the source code, download 'ASCA.zip' folder from the GitHub repository. 
Extract the contents of the 'ASCA.zip' folder into a suitable sub-directory eg: 'C:\Users\Vikram\Documents\ASCA'

::Running the simulation code::
1.To run the simulation navigate using 'cd' command to the unzipped ASCA directory. For help for 'cd' command refer - http://uk.mathworks.com/help/matlab/ref/cd.html to change current working directory.
2.Open the 'SimulationParameterOpt.m' script and press the Run button to run the simulation code or press F5.
3.Select the ASCA directory once the modal dialogue box appears.


::Running the Factorial design data analysis::
1.To run the Factorial design data analysis navigate to the unzipped ASCA directory. For help for 'cd' command refer - http://uk.mathworks.com/help/matlab/ref/cd.html to change current working directory.
2.Open the 'FacDesDataAnalysis.m' script and press the Run button to run the simulation code or press F5.
3.Select the ASCA directory once the modal dialogue box appears.

::TicVisualisation::
1.To produce the TIC figures shown in Figure 4 and attached as supporting information in the manuscript, navigate to the unzipped ASCA directory. For help for 'cd' command refer - http://uk.mathworks.com/help/matlab/ref/cd.html to change current working directory.
2.Open the 'TicVisualisation.m' script and press the Run button to run the simulation code or press F5.
3.Select the ASCA directory once the modal dialogue box appears.

File lists:
Readme.txt						this Readme file
SimulationParameterOpt.m		starting script to perform simulation and ASCA analysis of simulated data
FacDesDataAnalysis.m			starting script to perform ASCA analysis of experimental design data
simulation_.m					Simulation script to run from existing workspace after parameter optimisation
volcano_plot.m					matlab script to select peaks based on Volcano filtering (2 ind t-test and fold ratio)
asca.m							code for ASCA analysis
permutations.m					code to calculate SSQ significance based on permutation test
barwitherr.m					script to include error bars in barplot
rotateXLabels.m					script that include rotated X label in matlab plot
TIC_visualisation.m 			script to visualise EICs
FinalAnnotationVectorOk.mat		data to annotate LC-MS peaks with peptide sequence and proteins
FMatrix.txt						boolean matrix representing the distribution of low and high levels of the 7 analytical factors
Outstem_mzRadius=0.3_TRadius=1_Fraction=0.5_mzStart=100_rtStart=65_mzEnd=1500_rtEnd=135.mpks	quantitative data matrix of experimental design LC-MS files using TAPP pre-processing
peakList.mat					matlab data containing list of the most discriminating ASCA selected peaks in experimental design files
factorLabel.mat 				matlab data containing names of factor label
TicVisualisation				directory containing data for TIC visualisation