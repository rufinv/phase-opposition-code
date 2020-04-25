Phase Opposition Code  -  R. VanRullen (2016)

This Matlab code is intended to accompany the paper: VanRullen, R. (2016). How to evaluate phase differences between trial groups in ongoing electrophysiological signals. (submitted)
Please cite this paper whenever using this code in a published study.

This code is provided without explicit or implied guarantee, and without any form of technical support. The functions are free to use for any personal, scientific or commercial application.

Download PhaseOppositionCode.zip - a single file containing this entire folder.

Contents:
PhaseOpposition.m        =        a function to compute phase opposition (p-values) between two trial groups.
test_PhaseOpposition.m   =        an example of usage of the PhaseOpposition function. Just run the script.
ArtificialDataset.mat    =        the example dataset used in test_PhaseOpposition.
matrix_circ_wwtest.m     =        a matrix-compatible version of the circ_wwtest function (circstats toolbox)
ITCbaseline.m            =        a function to estimate uniform-distribution ITC as a function of trial number.
combine_pvalues.m        =        a function to combine results of PhaseOpposition across datasets (/observers).
To get started:

After downloading, unzipping and placing the folder in your Matlab path:

-run the script 'test_PhaseOpposition' to view an example. Open the script for details.
-type 'help PhaseOpposition' for details on function usage.
