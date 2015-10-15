Content:
========
This folder and its subfolders contain the Matlab routines constituting MEMO and all models and data 
needed to reproduce the analysis shown in the manuscript and supplement of:
'MEMO - Multi-experiment mixture model analysis of censored data' 



- ./auxiliary: MEMOs heart and soul. The fundamental Matlab functions constituting MEMOs functionality. 
               This folder furthermore contains PESTO, the parameter estimation toolbox developed by Jan Hasenauer (jan.hasenauer@helmholtz-muenchen.de)

- ./basic examples: two examples illustrating 1.) parameter estimation and 2.) model selection with MEMO

- ./SAC_data_analysis: routines to reproduce all analysis of SAC data analysis shown in the manuscript and Supplementary Material

- ./NGF-ERk_data_analysis: routines to reproduce analysis of NGF-Erk data shown in Supplementary Material. Application of MEMO considering a mechanistic model

- ./illustrating _examples_censoring: analysis of simulated data with regard to the effects of censoring as shown in manuscript Figure 1 + 2



Run install_toolboxes.m to add the required folders to your matlab path.

             
Requirements / System:
======================
The analysis has been performed on a Windows computer using MATLAB R2014a.
