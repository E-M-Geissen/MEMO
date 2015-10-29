Content:
========
This folder and its subfolders contain the Matlab routines constituting MEMO and all models and data 
needed to reproduce the analysis shown in the manuscript and supplement of:
'MEMO - Multi-experiment mixture model analysis of censored data' 



- ./auxiliary: MEMOs heart and soul. The fundamental Matlab functions constituting MEMOs functionality. 
               This folder furthermore contains PESTO, the parameter estimation toolbox developed by Jan Hasenauer (jan.hasenauer@helmholtz-muenchen.de)

- ./basic examples: two examples illustrating 1.) parameter estimation and 2.) model selection with MEMO

- ./SAC_data_analysis: routines to reproduce all analysis of SAC data analysis shown in the manuscript and Supplementary Material

                       Fig.4A: .\analysis\Mad3_datasets\delta\Mad3delta_cont
                       Fig.4B: .\analysis\Mad3_datasets\delta\Mad3delta_dis
                       Fig.4C: .\analysis\wild_type\WT_Mad2_right_cen_disregarded
                       Fig.4D: .\analysis\wild_type\Mad2_WT_MEMO

                       Fig.5BC: .\analysis\Mad2_datasets\log_normal\2_components_1_fixed_to_wt\with_Mc
                       Fig.5B:  .\analysis\Mad2_datasets\log_normal\2_components_free
                       Fig.6:   .\analysis\Mad2_Mad3_double_perturbations


- ./NGF-ERk_data_analysis: routines to reproduce analysis of NGF-Erk data shown in Supplementary Material. Application of MEMO considering a mechanistic model

- ./illustrating _examples_censoring: analysis of simulated data with regard to the effects of censoring as shown in manuscript Figure 1 + 2


Files: 
install_toolboxes.m: Run to add the required folders to your matlab path.
anonymous_run_estimation: Just insert the name of your own model and explore MEMO funtionality.

             
Requirements / System:
======================
The analysis has been performed on a Windows computer using MATLAB R2014a.
