# MSWD-PS
This repository includes all the MATLAB scripts used in our work entitled: Robust fMRI Time-Varying Functional Connectivity Analysis Using Multivariate Swarm Decomposition.

To replicate the results of our simulations, one can directly use the script simulation1_eval.m (modify for the rest of simulations accordingly). For the parameter sensitivity analysis, one can use the scripts PSA_simulation1.m. To generate new simulated data and re-run the analyses, one can just delete the files inside the "simulates" and "decomposed" folders and run the scripts simulation1.m (modify for the rest of simulations accordingly) prior to simulation1_eval.m

For the application to real rs-fMRI data from the Human Connectome Project, one can use the script real_data_application.m to decompose the signals and perform the time-varying phase synchronization analysis and real_data_application_states.m to extract the results. HCP data and t-maps used to organize time-series into brain networks are available through the links provided in our paper. Users must modify the paths to these data in order for real_data_application.m to run.

For the application to real rs-fMRI data from the ABIDEI dataset, one can use the script real_data_application_ABIDE.m to decompose the signals and perform the time-varying phase synchronization analysis and real_data_application_ABIDE_states.m to extract the results. To download these data, one can use "cyberduck". Following are the exact instructions:
1) Open Cyberduck
2) Click "Open Connection"
3) Select FTP -> Amazon S3
4) Access Key ID: anonymous
5) Path -> /fcp-indi/data/Projects/ABIDE/Outputs/fmriprep/fmriprep
6) From each subject download the following file: "...space-MNI152NLin2009cAsym_desc-smoothAROMAnonaggr_bold.nii"
In our work we used data from the following subjects:
"0050004"
"0050005"
"0050006"
"0050008"
"0050009"
"0050010"
"0050012"
"0050014"
"0050015"
"0050016"
"0050020"
"0050022"
"0050022"
"0050024"
"0050025"
"0050027"
"0050029"
"0050030"
"0050031"
"0050032"
