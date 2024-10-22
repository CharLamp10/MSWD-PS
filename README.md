# MSWD-PS
This repository includes all the MATLAB scripts used in our work: Robust fMRI Time-Varying Functional Connectivity Analysis Using Multivariate Swarm Decomposition.

To replicate the results of our simulations, one can directly use the scripts simulation1.m and simulation1_eval.m (modify for the rest of simulations accordingly). For the parameter sensitivity analysis, one can use the scripts PSA_simulation1.m

For the application to real rs-fMRI data from the Human Connectome Project, one can use the script real_data_application.m to decompose the signals and perform the time-varying phase synchronization analysis and real_data_application_states.m to extract the results. HCP data and t-maps used to organize time-series into brain networks are available through the links provided in our paper. Users must modify the paths to these data in order for real_data_application.m to run.
