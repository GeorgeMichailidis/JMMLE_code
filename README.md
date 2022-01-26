## Source code for implementation of the JMMLE algorithm in paper:
Joint Estimation and Inference for Data Integration Problems based on Multiple Multi-layered Gaussian Graphical Models
by Subhabrata Majumdar and George Michailidis

* JMLE.R: contains main functions for implementation: jmmle.1step() and jmmle(). Instructions inside.
* Generator.R: auxiliary functions for generating synthetic datasets.
* Objval.R: auxiliary functions for calculating objective function in each iteration.
* l1LS_Main.R: multitask regression wrapper for grpreg for calls inside each JMMLE teration.
* jsem.R: implementation of JSEM method (Ma, Michailidisims 2016) for lower-layer neighborhood coefficients estimation.
* sim_est_new.R: function to generate synthetic data, run JMMLE and generate performance metrics. Contains an example data setup.
