# ks_2Phenotype_phage
keller-segel type of solver (PDE), modified to include predation by 2 phages with different phenotype.
Parameters for cell and phage behavior are in 'parameters.m' 

Simulations with different conditions can be run using one of several wrappers, each titled 'Main_().m'
Each of these wrappers scans over multiple combinations of two parameters at a time and simulates 
two populations of bacteria migrating on a single shared resource pool. 

The actual simulations are performed by the 'simulateWave.m' scripts, which in most cases, should not need to be edited.

Simple analysis of the outputs (comparing fitness of cell populations with different susceptibility to phage) 
can be found in the 'Analysis_()/' directories. 

The two subdirectories contain copies of the same scripts, but the parameter for phage diffusion is 0 or 40, respectively

This work was done in collaboration with Mike Blazanin.
