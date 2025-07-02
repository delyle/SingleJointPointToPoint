# SingleJointPointToPoint

This repository provides code and data for the paper "The contribution of muscle parameters to ballistic performance
during a goal-directed tasK".

There are four categories of script
1. A script to run a single simulation (`MAIN_simpleAntagonist_trackO1.m`)
2. A script to run the large parameter sweep in the paper (`MAIN_parameterSweep.m` and `MAIN_parameterSweep_fineGrid.m`)
3. A script to perform linear regression on the parameter sweep data (`MAIN_linearRegression`)
4. Various scripts to generate the figures in the paper

Data used in the paper are included. The results from the large parameter sweep varying all 5 interoggated muscle parameters, as well as for a finer parameter sweep across changes in deactivation rate and stiffness, are contained in the directory "TrackO1ParameterSweep/VelActFmaxStiffSweep"

## Dependencies

This code is written in MATLAB and runs on MATLAB 2024b. Earlier versions may not work properly. The statistics and machine learning toolbox is required for the linear regression script.

The optimisation procedure requires [CasADi](https://web.casadi.org/).

## Running the code

With the MATLAB current folder set to the upper directory of the repository, run `addPaths`.

Then run `MAIN_simpleAntagonist_trackO1.m` to perform a simulation with baseline parameter settings.

# Single simulation script

`MAIN_simpleAntagonist_trackO1.m` sets up a single optimisation problem using baseline parameters. Each input is described in the script. The user may edit these inputs to explore how the output changes.

# Parameter sweep scripts

`MAIN_parameterSweep.m` initialises a matrix of parameter combinations and runs them sequentially through the optimisation routine. The solution of each combination is stored in individual `.mat` files.

