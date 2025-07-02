# SingleJointPointToPoint

This repository provides code and data for the paper "The contribution of muscle parameters to ballistic performance during a goal-directed tasK".

There are four categories of script
1. A script to run a single simulation (`MAIN_simpleAntagonist_trackO1.m`)
2. A script to run the large parameter sweep in the paper (`MAIN_parameterSweep.m` and `MAIN_parameterSweep_fineGrid.m`)
3. A script to perform linear regression on the parameter sweep data (`MAIN_linearRegression`)
4. Various scripts to generate the figures in the paper

Data used in the paper are included. The results from the large parameter sweep varying all 5 interrogated muscle parameters, as well as for a finer parameter sweep across changes in deactivation rate and stiffness, are contained in the directory "TrackO1ParameterSweep/VelActFmaxStiffSweep"

## Dependencies

This code is written in MATLAB and runs on MATLAB 2024b- it may not work properly on earlier versions of MATLAB. The statistics and machine learning toolbox is required for the linear regression script.

The optimisation procedure requires [CasADi](https://web.casadi.org/).

## Running the code

With the MATLAB current folder set to the upper directory of the repository, run `addPaths`.

Then run `MAIN_simpleAntagonist_trackO1.m` to perform a simulation with baseline parameter settings.

To generate figures, first unpack the `.zip` files into directories under the `TrackO1ParameterSweep/VelActFmaxStiffSweep` directory. These directories should keep the same name as their respective zip filename- if not, the figure scripts will not run correctly.

# Single simulation script

`MAIN_simpleAntagonist_trackO1.m` sets up a single optimisation problem using baseline parameters. Each input is described in the script. The user may edit these inputs to explore how the output changes.

# Parameter sweep scripts

`MAIN_parameterSweep.m` initialises a matrix of parameter combinations and runs them sequentially through the optimisation routine. The solution of each combination is stored in individual `.mat` files. If run as-is without unpacking the `.zip` data files, it will generate the same data as summarised in Figure 4- although this will take hours to complete on a good computer! 

If the `.zip` files have already been unpacked in the right place, then this script will simply skip over the simulations that have already been completed.

`MAIN_parameterSweep_fineGrid.m` is the same as above, but the parameter space is limited to only explore the stiffness and deactivation time dimensions (all other parameters are held at baseline. These data are summarised in Figure 5.

# Linear regression script

This script performs linear regression on the large parameter sweep data. It is divided into sections to showcase how removing different terms changes the regression outcomes.

This script requires the Statistics and Machine Learning toolbox.

# Figure-generating scripts

Various scripts are provided which generate the figures seen in the paper. Each begins with the filename `FigureX`, with X being the corresponding figure number in the manuscript. These are provided with the goal of being maximally transparent about the data we present.

# Further questions

Any questions about this repository can be directed to Delyle Polet: delylepolet@gmail.com

