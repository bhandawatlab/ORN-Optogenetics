# ORN-Optogenetics

ORN-Optogenetics contains the MATLAB code that accompanies the paper: "Sensorimotor transformation underlying odor-modulated locomotion in walking Drosophila."


## Getting Started

This code was made on MATLAB 2019b and tested on 2022 versions.
To get started, please download or clone the GitHub repository
```shell
$ git clone https://github.com/bhandawatlab/ORN-Optogenetics.git
```
## Usage on dataset presented in the paper

In MATLAB, navigate to the repository folder and run: `Run_Pipeline`

To toggle between using time averaged analysis or if considering adaptation/habituation, please change the adaptation in the Config file to 1(true)/0(false)

The code creates multipage postscript analysis figure files. To convert postscript to pdf, you can use adobe or [ghostscript](https://www.ghostscript.com/)

## Data location

A copy of the dataset is archived in [TBD](https://www.dropbox.com/s/qjyx5voz82onfic/Data.zip?dl=1)

When running `Run_Pipeline`, the dataset should automatically download from the repository. If not, please download directly from the TBD. Then place the zip file into MATLAB's path.



## General Information

Data folder structure example:
```
 |-- Experiment_20230218-08_36
    |-- Data
	    |-- Calibrations
	    |-- DataGen - Consolidated mat files where each file will contain information from all flies for a given genotype
	    |-- DataRaw - raw tracked files from the CircularArenaTrackingCode
		    |-- Genotype name - each genotype folder contains single fly mat files from the tracking code
	    |-- DataSpike - ORN firing rate data (extimated from linear filters)
	    |-- DataStim - Light intensity data
	    |-- SummationModel - Figure 6 summation model data
	    |-- Training Data - single sensillum recordings data
	|-- DataModel - Mat files (1 for each genotype) containing a fly object
	|-- DataRT - Synthetic fly files. Each genotype will have 2 files. The one ending with _flies contains the fly object.
	|-- Figures
	    |-- GeneralBehavior - General quantifications of behaviour such as probability inside, radial occupancy, etc.
		    |-- XY tracks - positions of each fly by genotype sorted according to attraction index
	    |-- f_dfInfo - general information about setting up the f/df space and other analysis
	    |-- Habituation/Time Averaged - KNN space plots
		    |-- KNN_Absolute_Heatmap - colormap is set to absolute values
		    |-- KNN_Relative_Heatmap - colormap is set to relative to before baseline
		    |-- KNN_Emp_KS_test - 1.) KS test between control and retinal 2.) Permutation test if looking at habituation
		    |-- Inhibition - how kinematics and turn optimality change over time when the firing rate is inhibited
	    |-- LinearFilterAnalysis - figures related to TTA, GLM, and linear filters from firing rate to kinematics/turn probability
	    |-- EmpSynth - Analysis figures for synthetic flies (e.g. KNN space) and for comparisons with empirical flies (e.g. radial occupancy)
	    |-- SummationAnalysis - analysis figures for how ORNs combine to influence behavior.
```
## Authors and Citation

Code by Liangyu Tao [lt532@drexel.edu](mailto:lt532@drexel.edu)

Tao, L., Wechsler, S. P., & Bhandawat, V. (2022). “Sensorimotor transformation underlying odor-modulated locomotion in walking Drosophila”, bioRxiv,  [doi:10.1101/2022.03.15.484478](https://doi.org/10.1101/2022.03.15.484478)
