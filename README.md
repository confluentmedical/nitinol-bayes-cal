[![License](https://img.shields.io/badge/License-Apache%202.0-yellowgreen.svg)](https://opensource.org/licenses/Apache-2.0)  

# A Probabilistic Calibration Framework for the Abaqus Superelastic Material Model


## Introduction

![Probabilistic Framework Abstract](/assets/graphical_abstract.jpeg)

We implement a Bayesian Inference approach to calibrate the material parameters of a model for the superelastic deformation of NiTi shape memory alloy. We specify a diamond-shaped specimen geometry that is suited to calibrate both tensile and compressive material parameters from a single test. We adopt the Bayesian Inference calibration scheme to take full-field strain measurements obtained using digital image correlation together with global load data as an input for calibration. The calibration itself is performed by comparing the experimentally measured quantities of interest -- strain data at select locations and global load -- with the corresponding results from a simulation library. We present a machine learning based approach to enrich the simulation library and improve the calibration accuracy. This approach is versatile and can be used to calibrate other models of superelastic deformation from data obtained using various modalities. This probabilistic calibration approach can become an integral part of a framework to assess and communicate the risk-informed credibility of simulations performed in the design of superelastic NiTi articles such as medical devices.

In the specific example presented, following six material parameters in the [superelastic material model](https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-superelasticity.htm) implemented in Abaqus finite element modeling framework are calibrated.
* ![E_\textrm{A}](https://render.githubusercontent.com/render/math?math=E_%5Ctextrm%7BA%7D): Austenite stiffness.
* ![E_\textrm{M}](https://render.githubusercontent.com/render/math?math=E_%5Ctextrm%7BM%7D): Martensite stiffness.
* ![\varepsilon_\textrm{tr}](https://render.githubusercontent.com/render/math?math=%5Cvarepsilon_%5Ctextrm%7Btr%7D): Maximum transformation strain.
* ![\sigma_\textrm{UPS}](https://render.githubusercontent.com/render/math?math=%5Csigma_%5Ctextrm%7BUPS%7D): Upper plateau stress.
* ![\sigma_\textrm{LPS}](https://render.githubusercontent.com/render/math?math=%5Csigma_%5Ctextrm%7BLPS%7D): Lower plateau stress.
* ![\sigma_\textrm{CPS}](https://render.githubusercontent.com/render/math?math=%5Csigma_%5Ctextrm%7BCPS%7D): Compression plateau stress.

Matlab code to execute the full calibration pipeline is presented here. For details see the research publication describing the method :point_down:.

### Publication

Harshad M. Paranjape, Kenneth I. Aycock, Craig Bonsignore, Jason D. Weaver, Brent A. Craven, Thomas W. Duerig. "A Probabilistic Approach with Built-in Uncertainty Quantification for the Calibration of a Superelastic Constitutive Model from Full-field Strain Data." Under review.

## Prerequisites

Below are pre-requisite tools to perform the material parameter calibration described here.
* Matlab with Image Processing Toolbox and Statistics and Machine Learning Toolbox.
* Matlab source code in this repository.
* [GWMCMC sampler](https://github.com/grinsted/gwmcmc) by Aslak Grinsted. Add this to the Matlab path.
* A diamond specimen with the geometry specified in the engineering drawing `data/NDC-56-03018_rev3_25102018_00.pdf`. This specimen should be laser-cut and heat treated from the material to which the model is to be calibrated.
* Instron or similar load frame to perform tension experiments on the specimen.
* A 2D digital image correlation (DIC) setup including tools for speckle pattern application on the diamond, digital camera, lighting, and other supporting equipment.
* [NCORR](https://ncorr.com/) DIC data processing software for Matlab.

## Performing the Calibration

### Step 1: Process Experimental Data

![Diamond with speckle pattern](/assets/fig_diamond.jpeg)

The purpose of this step is to obtain full-field surface strain data during a tensile test. Experimental quantities of interest (QoI) will be extracted from this data and calibration will be subsequently performed against these QoI.

Cover a diamond specimen with [speckle pattern](https://www.sciencedirect.com/science/article/pii/S135964621930733X) as shown in the figure. Mount the speckled specimen on an Instron load frame at appropriate temperature and perform a tension test with two steps: Start from zero and load to 0.8 mm, then unload to 0.4 mm. Take images with a DIC setup. Process the DIC images using NCORR. Process images corresponding to those displacement increments at which you will be obtaining simulation output in Step 2 below. Save the results as `data/QoI_expt_DIC_strains.mat`. Copy the Instron load-displacement file to `QoI_expt_load.csv`.

Set `process_mode = 1;` in `NitinolBayesCal.m`. Specify the input values in the first `process_mode` block. Comments related to the inputs are provided in the program file.

### Step 2: Create a Library of Simulations

The goal of this step is to create a library of simulations with material property inputs spanning a reasonable parameter space.

Set `process_mode = 2;` in `NitinolBayesCal.m` and run the script. It will create a number of Abaqus simulation input files in directory `private/SimLib` using the template `private/SimLib_ftgdia-0p8_0p4.template`. The simulation is with the identical specimen geometry and boundary conditions as in the experiment.

Copy the simulation input files to a suitable location and run them using Abaqus. Post-process each resultant `odb` file using two post-processing scripts in the `SimLib` folder: `fem-get_s_e.py` and `fem-get_rf_u.py`. This should create several output files with the extension `data` for each simulation. Copy the `data` files to `SimLib` directory. There is no need to copy the `odb` files.

### Step 3: Batch-process Simulation Results

In this step, the script will batch-process all simulation outputs (i.e., the `data` files) and store simulated QoIs. These will be later compared with the experimental QoI to perform the material property calibration.

Set `process_mode = 3;` in `NitinolBayesCal.m` and run the script.

### Step 4: Perform Least-squares Calibration

![Error in least-squares calibration](/assets/fig_diamond_error.jpeg)

This step in the analysis will provide the results of a least-squares calibration. An objective function is defined that depends on the experimental QoI and simulated QoI from each simulation in the library. The material parameters in the simulation library that minimize the objective function are reported.

Set `process_mode = 4;` in `NitinolBayesCal.m` and run the script. The calibration results will be printed and a plot of objective function vs. each of the material parameter value as shown above will be printed.

### Step 5: Perform Calibration using Bayesian Inference

### Step 6: Fit a Surrogate Model to Simulation QoI using Machine Learning

### Step 7: Perform Calibration using the Surrogate Model and Bayesian Inference

## License

    Copyright 2020 Harshad M. Paranjape

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
