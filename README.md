[![License](https://img.shields.io/badge/License-Apache%202.0-yellowgreen.svg)](https://opensource.org/licenses/Apache-2.0)  

# A Probabilistic Calibration Framework for the Abaqus Superelastic Material Model


## Introduction

![Probabilistic Framework Abstract](/assets/graphical_abstract.png){:width="60%"}

We implement a Bayesian Inference approach to calibrate the material parameters of a model for the superelastic deformation of NiTi shape memory alloy. We specify a diamond-shaped specimen geometry that is suited to calibrate both tensile and compressive material parameters from a single test. We adopt the Bayesian Inference calibration scheme to take full-field strain measurements obtained using digital image correlation together with global load data as an input for calibration. The calibration itself is performed by comparing the experimentally measured quantities of interest -- strain data at select locations and global load -- with the corresponding results from a simulation library. We present a machine learning based approach to enrich the simulation library and improve the calibration accuracy. This approach is versatile and can be used to calibrate other models of superelastic deformation from data obtained using various modalities. This probabilistic calibration approach can become an integral part of a framework to assess and communicate the risk-informed credibility of simulations performed in the design of superelastic NiTi articles such as medical devices.

Matlab code to execute the full calibration pipeline is presented here. For details see the research publication describing the method :point_down:.

### Publication

Harshad M. Paranjape, Kenneth I. Aycock, Craig Bonsignore, Jason D. Weaver, Brent A. Craven, Thomas W. Duerig. "A Probabilistic Approach with Built-in Uncertainty Quantification for the Calibration of a Superelastic Constitutive Model from Full-field Strain Data." Under review.

## Prerequisites

Below are pre-requisite tools to perform the material parameter calibration described here.
* Matlab with Image Processing Toolbox and Statistics and Machine Learning Toolbox.
* A diamond specimen with the geometry specified in the engineering drawing `data/NDC-56-03018_rev3_25102018_00.pdf`. This specimen should be laser-cut and heat treated from the material to which the model is to be calibrated.
* Instron or similar load frame to perform tension experiments on the specimen.
* A 2D digital image correlation (DIC) setup including tools for speckle pattern application on the diamond, digital camera, lighting, and other supporting equipment.
* [NCORR](https://ncorr.com/) DIC data processing software for Matlab.
* Matlab source code in this repository.

## Performing the Calibration

### Step 1: Process Experimental Data

Cover a diamond specimen with [speckle pattern](https://www.sciencedirect.com/science/article/pii/S135964621930733X). Mount the speckled specimen on an Instron load frame at appropriate temperature and perform a tension test with two steps: Start from zero and load to 0.8 mm, then unload to 0.4 mm. Take images with a DIC setup. Process the DIC images using NCORR. Process images corresponding to those displacement increments at which you will be obtaining simulation output in Step 2 below. Save the results as `data/QoI_expt_DIC_strains.mat`. Copy the Instron load-displacement file to `QoI_expt_load.csv`.

Set `process_mode = 1;` in `NitinolBayesCal.m`. Specify the input values in the first `process_mode` block. Comments related to the inputs are provided in the program file.

### Step 2: Create a Library of Simulations

### Step 3: Batch-process Simulation Results

### Step 4: Perform Least-squares Calibration

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
