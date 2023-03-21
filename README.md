# SANDI Matlab Toolbox v1.0

This is the very first version of the SANDI Matlab Toolbox, no longer updated. **The most recent and up-to-date version of the SANDI Matlab Toolbox can be found [here](https://github.com/palombom/SANDI-Matlab-Toolbox-Latest-Release)**. We recommend users to use the most recent version. 

The "***SANDI (Soma And Neurite Density Imaging) Matlab Toolbox***" enables model-based estimation of MR signal fraction of brain cell bodies (of all cell types, from neurons to glia, namely soma) and cell projections (of all cell types, from dentrites and myelinated axons to glial processes, namely neurties ) as well as apparent MR cell body radius and intraneurite and extracellular apparent diffusivities from a suitable diffusion-weighted MRI acquisition using Machine Learning (see the original SANDI paper for model assumptions and acquisition requirements DOI: https://doi.org/10.1016/j.neuroimage.2020.116835).

The current implementation of SANDI Matlab Toolbox allows to choose between Random Forest regressor (option 'RF', default); Multi Layer Perceptron (option 'MLP') or Generalized Regression Neural Network (option 'GNRR') to fit the SANDI model to given diffusion-weighted data.

For queries about the toolbox or suggestions on how to improve it, please email: palombom@cardiff.ac.uk

## Dependencies
To use SANDI Matlab Toolbox you will need a MATLAB distribution with the Parallel Computing Toolbox, the Statistics and Machine Learning Toolbox and the Optimization Toolbox. Additionally, you will also need an external repository that is already included in the SANDI Matlab Toolbox:
* [Tools for NIfTI and ANALYZE image] Jimmy Shen (2022). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange. Retrieved April 16, 2022.

## Download 
To get the SANDI Matlab Toolbox clone this repository. The tools include all the necessary dependencies and should be ready for you to run.

If you use Linux or MacOS:

1. Open a terminal;
2. Navigate to your destination folder;
3. Clone SANDI Matlab Toolbox:
```
$ git clone https://github.com/palombom/SANDI-Matlab-Toolbox-v1.0.git 
```
4. The main script is called "SANDI_batch_analysis" and it explains how to analyse one or more datasets with the SANDI model. 
5. A very simple Matlab App is also available via the installer "SANDI Matlab Toolbox.mlappinstall" for simple single-subject analysis or for quick tests of performances. This is an experimental version, tested on MacOS Catalina (10.15.7) and MATLAB R2019b Update 4. It can potantially be used to process multiple datasets by performing Step 1 and 3 for each dataset and Step 2 only for the first dataset (i.e. the Machine Learning model is trained only once, stored in the Matlab environment, and can be used to analyse any other datasets until the App is closed or the training Step 2 rerun).
6. You should now be able to use the code. 

## Usage
The script "SANDI_batch_analysis" shows how to use the toolbox to perform the SANDI analysis on one or more datasets. To use the App, just open Matlab and run the "SANDI Matlab Toolbox.mlappinstall" installer (this may work with only MATLAB >=R2019b). 

## Citation
If you use SANDI Matlab Toolbox, please remember to cite our main SANDI work:

1. Marco Palombo, Andrada Ianus, Michele Guerreri, Daniel Nunes, Daniel C. Alexander, Noam Shemesh, Hui Zhang, "SANDI: A compartment-based model for non-invasive apparent soma and neurite imaging by diffusion MRI", NeuroImage 2020, 215: 116835, ISSN 1053-8119, DOI: https://doi.org/10.1016/j.neuroimage.2020.116835. 

and our preclinical optimization:

2. Andrada Ianuş, Joana Carvalho, Francisca F. Fernandes, Renata Cruz, Cristina Chavarrias, Marco Palombo, Noam Shemesh, "Soma and Neurite Density MRI (SANDI) of the in-vivo mouse brain and comparison with the Allen Brain Atlas", NeuroImage 2022, 254: 119135, ISSN 1053-8119, DOI: https://doi.org/10.1016/j.neuroimage.2022.119135.


## License
SANDI Matlab Toolbox is distributed under the BSD 2-Clause License (https://github.com/palombom/SANDI-Matlab-Toolbox/blob/main/LICENSE), Copyright (c) 2022 Cardiff University and University College London. All rights reserved.

**The use of SANDI Matlab Toolbox MUST also comply with the individual licenses of all of its dependencies.**

## Acknowledgements
The development of SANDI was supported by EPSRC (EP/G007748, EP/I027084/01, EP/L022680/1, EP/M020533/1, EP/N018702/1, EP/M507970/1) and European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Starting Grant, agreement No. 679058). Dr. Marco Palombo is currently supported by the UKRI Future Leaders Fellowship MR/T020296/2.
