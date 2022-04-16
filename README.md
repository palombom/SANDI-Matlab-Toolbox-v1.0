# SANDI Matlab Toolbox
The "***SANDI (Soma And Neurite Density Imaging) Matlab Toolbox***" enables model-based estimation of MR signal fraction of brain cell bodies (of all cell types, from neurons to glia, namely soma) and cell projections (of all cell types, from dentrites and myelinated axons to glial processes, namely neurties ) as well as apparent MR cell body radius and intraneurite and extracellular apparent diffusivities from a suitable diffusion-weighted MRI acquisition (see the original SANDI paper for model assumptions and acquisition requirements DOI: https://doi.org/10.1016/j.neuroimage.2020.116835).

The current implementation of SANDI Matlab Toolbox uses either Random Forest regressor (option 'RF', default) or Multi Layer Perceptron (option 'MLP') to fit the SANDI model to given diffusion-weighted data.

## Dependencies
To use SANDI Matlab Toolbox you will need a MATLAB distribution with the Parallel Computing Toolbox, the Statistics and Machine Learning Toolbox and the Optimization Toolbox. Additionally, you will also need external repositories that are already included in the SANDI Matlab Toolbox, such as:
* [Tools for NIfTI and ANALYZE image] Jimmy Shen (2022). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange. Retrieved April 16, 2022.

## Download 
To get the SANDI Matlab Toolbox clone this repository. The tools includes all the necessary dependencies and should be ready for you to run.

If you use Linux or MacOS:

1. Open a terminal;
2. Navigate to your destination folder;
3. Clone SANDI Matlab Toolbox:
```
$ git clone https://github.com/palombom/SANDI-Matlab-Toolbox.git 
```
4. SARDU-Net tools (i.e. a python class and command-line tools to train and use objects of that class) are now available in the [`./sardunet/ainas`](https://github.com/fragrussu/sardunet/tree/master/ainas) folder (*ainas* means *tools* in [**Sardinian language**](http://sc.wikipedia.org/wiki/Limba_sarda)), while [`./sardunet/tutorials`](https://github.com/fragrussu/sardunet/tree/master/tutorials) contains a number of tutorials. 
5. You should now be able to use the code. 

## User guide and tutorials
SARDU-Net is implemented in the `sardunet_v1` python class, defined in the [`sardunet.py`](https://github.com/fragrussu/sardunet/blob/master/ainas/sardunet.py) file. Details on `sardunet_v1` methods can be found in the [user guide](https://github.com/fragrussu/sardunet/blob/master/tutorials/README.md). 

Command line tools are also provided to help you train and use `sardunet_v1` objects. These tutorials demonstrate how to use the tools:  

* [**Tutorial 1**](https://github.com/fragrussu/sardunet/tree/master/tutorials/tutorial1.md) shows how to extract voxels for SARDU-Net training; 

* [**Tutorial 2**](https://github.com/fragrussu/sardunet/tree/master/tutorials/tutorial2.md) shows how to train a SARDU-Net, choosing the learning options and accessing qMRI sub-protocols selected by the trained SARDU-Net; 

* [**Tutorial 3**](https://github.com/fragrussu/sardunet/tree/master/tutorials/tutorial3.md) shows how to use a trained SARDU-Net to downsample or upsample a qMRI scan.

## Citation
If you use SANDI Matlab Toolbox, please remember to cite our main SANDI work:

1. Marco Palombo, Andrada Ianus, Michele Guerreri, Daniel Nunes, Daniel C. Alexander, Noam Shemesh, Hui Zhang, "SANDI: A compartment-based model for non-invasive apparent soma and neurite imaging by diffusion MRI", NeuroImage 2020, 215: 116835, ISSN 1053-8119, DOI: https://doi.org/10.1016/j.neuroimage.2020.116835. 

and our preclinical optimization:

2. Andrada Ianuş, Joana Carvalho, Francisca F. Fernandes, Renata Cruz, Cristina Chavarrias, Marco Palombo, Noam Shemesh, "Soma and Neurite Density MRI (SANDI) of the in-vivo mouse brain and comparison with the Allen Brain Atlas", NeuroImage 2022, 254: 119135, ISSN 1053-8119, DOI: https://doi.org/10.1016/j.neuroimage.2022.119135.


## License
SANDI Matlab Toolbox is distributed under the BSD 2-Clause License (https://github.com/palombom/SANDI-Matlab-Toolbox/LICENSE.md), Copyright (c) 2022 Cardiff University. All rights reserved.

**The use of SANDI Matlab Toolbox MUST also comply with the individual licenses of all of its dependencies.**

## Acknowledgements
The development of SANDI was supported by EPSRC (EP/G007748, EP/I027084/01, EP/L022680/1, EP/M020533/1, EP/N018702/1, EP/M507970/1) and European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Starting Grant, agreement No. 679058). Dr. Marco Palombo is currently supported by the UKRI Future Leaders Fellowship MR/T020296/2.
