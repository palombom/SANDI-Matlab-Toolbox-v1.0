% Main script to perform the SANDI analysis (on one or more datasets) using machine learning, as in Palombo M. et al. Neuroimage 2020: https://doi.org/10.1016/j.neuroimage.2020.116835

% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

clear all
close all
clc

%% Add the path to main and support functions used for SANDI analysis

addpath(genpath(fullfile(pwd, 'functions')));

%% STEP 1 - Load b values and acquisition parameters to train the Machine Learning (ML) model

%%%%%%%%%%%%%%%%%%%%%%%%%%% USER DEFINED INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bvalues_filename = ; % in s/mm^2, FSL-like format
delta = ; % in ms, the diffusion gradient separation in a PGSE sequence
smalldelta = ; % in ms, the diffusion gradient duration in a PGSE sequence
SNR = ; % Average signal-to-noise ratio (SNR) of an individual b=0 image for the whole brain

% Define the parameters to build the training set. This can (and should) change according to the acquisition!

% UB stands for upper bound. 1E5 model parameters values will be chosen randomly between the intervals:
% fsoma and fneurite = [0 1]
% Din and De = [Dsoma/12 Din_UB] and [Dsoma/12 De_UB],
% Rsoma = [1 Rsoma_UB]

Dsoma =  ; % in micrometers^2/ms, e.g. 2
Din_UB =  ; % in micrometers^2/ms, e.g. 3
Rsoma_UB = ; % in micrometers, e.g. 12
De_UB = ; % in micrometers^2/ms, e.g. 3
seed_rng = 1; % for reproducibility
MLmodel = 'RF'; % can be 'RF' for Random Forest (default); 'MLP' for multi-layer perceptron; 'GRNN' for Generalized Regression NEural Network

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('*****   SANDI analysis using Machine Learning based fitting method   ***** ')

outputfolder = fullfile(pwd, 'Acq_Params'); % Folder where the schemefile containing the acquisition parameters will be saved
mkdir(outputfolder);
load_acquisition_parameters(bvalues_filename, delta, smalldelta, outputfolder); % function that loads and saves the acquisition parameters in the file ~/Acq_Params/diravg.scheme

% Train the ML model
schemefile = fullfile(pwd, 'Acq_Params', 'diravg.scheme');
[Mdl, train_perf] = setup_and_run_model_training(schemefile, SNR, outputfolder, Dsoma, Din_UB, Rsoma_UB, De_UB, seed_rng, MLmodel); % Funciton that train the ML model. For options 'RF' and 'MLP', the training performance (mean squared error as a function of #trees or training epochs, respectively) are also saved in 'train_perf'

%% STEP 2 - PreProcess each subject

% Here each subject can be preprocessed to compute the direction averaged
% signal by inserting the code below within a for loop over each subject

%%%%%%%%%%%%%%%%%%%%%%%%%%% USER DEFINED INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_folder = ; % output folder where the direction-averaged signal and SANDI fit results for each subject will be stored
data_filename = ; % each subject DWI. HIGHLY RECOMMENDED: first process the raw data to correct them from motion and artifacts. Suggested pipeline incluide, in this order: 1 - denoising, 2 - Gibbs ringing correction, 3 - Outliers and Eddy current distortions and motion corrections, 4 - gradient non-linearities correction, 5 - B1 inhomogeneity correction, 6 - Rician bias / noise floor correction.
mask_data = ; % estimated mask of the brain, in nifti format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_direction_average(output_folder, data_filename, bvalues_filename, delta, smalldelta)

%% Step 3 - SANDI fit

% Here each subject can be analysed using the SANDI model by inserting the code below within a for loop over each subject

img_data = fullfile(output_folder, 'diravg_signal.nii.gz'); % data to process: direction-averaged signals for each subject
schemefile = fullfile (output_folder,'diravg.scheme'); % corresponding acquisition scheme
run_model_fitting(img_data, mask_data, schemefile, output_folder, Mdl, MLmodel);