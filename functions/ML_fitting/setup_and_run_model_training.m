function [Mdl, train_perf, model] = setup_and_run_model_training(schemefile, SNR, output_folder, Dsoma, Din_UB, Rsoma_UB, De_UB, seed, MLmodel)

% Main script to setup and train the Random
% Forest (RF) or multi-layers perceptron (MLP) regressors used to fit the
% SANDI model
%
% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

addpath(genpath([pwd '/functions']));

rng(seed) % for reproducibility
mkdir(output_folder) % create the output folder

protocol = SchemeToProtocol(schemefile);

%% Load data

% Protocol definition
protocol.roots_sphere = BesselJ_RootsSphere(100); % Calculates the roots for the calculation of the signal for diffusion restricted in spheres under GPD approximation 

%% Build the model and the corresponding training set

% Build the SANDI model to fit (as in the Neuroimage 2020) using the parameter transformation introduced in
% doi:10.1016/j.neuroimage.2011.09.081 to ensure that the three signal fractions
% sum up to unity.

f = @(p,prot) (cos(p(1)).^2).*SynthMeasAstroSticks(p(3).*1e-9,prot)' + ...
              (1 - (cos(p(1)).^2) - ((1-(cos(p(1)).^2)).*cos(p(2)).^2)).*SynthMeasSphere([Dsoma.*1e-9 p(4).*1e-6],prot)' + ...
              ((1-(cos(p(1)).^2)).*cos(p(2)).^2).*SynthMeasBall(p(5).*1e-9,prot)';

model = struct;
model.protocol = protocol;
model.noise = 'gaussian'; %'gaussian' to add Gaussian noise; 'rician' to add Rician noise
model.paramsrange = [0 1; 0 1; Dsoma/12 Din_UB; 1 Rsoma_UB; Dsoma/12 De_UB];
model.Nparams = size(model.paramsrange,1);
model.function = f;
model.SNR = SNR;
model.Nset = 1e5;

save(fullfile(output_folder, 'model.mat'), 'model') % Save the model used

% Build the training set

[database_train, database_train_noisy, params_train] = build_training_set(model, output_folder);

disp('Training will be performed with:')
disp(['   - Dsoma = ' num2str(Dsoma) ' um^2/ms'])
disp(['   - Rsoma = [' num2str(model.paramsrange(4,1)) ', ' num2str(model.paramsrange(4,2)) '] um'])
disp(['   - De = [' num2str(model.paramsrange(5,1)) ', ' num2str(model.paramsrange(5,2)) '] um^2/ms'])
disp(['   - Din = [' num2str(model.paramsrange(3,1)) ', ' num2str(model.paramsrange(3,2)) '] um^2/ms'])

if isempty(MLmodel), MLmodel = 'RF'; end

switch MLmodel
    
    case 'RF'
        %% RF train
        
        disp('Training using a Random Forest regressor implemented in matlab ')
        
        % --- Using Matlab
        n_trees = 200;
        disp(['Training the Random Forest with ' num2str(n_trees) ' trees ...'])
        
        Mdl = train_RF_matlab(database_train_noisy, params_train, n_trees);
        train_perf = cell(5,1);
        
        for i=1:5
        train_perf{i} = oobError(Mdl{i});
        end
        
        % save([output_folder '/trained_RFmodel.mat'], 'Mdl', '-v7.3')
        
    case 'MLP'
        %% MLP train
        
        disp('Training using a MLP regressor implemented in matlab ')
        
        % --- Using Matlab
        n_layers = 3;
        n_neurons = 2*min(size(database_train,1),size(database_train,2));

        disp(['Training the MLP with ' num2str(n_layers) ' hidden layers and ' num2str(n_neurons) ' units per layer ...'])
        
        [Mdl, train_perf] = train_MLP_matlab(database_train_noisy, params_train, n_layers, n_neurons);
        % save([output_folder '/trained_MLPmodel.mat'], 'Mdl', '-v7.3')

    case 'GRNN'
        %% GRNN setup
        
        disp('Training using a Generalized Regression Neural Network implemented in matlab ')
        
        % --- Using Matlab
        
        disp('Training the Generalized Regression Neural Network...')
        Mdl = newgrnn(database_train', params_train', 1./model.SNR);
        train_perf = [];
end

end
