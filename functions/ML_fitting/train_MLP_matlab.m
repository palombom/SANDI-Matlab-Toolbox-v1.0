function [Mdl, training_performances] = train_MLP_matlab(database_train, params_train, n_layers, n_neurons)

% Train a Multi Layer Perceptron regressor for SANDI fitting
%
% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

tic

Mdl = cell(size(params_train,2),1);
training_performances = cell(size(params_train,2),1);

net_structure = zeros(1,n_layers);

for i=1:n_layers
    net_structure(i) = n_neurons;
end

rng(1);

for i = 1:size(params_train,2)
    Mdl{i} = feedforwardnet(net_structure);
    Mdl{i}.trainParam.showWindow = false;
    Mdl{i}.trainParam.showCommandLine = false;
end

parfor i=1:size(params_train,2)
    
    [Mdl{i}, training_performances{i}] = train(Mdl{i}, database_train', params_train(:,i)');

end

tt = toc;

fprintf('DONE - RF trained in %3.0f sec.\n', tt)

end