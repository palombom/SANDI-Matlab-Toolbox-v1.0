function [database_train, database_train_noisy, params_train] = build_training_set(model, output_folder)

% Builds the dataset for supervised training of the machine learning models
% for SANDI fitting.
%
% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

disp('Building training set...')

tic

%% Build the training set

f = model.function;
protocol = model.protocol;
Nparams = model.Nparams;
SNR = model.SNR;
Nset = model.Nset;

params = zeros(Nset, Nparams);

 for i = 1:2
     params(:,i) = acos( sqrt((model.paramsrange(i,2)-model.paramsrange(i,1)).*rand(Nset,1) + model.paramsrange(i,1)));
 end

for i = 3:Nparams
    params(:,i) = (model.paramsrange(i,2)-model.paramsrange(i,1)).*rand(Nset,1) + model.paramsrange(i,1);
end
database = zeros(size(params,1), length(protocol.G));


parfor i=1:size(params,1)
    
    database(i,:) = f(params(i,:), protocol);
    
end

database_train = database;
params_train = params;

if ~isfield(model,'noise') || strcmp(model.noise,'gaussian')

        % Add experimental Gaussian noise at given SNR for training 

    database_train_noisy = squeeze(database_train + 1./SNR(1).*randn(size(squeeze(database_train))));

    for i= 2 : numel(SNR)
        database_train_noisy = [database_train_noisy; squeeze(database_train + 1./SNR(i).*randn(size(squeeze(database_train))))];
        params_train = [params_train; params_train];
    end

elseif strcmp(model.noise,'rician')
    
% Add experimental Rician noise at given SNR for training
    
database_train_noisy = sqrt(squeeze(database_train + 1./SNR(1).*randn(size(squeeze(database_train)))).^2 + (1./SNR(1).*randn(size(squeeze(database_train)))).^2);

for i= 2 : numel(SNR)
    database_train_noisy = [database_train_noisy; sqrt(squeeze(database_train + 1./SNR(i).*randn(size(squeeze(database_train)))).^2 + (1./SNR(i).*randn(size(squeeze(database_train)))).^2)];
    params_train = [params_train; params_train];
end

end

save([output_folder '/database_training_set'], 'database_train', 'params_train')

tt = toc;

disp(['DONE - Training set built in ' num2str(round(tt)) ' sec.'])

end