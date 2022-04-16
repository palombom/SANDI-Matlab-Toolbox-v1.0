function [] = run_model_fitting(img_data, mask_data, schemefile, output_folder, Mdl, MLmodel)

% Main script to fit the SANDI model using Random
% Forest (RF) and/or multi-layers perceptron (MLP) 
%
% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

addpath(genpath([pwd '/functions']));

% tmp_model = load([output_folder '/trained_RFmodel.mat']);
% Mdl = tmp_model.Mdl;

%% Test Set

protocol = SchemeToProtocol(schemefile);


%% Load data

if ~isempty(mask_data)
    tmp = load_untouch_nii(mask_data);
    mask = double(tmp.img);
end

tmp = load_untouch_nii(img_data);
nifti_struct = tmp;

I = double(tmp.img);
[sx, sy, sz, vol] = size(I);

disp(['Data ' img_data ' loaded:'])
disp(['  - matrix size = ' num2str(sx) ' x ' num2str(sy) ' x ' num2str(sz)])
disp(['  - volumes = ' num2str(vol) ])

disp('Protocol loaded:')
disp(['  - gradient pulse duration ~ ' num2str(round(protocol.smalldel(1)*1e3)) ' ms'])
disp(['  - gradient pulse separation ~ ' num2str(round(protocol.delta(1)*1e3)) ' ms'])
disp(['  - diffusion time ~ ' num2str(round(protocol.delta(1)*1e3  - protocol.smalldel(1)*1e3/3)) ' ms'])
disp(['  - b values ~ ' num2str(round(GetBvalues(protocol)./1e6)) ' s/mm^2'])
disp(['  - #' num2str(sum(round(GetBvalues(protocol)./1e6)<100)) ' b0s and #' num2str(sum(round(GetBvalues(protocol)./1e6)>=100)) ' b shells'])

if isempty(mask_data), mask = ones(sx,sy,sz); end

% Prepare ROI for fitting
ROI = reshape(I, [sx*sy*sz vol]);
m = reshape(mask, [sx*sy*sz 1]);
signal = (ROI(m==1,:));

% Remove nan or inf and impose that the normalised signal is >= 0
signal(isnan(signal)) = 0; signal(isinf(signal)) = 0; signal(signal<0) = 0;

%% Fitting the model to the data using pretrained models

if isempty(MLmodel), MLmodel='RF'; end

switch MLmodel
    
    case 'RF'
%% RF fit

disp('Fitting using a Random Forest regressor implemented in matlab ')

% --- Using Matlab

disp('Applying the Random Forest...')
mpgMean = apply_RF_matlab(signal, Mdl);

    case 'MLP'
        
%% MLP fit

disp('Fitting using a MLP regressor implemented in matlab ')

% --- Using Matlab

disp('Applying the MLP...')
mpgMean = apply_MLP_matlab(signal, Mdl);   

    case 'GRNN'
        
%% MLP fit

disp('Fitting using a GRNN implemented in matlab ')

% --- Using Matlab

disp('Applying the GRNN...')
mpgMean = apply_GRNN_matlab(signal, Mdl); 

end
%% Calculate and save SANDI parametric maps

names = {'fneurite', 'fsoma', 'Din', 'Rsoma', 'De', 'fextra'};

fneu = cos(mpgMean(:,1)).^2;
fe = (1-cos(mpgMean(:,1)).^2).*cos(mpgMean(:,2)).^2;

fsom = 1 - fneu - fe;

fneurite = fneu ./ (fneu + fsom + fe);
fsoma= fsom ./ (fneu + fsom + fe);
fextra = fe ./ (fneu + fsom + fe);

disp('Saving SANDI parametric maps')

for i=1:size(mpgMean,2)+1
    
    itmp = zeros(sx*sy*sz,1);
    
    if i==1
        itmp(mask==1) = fneurite;
    elseif i==2
        itmp(mask==1) = fsoma;
    elseif i==size(mpgMean,2)+1
        itmp(mask==1) = fextra;
    else
        itmp(mask==1) = mpgMean(:,i);
    end
    
    itmp = reshape(itmp,[sx sy sz]);
    
    nifti_struct.img = itmp;
    nifti_struct.hdr.dime.dim(5) = size(nifti_struct.img,4);
    if size(nifti_struct.img,4)==1
        nifti_struct.hdr.dime.dim(1) = 3;
    else
        nifti_struct.hdr.dime.dim(1) = 4;
    end
    nifti_struct.hdr.dime.datatype = 16;
    nifti_struct.hdr.dime.bitpix = 32;
    
    save_untouch_nii(nifti_struct,[output_folder '/SANDI-fit_' names{i} '.nii.gz']);
    disp(['  - ' output_folder '/SANDI-fit_' names{i} '.nii.gz'])
    
end

end

